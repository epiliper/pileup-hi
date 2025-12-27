use crate::alignment::PileupAlignment;
use crate::bamio::OutputDataDest;
use crate::engine::BUFWRITER_CAP;
use crate::utils::{get_writer, temp_fname};
use anyhow::Error;
use log::warn;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::sync::Mutex;

pub static FILE_MERGE_SINGLETON: Mutex<OutputFileMerge> = Mutex::new(OutputFileMerge {
    outfile: OutputDataDest::Stdout,
    subfiles: vec![],
});

/// The interface requirements for a pileup output. It needs to give ref information,
/// intake pileup alignments, update current ref info, display depth, and write itself.
pub trait OrderedPileupOutput: Send + Sync + Clone + std::fmt::Debug {
    /// Get the reference of the pileup
    #[allow(dead_code)]
    fn tid(&self) -> i32;
    /// Get the coordinate of the pileup
    #[allow(dead_code)]
    fn pos(&self) -> i64;
    /// Update internal data with pileup alignment
    fn intake(&mut self, p: &PileupAlignment, refseq: Option<&[u8]>) -> Result<(), Error>;
    /// Update reference data given ref num, pos, name, and sequence
    fn set_ref_info(&mut self, tid: i32, pos: i64, ref_name: &str, ref_seq: Option<&[u8]>);
    fn write<W: std::io::Write>(&mut self, writer: &mut W) -> Result<(), Error>;
    fn depth(&self) -> u32;
    fn clear(&mut self);
    fn new() -> Self;
}

/// Used to keep track of our main output file and the subfiles we want to merge into it. Subfiles
/// are ordered by thread ID.
#[derive(Clone)]
pub struct OutputFileMerge {
    pub outfile: OutputDataDest,
    pub subfiles: Vec<OutputDataDest>,
}

impl OutputFileMerge {
    /// Copy the data in the subfiles over to the main file
    pub fn merge<W: std::io::Write>(&self, mut dest: W) -> Result<(), Error> {
        for s in &self.subfiles {
            match s {
                OutputDataDest::Stdout => anyhow::bail!("cannot merge from stdout! Critical error"),
                OutputDataDest::File(ref f) => {
                    let fhandle = File::open(f)?;
                    let mut reader = BufReader::with_capacity(2 * 1024 * 1024, fhandle);
                    std::io::copy(&mut reader, &mut dest)?;
                    std::fs::remove_file(f)?;
                }
            }
        }

        Ok(())
    }

    /// Delete all files at once, useful if we abruptly exit
    pub fn cleanup(&mut self) -> Result<(), Error> {
        for s in &self.subfiles {
            if let OutputDataDest::File(f) = s {
                if let Err(e) = std::fs::remove_file(f) {
                    match e.kind() {
                        std::io::ErrorKind::NotFound => (),
                        _ => anyhow::bail!(e),
                    }
                }
            }
        }

        Ok(())
    }

    /// If main thread, return the writer to the final output file
    pub fn get_writer(&self, thread_idx: usize) -> Result<TempOutputHandle, Error> {
        if thread_idx == 0 {
            get_writer(&self.outfile, BUFWRITER_CAP, true, true)
        } else {
            get_writer(&self.subfiles[thread_idx - 1], BUFWRITER_CAP, true, false)
        }
    }
}

pub fn generate_subfile_dests(outprefix: &str, n: usize, suffix: &str) -> Vec<OutputDataDest> {
    let mut ret = Vec::with_capacity(n);
    for i in 0..n {
        let temp = temp_fname(outprefix, &i.to_string(), suffix);
        ret.push(OutputDataDest::File(temp));
    }

    ret
}

pub struct TempOutputHandle {
    pub writer: BufWriter<Box<dyn Write>>,
}

impl TempOutputHandle {
    pub fn write(&mut self, data: &[u8]) {
        let _ = self.writer.write_all(data);
    }
}

pub fn setup_exit_handler() {
    ctrlc::set_handler(|| {
        warn!("Received termination signal. Cleaning up intermediate files...");
        if let Ok(mut outputs) = FILE_MERGE_SINGLETON.lock() {
            outputs
                .cleanup()
                .expect("Failed to cleanup temp files during termination sequence");
        }

        std::process::exit(130);
    })
    .expect("Failed to set exit handler")
}

/// A chunked dynamic array used for batching data writes and reducing allocations, intended for
/// multithreading where a worker thread also owns its writer. Chunks span a range of coordinates,
/// each of which should be assigned its output or None if the coordinate failed to meet an output
/// criterion (e.g. depth). Each chunk contains a sub-chunk of size write_batch_size. Once this
/// chunk is filled, it is written to the output. Once all sub-chunks of a chunk have been written,
/// a new chunk is allocated.
pub struct PileupOutputArray<T: OrderedPileupOutput> {
    data: Vec<Vec<Option<T>>>,
    capacity: usize,
    pub cur_entry: usize,
    cur_chunk: usize,
    write_batch_size: usize,
    output: TempOutputHandle,
    #[allow(dead_code)]
    pub id: usize, // keeping this field in case we want to identify logs by thread
    outbuf: Vec<u8>,
}

impl<T: OrderedPileupOutput> PileupOutputArray<T> {
    pub fn alloc_chunk(&mut self) {
        let n_chunks = self.capacity / self.write_batch_size;
        let remainder = self.capacity % self.write_batch_size;

        self.data = Vec::with_capacity(n_chunks);

        for _ in 0..n_chunks - 1 {
            self.data.push(vec![Some(T::new()); self.write_batch_size]);
        }

        let final_size = remainder + self.write_batch_size;

        self.data.push(vec![Some(T::new()); final_size]);

        self.cur_entry = 0;
        self.cur_chunk = 0;
    }

    pub fn new(capacity: usize, write_batch_size: usize, id: usize, output: TempOutputHandle) -> Result<Self, Error> {
        let outbuf = Vec::with_capacity(write_batch_size * size_of::<T>());
        let mut s = Self {
            data: Vec::new(),
            capacity,
            cur_entry: 0,
            cur_chunk: 0,
            output,
            write_batch_size,
            outbuf,
            id,
        };

        s.alloc_chunk();
        Ok(s)
    }

    pub fn get_current_mut(&mut self) -> &mut T {
        self.data[self.cur_chunk][self.cur_entry].as_mut().unwrap()
    }

    pub fn tombstone(&mut self) {
        self.data[self.cur_chunk][self.cur_entry] = None;
        self.advance();
    }

    pub fn advance(&mut self) {
        self.cur_entry += 1;

        // have enough items to write a batch.
        if self.cur_entry >= self.data[self.cur_chunk].len() {
            self.yield_data_chunk();
        }

        // we wrote the last batch of the chunk, so make a new one.
        if self.cur_chunk >= self.data.len() {
            self.alloc_chunk();
        }
    }

    pub fn yield_data_chunk(&mut self) {
        let batch = std::mem::take(&mut self.data[self.cur_chunk]);

        for mut item in batch.into_iter().flatten() {
            let _ = item.write(&mut self.outbuf);
        }

        self.output.write(&self.outbuf);
        self.outbuf.clear();

        self.cur_chunk += 1;
        self.cur_entry = 0;
    }

    pub fn flush(&mut self) {
        let batch = std::mem::take(&mut self.data[self.cur_chunk]);

        for (i, entry) in batch.into_iter().enumerate() {
            if i >= self.cur_entry {
                break;
            }

            if let Some(mut dat) = entry {
                let _ = dat.write(&mut self.outbuf);
            }
        }

        self.output.write(&self.outbuf);
        self.outbuf.clear();

        self.cur_chunk += 1;
        self.cur_entry = 0;
    }
}

/// Defines how to get output data from iterators from a thread. If using a single thread, we can just print directly and not waste memory queueing output.
pub enum OutputMethod<T: OrderedPileupOutput> {
    WriteDirectly(Box<dyn Write>),
    QueueForOutput(PileupOutputArray<T>),
}
