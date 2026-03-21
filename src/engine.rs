use crate::{
    bamio::{BamDataSource, BamReader, OutputDataDest},
    output::{
        generate_subfile_dests, OrderedPileupOutput, OutputFileMerge, OutputMethod, PileupOutputArray,
        FILE_MERGE_SINGLETON,
    },
    params::{InputParams, PileupParams},
    pileup_iterator::PileupIterator,
    position_queue::{create_region_queue, intervals_from_header, GenomeInterval},
    refseq::{RefSeq, RefSeqHandle},
    utils::OutputWriter,
};

use std::collections::VecDeque;
use std::sync::{Arc, Condvar, Mutex};
use std::time::Instant;

const OUTPUT_ARRAY_YIELD_SIZE: i64 = 2000;
pub const BUFWRITER_CAP: usize = 2 * 1024 * 1024;
pub const MIN_BAM_READ_THREADS: usize = 2;

/// The default minimum number of coordinates to give each thread for processing.
/// This basically exists to prevent doing unnecessary work for very small regions.
/// Can be overridden if you need more horsepower for, say, high-depth regions.
pub const MIN_COORDS_PER_THREAD: i64 = 250_000;

use anyhow::Error;
use log::{info, warn};
use std::io::BufWriter;

pub struct ThreadSignal {
    lock: Mutex<usize>,
    cvar: Condvar,
}

impl ThreadSignal {
    pub fn wait_while(&self) {
        std::mem::drop(
            self.cvar
                .wait_while(self.lock.lock().unwrap(), |free| *free == 0)
                .unwrap(),
        )
    }

    pub fn mark_running(&self) {
        *self.lock.lock().unwrap() -= 1;
    }

    pub fn mark_done(&self) {
        *self.lock.lock().unwrap() += 1;
        self.cvar.notify_one();
    }
}

pub struct PileupWorker {
    id: usize,
    jobid: usize,
    handle: Option<std::thread::JoinHandle<()>>,
    notify: Arc<ThreadSignal>,
}

impl PileupWorker {
    pub fn new(id: usize, notify: Arc<ThreadSignal>) -> Self {
        Self {
            id,
            jobid: 0,
            handle: None,
            notify,
        }
    }

    pub fn is_finished(&self) -> bool {
        self.handle
            .as_ref()
            .and_then(|h| if h.is_finished() { Some(()) } else { None })
            .is_some()
    }

    pub fn join(&mut self) {
        assert!(self.is_finished() && self.handle.is_some());
        self.handle.take().unwrap().join().unwrap();
    }

    pub fn run<T>(
        &mut self,
        id: usize,
        params: PileupParams,
        interval: GenomeInterval,
        src: BamDataSource,
        o: T,
        out: OutputWriter,
        refseq: RefSeqHandle,
    ) where
        T: OrderedPileupOutput + 'static,
    {
        self.jobid = id;
        let notify = Arc::clone(&self.notify);

        self.handle = Some(std::thread::spawn(move || {
            notify.mark_running();

            let mut iterator = PileupIterator::new(
                &src,
                refseq,
                &params,
                o,
                OutputMethod::QueueForOutput(PileupOutputArray::new(
                    std::cmp::min((interval.len() / 10).max(1), OUTPUT_ARRAY_YIELD_SIZE) as usize,
                    out,
                )),
            )
            .unwrap();

            iterator.auto_loop2(&interval).unwrap();

            notify.mark_done();
        }));
    }
}

pub struct ThreadPool {
    workers: Vec<PileupWorker>,
    notify: Arc<ThreadSignal>,
}

impl ThreadPool {
    pub fn new(n_threads: usize) -> Self {
        let notify = Arc::new(ThreadSignal {
            lock: Mutex::new(n_threads),
            cvar: Condvar::new(),
        });

        let mut s = Self {
            notify,
            workers: Vec::with_capacity(n_threads),
        };

        (0..n_threads).for_each(|id| s.workers.push(PileupWorker::new(id, Arc::clone(&s.notify))));

        s
    }

    pub fn get_available(&mut self) -> &mut PileupWorker {
        self.notify.wait_while();

        for worker in self.workers.iter_mut() {
            if worker.is_finished() {
                worker.join();
                return worker;
            }
        }
        unreachable!();
    }
}

pub struct PileupEngine<T: OrderedPileupOutput> {
    intervals: Vec<GenomeInterval>,
    plp_params: PileupParams,
    src: BamDataSource,
    output: T,
    dest: OutputDataDest,
    refseq: RefSeq,
}

impl<T: OrderedPileupOutput + 'static> PileupEngine<T> {
    pub fn initialize(in_params: InputParams, plp_params: PileupParams, output: T) -> Result<Self, Error> {
        let src = BamDataSource::from_string(&in_params.file)?;
        let dest = OutputDataDest::from_string(&plp_params.output);

        let tempreader = BamReader::new(&src, 1)?;
        let header = &tempreader.header;

        let intervals = if let Some(region) = in_params.region {
            create_region_queue(&region, header)?
        } else {
            intervals_from_header(header)?
        };

        let refseq = RefSeq::new();

        Ok(Self {
            intervals,
            plp_params,
            src,
            output,
            dest,
            refseq,
        })
    }

    pub fn run(self) -> Result<(), Error> {
        if self.intervals.is_empty() {
            return Ok(());
        }

        // remove old output file if it exists.
        if let OutputDataDest::File(ref f) = self.dest {
            if std::fs::exists(f)? {
                warn!("Output file {} already exists! Overwriting...", f);

                if let Err(e) = std::fs::remove_file(f) {
                    warn!("Failed to remove file {f}; {e}. Output will be appended...");
                };
            }
        }

        if self.src.has_index()? {
            info!("Found index for for input file {}", self.src.fname()?);
        }

        if self.plp_params.threads == 1 {
            self.run_single()
        } else if !self.src.has_index()? {
            warn!(
                "User asked for more than {} threads but file is unindexed. Running in single-thread mode...",
                self.plp_params.threads
            );
            self.run_single()
        } else {
            info!("Running with {} threads...", self.plp_params.threads);
            self.run_multi()
        }
    }

    /// Use a single thread for both processing and writing.
    pub fn run_single(self) -> Result<(), Error> {
        for interval in self.intervals.iter() {
            let lock = Box::new(BufWriter::with_capacity(2 * 1024 * 1024, std::io::stdout().lock()));

            let refseq_handle = self
                .refseq
                .yield_handle(&interval.name, self.plp_params.refseq.as_deref())?;

            let mut iterator = PileupIterator::new(
                &self.src,
                refseq_handle,
                &self.plp_params,
                self.output.clone(),
                OutputMethod::WriteDirectly(self.output.clone(), lock),
            )?;

            iterator.auto_loop2(interval)?;
            self.refseq.decrement_ref_usage(&interval.name);
        }
        Ok(())
    }

    /// Use separate threads for processing and writing. Each processing thread owns its IO readers for input BAM, index, and any other files.
    /// The problem with this: all threads block until the last per-ref thread finishes.
    /// We need this to dynamically determinate how many chunks across ALL threads, dole out chunks
    /// in order to threads, and make sure file manager respects order.
    ///
    /// Thoughts: we need to be finish as many parts of the same reference before moving on to the
    /// next one. if we are processing eukaryotic genomes and have multiple chromosomes in memory, that
    /// will get ugly very fast; better to have one or two in memory at a time.
    pub fn run_multi(self) -> Result<(), Error> {
        let outprefix = self.src.fname()?;

        // create a list of all jobs based on input references and number of available threads.
        // Note that this list is ordered, and is assumed to be so when merging output files.
        let mut chunks = self
            .intervals
            .into_iter()
            .flat_map(|i| {
                if i.len() < self.plp_params.coords_per_thread {
                    i.chunks(self.plp_params.coords_per_thread)
                        .collect::<Vec<GenomeInterval>>()
                } else {
                    i.n_chunks(self.plp_params.threads as i64)
                        .collect::<Vec<GenomeInterval>>()
                }
            })
            .enumerate()
            .collect::<VecDeque<(usize, GenomeInterval)>>();

        let mut output_merge_lock = FILE_MERGE_SINGLETON.lock().expect("Failed to lock output file mutex");
        *output_merge_lock = OutputFileMerge {
            outfile: self.dest.clone(),
            subfiles: generate_subfile_dests(&outprefix, chunks.len() - 1, "temp.txt"),
        };

        // we use thread-local copy so we can drop the mutex lock
        let local_outputs = output_merge_lock.clone();
        drop(output_merge_lock);

        let src = &self.src.clone();
        let mut pool = ThreadPool::new(self.plp_params.threads);

        while !chunks.is_empty() {
            let worker = pool.get_available();
            let (id, chunk) = chunks.pop_front().unwrap();
            let refseq_handle = self
                .refseq
                .yield_handle(&chunk.name, self.plp_params.refseq.as_deref())?;

            worker.run(
                id,
                self.plp_params.clone(),
                chunk.clone(),
                src.clone(),
                self.output.clone(),
                local_outputs.get_writer(id)?,
                refseq_handle,
            );
        }

        // threadpool.install(|| {
        //     chunks.par_iter().enumerate().for_each(|(i, chunk)| {
        //         let refhandle = self
        //             .refseq
        //             .yield_handle(&chunk.name, self.plp_params.refseq.as_deref())
        //             .unwrap();

        //         let before = Instant::now();
        //         let mut worker = PileupWorker::new(self.plp_params.clone(), chunk.clone(), src.clone());
        //         let writer = local_outputs.get_writer(i).expect("failed to get writer");
        //         worker.run(self.output.clone(), writer, refhandle);

        //         info!(
        //             "Chunk {} / {} completed in {} seconds...",
        //             i,
        //             chunks.len(),
        //             before.elapsed().as_secs()
        //         );

        //         self.refseq.decrement_ref_usage(&chunk.name);
        //     });
        // });

        let main_writer = local_outputs.get_writer(0)?;
        local_outputs.merge(main_writer)?;

        Ok(())
    }
}
