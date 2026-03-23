use anyhow::Error;
use bio::io::{
    fasta,
    fasta::{FastaRead, Record},
};
use log::warn;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::sync::{Arc, Mutex, MutexGuard};

struct FastaIndexedReader {
    inner: fasta::IndexedReader<File>,
}

struct FastaReader {
    inner: fasta::Reader<BufReader<File>>,
}

pub type RefSeqHandle = Option<Arc<Vec<u8>>>;

pub trait ReadsFasta {
    fn read_to_bytes(&mut self, refname: &str) -> Result<Option<Vec<u8>>, Error>;
}

impl ReadsFasta for FastaIndexedReader {
    fn read_to_bytes(&mut self, refname: &str) -> Result<Option<Vec<u8>>, Error> {
        // ref not found
        if self.inner.fetch_all(refname).is_err() {
            warn!("Unable to find ref {refname} in fasta. Proceeding without reference...");
            return Ok(None);
        };

        let mut output = Vec::new();
        self.inner.read(&mut output)?;

        Ok(Some(output))
    }
}

impl ReadsFasta for FastaReader {
    fn read_to_bytes(&mut self, refname: &str) -> Result<Option<Vec<u8>>, Error> {
        let mut record: Record = Default::default();

        loop {
            self.inner.read(&mut record)?;
            if record.id() == refname {
                return Ok(Some(record.seq().to_vec())); // found it
            } else if record.seq().is_empty() {
                warn!("Unable to find ref {refname} in fasta. Proceeding without reference...");
                return Ok(None); // read through all refs without finding one matching the given id
            }
        }
    }
}

//////////////////////////////////////////

enum RefSeqSlotMarker {
    // We've loaded this reference sequence before and freed it.
    ReadBefore,

    // A number of threads are currently using this reference sequence.
    InUse(usize),
}

struct RefSeqSlot {
    marker: RefSeqSlotMarker,
    data: RefSeqHandle,
}

// Holds reference information requested by any number of processing threads, giving threads
// read-only access to a reference on demand. Responsible for freeing unused references and loading
// new ones.
pub struct RefSeq {
    data: Arc<Mutex<HashMap<String, RefSeqSlot>>>,
}

impl RefSeq {
    pub fn get_reader(file: &str) -> Result<Box<dyn ReadsFasta>, Error> {
        // TODO: use regular fasta reader to avoid using rust_bio
        let idx_name = format! {"{file}.fai"};
        let faidx = Path::new(&idx_name);

        let reader: Box<dyn ReadsFasta> = if !faidx.exists() {
            Box::new(FastaReader {
                inner: fasta::Reader::from_file(Path::new(&file))?,
            })
        } else {
            Box::new(FastaIndexedReader {
                inner: fasta::IndexedReader::from_file(&Path::new(&file))?,
            })
        };
        Ok(reader)
    }

    pub fn load_seq(file_name: &str, ref_name: &str) -> Result<RefSeqHandle, Error> {
        let mut reader = RefSeq::get_reader(file_name)?;
        let read = reader.read_to_bytes(ref_name)?;
        if let Some(bytes) = read {
            Ok(Some(Arc::new(bytes)))
        } else {
            Ok(None)
        }
    }

    pub fn new() -> Self {
        Self {
            data: Arc::new(Mutex::new(HashMap::new())),
        }
    }

    fn lock_data(&self) -> MutexGuard<'_, HashMap<String, RefSeqSlot>> {
        self.data
            .lock()
            .expect("failed to acquire mutex lock for reference data")
    }

    pub fn yield_refseq(&self, ref_name: &str) -> RefSeqHandle {
        let lock = self.lock_data();

        if let Some(slot) = lock.get(ref_name) {
            if let Some(ref data) = slot.data {
                return Some(Arc::clone(data));
            }
        }
        None
    }

    pub fn yield_handle(&self, ref_name: &str, file_name: Option<&str>) -> Result<RefSeqHandle, Error> {
        if let Some(fname) = file_name {
            let mut lock = self.lock_data();

            if let Some(slot) = lock.get_mut(ref_name) {
                match slot.marker {
                    RefSeqSlotMarker::ReadBefore => {
                        warn!("Loading refseq {} freed previously...", ref_name);
                        slot.marker = RefSeqSlotMarker::InUse(1);
                        slot.data = RefSeq::load_seq(fname, ref_name)?;
                    }

                    RefSeqSlotMarker::InUse(ref mut x) => {
                        *x += 1;
                    }
                }
            } else {
                // haven't loaded
                lock.entry(ref_name.to_string()).insert_entry(RefSeqSlot {
                    data: RefSeq::load_seq(fname, ref_name)?,
                    marker: RefSeqSlotMarker::InUse(1),
                });
            }

            drop(lock);
            Ok(self.yield_refseq(ref_name))
        } else {
            Ok(None)
        }
    }
}
