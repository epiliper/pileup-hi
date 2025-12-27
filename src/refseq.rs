use anyhow::Error;
use bio::io::{
    fasta,
    fasta::{FastaRead, Record},
};
use log::warn;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

pub struct FastaIndexedReader {
    inner: fasta::IndexedReader<File>,
    data: Vec<u8>,
}

pub struct FastaReader {
    inner: fasta::Reader<BufReader<File>>,
    data: Record,
}

pub trait ReadsFasta {
    fn read_to_bytes(&mut self, refname: &str) -> Result<Option<i64>, Error>;
    fn yield_seq(&self) -> Option<&[u8]>;
}

impl ReadsFasta for FastaIndexedReader {
    fn read_to_bytes(&mut self, refname: &str) -> Result<Option<i64>, Error> {
        // ref not found
        if self.inner.fetch_all(refname).is_err() {
            warn!("Unable to find ref {refname} in fasta. Proceeding without reference...");
            return Ok(None);
        };

        self.inner.read(&mut self.data)?;

        Ok(Some(self.data.len() as i64))
    }

    fn yield_seq(&self) -> Option<&[u8]> {
        if self.data.is_empty() {
            None
        } else {
            Some(self.data.as_slice())
        }
    }
}

impl ReadsFasta for FastaReader {
    fn read_to_bytes(&mut self, refname: &str) -> Result<Option<i64>, Error> {
        loop {
            self.inner.read(&mut self.data)?;
            if self.data.id() == refname {
                return Ok(Some(self.data.seq().len() as i64)); // found it
            } else if self.data.seq().is_empty() {
                return Ok(None); // read through all refs without finding one matching the given id
            }
        }
    }

    fn yield_seq(&self) -> Option<&[u8]> {
        if self.data.seq().is_empty() {
            None
        } else {
            Some(self.data.seq())
        }
    }
}

pub struct RefSeq {
    reader: Box<dyn ReadsFasta>,
}

impl RefSeq {
    // TODO: use regular fasta reader to avoid using rust_bio
    pub fn from_file(file: &str) -> Result<Self, Error> {
        let idx_name = format! {"{file}.fai"};
        let faidx = Path::new(&idx_name);

        if !faidx.exists() {
            let reader = FastaReader {
                inner: fasta::Reader::from_file(Path::new(&file))?,
                data: Record::new(),
            };

            Ok(Self {
                reader: Box::new(reader),
            })
        } else {
            let reader = FastaIndexedReader {
                inner: fasta::IndexedReader::from_file(&Path::new(&file))?,
                data: vec![],
            };

            Ok(Self {
                reader: Box::new(reader),
            })
        }
    }

    pub fn load_seq(&mut self, t_name: &str) -> Result<(), Error> {
        self.reader.read_to_bytes(t_name)?.unwrap_or(0);
        Ok(())
    }

    pub fn yield_seq(&self) -> Option<&[u8]> {
        self.reader.yield_seq()
    }
}
