use crate::params::InputParams;
use crate::utils::has_index;

use anyhow::{Context, Error};
use rust_htslib::bam::{HeaderView, IndexedReader, Read, Reader, Record};

pub struct BamReader {
    inner: Box<dyn BamRead>,
    pub header: HeaderView,
    pub max_tid: i32,
    pub cur_ref: String,
}

impl BamReader {
    pub fn new(params: &InputParams) -> Result<Self, Error> {
        match has_index(&params.input)? {
            true => {
                // println! {"Found index for {}.", &params.input}
                let inner = IndexedReader::new(&params.input, params.threads)?;
                let header = inner.header().clone();
                let max_tid = header.target_count() as i32;
                let cur_ref = "UNINIT".to_string();
                return Ok(Self {
                    inner,
                    header,
                    max_tid,
                    cur_ref,
                });
            }

            false => {
                // println! {"No index found for {}. Using slower iteration...", &params.input}
                let inner = Reader::new(&params.input, params.threads)?;
                let header = inner.header().clone();
                let max_tid = header.target_count() as i32;
                let cur_ref = "UNINIT".to_string();
                return Ok(Self {
                    inner,
                    header,
                    max_tid,
                    cur_ref,
                });
            }
        }
    }

    pub fn read_no_alloc(&mut self, stored_read: &mut Record) -> Option<Result<(), Error>> {
        self.inner.read_no_alloc(stored_read)
    }

    pub fn init_to_ref(&mut self, tid: u32) -> Result<(), Error> {
        self.cur_ref = std::str::from_utf8(self.header.tid2name(tid))?.to_string();
        self.inner.init_to_ref(tid)
    }
}

pub trait BamRead {
    fn init_to_ref(&mut self, tid: u32) -> Result<(), Error>;
    fn new(input_file: &str, threads: usize) -> Result<Box<Self>, Error>
    where
        Self: Sized;
    fn read_no_alloc(&mut self, stored_read: &mut Record) -> Option<Result<(), Error>>;
}

impl BamRead for Reader {
    fn init_to_ref(&mut self, _tid: u32) -> Result<(), Error> {
        Ok(())
    }

    fn new(input_file: &str, threads: usize) -> Result<Box<Self>, Error>
    where
        Self: Sized,
    {
        let mut ret = Self::from_path(input_file)?;
        ret.set_threads(threads)?;
        Ok(Box::new(ret))
    }

    fn read_no_alloc(&mut self, stored_read: &mut Record) -> Option<Result<(), Error>> {
        self.read(stored_read)
            .map(|e| e.context("Failed to retrieve read"))
    }
}

impl BamRead for IndexedReader {
    fn init_to_ref(&mut self, tid: u32) -> Result<(), Error> {
        self.fetch((tid, 0, u32::MAX)).context("Failed to fetch")
    }

    fn new(input_file: &str, threads: usize) -> Result<Box<Self>, Error>
    where
        Self: Sized,
    {
        let mut ret = Self::from_path(input_file)?;
        ret.set_threads(threads)?;
        Ok(Box::new(ret))
    }

    fn read_no_alloc(&mut self, stored_read: &mut Record) -> Option<Result<(), Error>> {
        self.read(stored_read)
            .map(|e| e.context("Failed to retrieve read"))
    }
}
