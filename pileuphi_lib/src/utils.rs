use std::io::Write;
use std::{fs::OpenOptions, io::BufWriter};

use crate::errors::Error;
use crate::{alignment::PileupAlignment, bamio::OutputDataDest};

pub fn read_ends_before_pos(a: &PileupAlignment, pos: i64) -> bool {
    a.rec.pos() + a.cstate.read_len_from_cigar - 1 < pos
}

pub fn temp_fname(prefix: &str, suffix: &str, ext: &str) -> String {
    format!("{prefix}_{suffix}.{ext}")
}

pub struct OutputWriter {
    #[allow(dead_code)]
    dest: OutputDataDest,
    writer: BufWriter<Box<dyn std::io::Write>>,
}

impl OutputWriter {
    pub fn new(handle: &OutputDataDest, writer_cap: usize, lock: bool, append: bool) -> Result<Self, Error> {
        let writer: BufWriter<Box<dyn std::io::Write>> = match handle {
            OutputDataDest::File(file) => {
                let mut o = OpenOptions::new();
                let file = o.write(true).create(true).append(append).open(file)?;

                if lock {
                    file.lock()?;
                }

                BufWriter::with_capacity(writer_cap, Box::new(file))
            }

            OutputDataDest::Stdout => {
                if lock {
                    BufWriter::with_capacity(writer_cap, Box::new(std::io::stdout().lock()))
                } else {
                    BufWriter::with_capacity(writer_cap, Box::new(std::io::stdout()))
                }
            }
        };

        Ok(Self {
            dest: handle.clone(),
            writer,
        })
    }

    pub fn get(&mut self) -> &mut BufWriter<Box<dyn std::io::Write>> {
        &mut self.writer
    }

    pub fn flush(&mut self) -> Result<(), Error> {
        self.writer.flush().map_err(Error::from)
    }
}

pub fn has_index(bam_file: &str) -> Result<bool, Error> {
    let potential_index = format! {"{bam_file}.bai"};

    std::fs::exists(&potential_index).map_err(Error::from)
}
