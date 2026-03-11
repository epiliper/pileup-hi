use anyhow::{Context, Error};
use rust_htslib::bam::{
    // record::{Cigar, CigarString},
    Record,
};

use crate::bamio::{BamDataSource, BamReader};
use crate::baq::realign_record;
use crate::engine::BUFWRITER_CAP;
use crate::params::RealignParams;
use std::io::{BufWriter, Write};

const N_CIGAR_OPS: usize = 9;

// const CMATCH: usize = 0;
// const CINS: usize = 1;
// const CDEL: usize = 2;
// const CREFSKIP: usize = 3;
// const CSOFTCLIP: usize = 4;
// const CHARDCLIP: usize = 5;
// const CPAD: usize = 6;
// const CEQUAL: usize = 7;
// const CDIFF: usize = 8;

pub type CigarOpCounter = [u32; N_CIGAR_OPS];

const POS_UNINIT: i64 = -1;

// ==================
// Realignment report
// ==================

// Current format for report:
// TID       POS        NUM_READS_TOTAL     NUM_READS_REALIGNED     FRAC_READS_REALIGNED     AVG_POSITION_SHIFT
//
// Followed by counts of each cigar string operator between original and new.
// Idea is to make patterns in read alignment drift more noticeable.
//

// Track read metrics before and after alignment
struct RealignmentReport {
    old_ops: CigarOpCounter,
    new_ops: CigarOpCounter,
    total_reads: usize,
    realigned_reads: usize,
    intbuf: itoa::Buffer,
    floatbuf: zmij::Buffer,
    writer: BufWriter<std::fs::File>,
    realign_report_thres: f32,
}

impl RealignmentReport {
    fn new(output_path: std::path::PathBuf, realign_report_thres: f32) -> Result<Self, Error> {
        Ok(Self {
            old_ops: [0; N_CIGAR_OPS],
            new_ops: [0; N_CIGAR_OPS],
            total_reads: 0,
            realigned_reads: 0,
            intbuf: itoa::Buffer::new(),
            floatbuf: zmij::Buffer::new(),
            writer: BufWriter::with_capacity(BUFWRITER_CAP, std::fs::File::create(output_path)?),
            realign_report_thres,
        })
    }

    fn clear(&mut self) {
        self.old_ops.fill(0);
        self.new_ops.fill(0);
        self.total_reads = 0;
        self.realigned_reads = 0;
    }

    // fn count_cigar_ops(&mut self, cig: &CigarString, old: bool) {
    //     let counter = if old { &mut self.old_ops } else { &mut self.new_ops };
    //     counter.fill(0);

    //     for c in cig {
    //         match c {
    //             Cigar::Match(_) => counter[CMATCH] += 1,
    //             Cigar::Ins(_) => counter[CINS] += 1,
    //             Cigar::Del(_) => counter[CDEL] += 1,
    //             Cigar::RefSkip(_) => counter[CREFSKIP] += 1,
    //             Cigar::SoftClip(_) => counter[CSOFTCLIP] += 1,
    //             Cigar::HardClip(_) => counter[CHARDCLIP] += 1,
    //             Cigar::Pad(_) => counter[CPAD] += 1,
    //             Cigar::Equal(_) => counter[CEQUAL] += 1,
    //             Cigar::Diff(_) => counter[CDIFF] += 1,
    //         }
    //     }
    // }

    fn write(&mut self, refname: &str, pos: i64) -> Result<(), Error> {
        // eprintln!("{} {}", self.realigned_reads, self.total_reads);
        if self.realigned_reads > 0 {
            std::process::exit(1);
        }
        if (self.realigned_reads as f32 / self.total_reads as f32) >= self.realign_report_thres {
            // TID
            self.writer.write_all(refname.as_bytes())?;
            self.writer.write_all(b"\t")?;

            // POS
            self.writer.write_all(self.intbuf.format(pos).as_bytes())?;
            self.writer.write_all(b"\t")?;

            // NUM_READS_TOTAL
            self.writer.write_all(self.intbuf.format(self.total_reads).as_bytes())?;
            self.writer.write_all(b"\t")?;

            // NUM_READS_REALIGNED
            self.writer
                .write_all(self.intbuf.format(self.realigned_reads).as_bytes())?;
            self.writer.write_all(b"\t")?;

            // FRAC READS REALIGNED
            self.writer.write_all(
                self.floatbuf
                    .format(self.realigned_reads as f32 / self.total_reads as f32)
                    .as_bytes(),
            )?;
            self.writer.write_all(b"\n")?;
        }

        self.clear();
        Ok(())
    }
}

// =============================================
// Realignment Output Types and Implementations
// =============================================
pub struct RealignmentReportDesc {
    output_path: std::path::PathBuf,
    realign_report_thres: f32,
}

pub struct RealignmentBamDesc<'a> {
    output_path: std::path::PathBuf,
    header: &'a rust_htslib::bam::HeaderView,
}

// Interface for a record realigner
pub trait RealignsBamRecords {
    fn realign_record(&mut self, r: &mut Record, ref_seq: &[u8], ref_len: i64) -> Result<(), Error>;
    fn default() -> Self;
}

// Create a default BAQ realigner
pub struct BAQRealigner {
    pub redo_baq: bool,
}

impl RealignsBamRecords for BAQRealigner {
    fn realign_record(&mut self, r: &mut Record, ref_seq: &[u8], ref_len: i64) -> Result<(), Error> {
        let flag = if self.redo_baq { 7 } else { 3 };
        realign_record(r, ref_seq, ref_len, flag)?;
        Ok(())
    }

    fn default() -> Self {
        Self { redo_baq: false }
    }
}

pub struct Realigner<A: RealignsBamRecords> {
    pub core: A,
    enabled: bool,
    report: Option<RealignmentReport>,
    bam_output: Option<rust_htslib::bam::Writer>,
    pos: i64,
}

impl<A: RealignsBamRecords> Realigner<A> {
    pub fn construct_from_args(params: &RealignParams, input_file: &BamDataSource) -> Result<Option<Self>, Error> {
        if params.no_baq {
            return Ok(None);
        }

        let header = BamReader::new(input_file, 1)?.header;

        let report_desc = if let Some(report_file) = params.report_file.as_ref() {
            assert!(params.report_realign_thres > 0.0);
            Some(RealignmentReportDesc {
                output_path: report_file.to_path_buf(),
                realign_report_thres: params.report_realign_thres,
            })
        } else {
            None
        };

        let bam_desc = if let Some(bam_output) = params.realign_bam_out.as_ref() {
            Some(RealignmentBamDesc {
                header: &header,
                output_path: bam_output.to_path_buf(),
            })
        } else {
            None
        };

        Some(Self::new(report_desc, bam_desc)).transpose()
    }

    pub fn enable(&mut self) {
        self.enabled = true;
    }

    pub fn disable(&mut self) {
        self.enabled = false;
    }

    pub fn new(
        report_desc: Option<RealignmentReportDesc>,
        bam_desc: Option<RealignmentBamDesc>,
    ) -> Result<Self, Error> {
        let report = if let Some(desc) = report_desc {
            Some(RealignmentReport::new(desc.output_path, desc.realign_report_thres)?)
        } else {
            None
        };

        let bam_output = if let Some(desc) = bam_desc {
            Some(rust_htslib::bam::Writer::from_path(
                desc.output_path,
                &rust_htslib::bam::Header::from_template(desc.header),
                rust_htslib::bam::Format::Bam,
            )?)
        } else {
            None
        };

        Ok(Self {
            core: A::default(),
            enabled: false,
            report,
            bam_output,
            pos: POS_UNINIT,
        })
    }

    fn write_report(&mut self, refname: &str) -> Result<(), Error> {
        if let Some(report) = self.report.as_mut() {
            report.write(refname, self.pos)
        } else {
            Ok(())
        }
    }

    fn write_bam(&mut self, r: &Record) -> Result<(), Error> {
        if let Some(writer) = self.bam_output.as_mut() {
            writer.write(r).context("failed to write bam record")
        } else {
            Ok(())
        }
    }

    fn count_read(&mut self, is_realigned: bool) {
        if let Some(report) = self.report.as_mut() {
            report.realigned_reads += is_realigned as usize;
            report.total_reads += 1;
        }
    }

    pub fn realign(&mut self, r: &mut Record, refname: &str, refseq: &[u8]) -> Result<(), Error> {
        if !self.enabled {
            return Ok(());
        }

        let prev_cig = fxhash::hash64(r.raw_cigar());
        let prev_pos = r.pos();

        if self.pos != prev_pos && prev_pos != POS_UNINIT {
            self.write_report(refname)?;
        }

        self.pos = prev_pos;

        self.core.realign_record(r, refseq, refseq.len() as i64)?;

        let new_cig = fxhash::hash64(r.raw_cigar());
        let new_pos = r.pos();

        // eprintln!("{} {} | {} {}", prev_cig, new_cig, prev_pos, new_pos);
        let is_realigned = (prev_pos != new_pos) || (prev_cig != new_cig);

        if is_realigned {
            self.write_bam(r)?;
        }

        self.count_read(is_realigned);

        Ok(())
    }
}
