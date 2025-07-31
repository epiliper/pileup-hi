#![allow(dead_code)]

use anyhow::Error;
use minimap2::{Aligner, Built, Mapping};
use rust_htslib::bam::{
    record::{Cigar, CigarString},
    Record,
};

use crate::pileup::PileupRef;

pub type Remapper = Aligner<Built>;

const DUMMY_REFERENCE: &[u8; 3] = b"ACT";

#[derive(Debug)]
pub struct AlignmentPayload {
    pos: i32,
    mapq: u32,
    cigarstr: String,
    tid: i32,
}

pub struct Realigner {
    aligner: Remapper,
    refname: Vec<u8>,
}

/// Adjust an existing read with updated fields aligned from a realignment.
pub fn set_read_to_realign(r: &mut Record, alp: AlignmentPayload) {
    // println! {"BEFORE ------------ MAPQ: {} START: {} CIGAR{}", r.mapq(), r.pos(), r.cigar()}
    r.set_mapq(alp.mapq as u8);
    r.set_pos(alp.pos as i64);
    r.set_cigar(Some(&parse_cigar_string(&alp.cigarstr)));
    r.set_tid(alp.tid);
    r.unset_unmapped();
    // println! {"AFTER ------------- MAPQ: {} START: {} CIGAR{}", r.mapq(), r.pos(), r.cigar()}
}

/// remove supplementary maps, get only one alignment per read.
/// this is meant to be called on a [Vec<Mapping>] from a single read
pub fn filter_maps(maps: &mut Vec<Mapping>) -> AlignmentPayload {
    let mut aln: Vec<AlignmentPayload>;
    // println!{"{:?}", maps}

    aln = maps
        .drain(..)
        .filter_map(|m| create_alignment_payload(m))
        .collect();

    assert_eq!(aln.len(), 1);
    aln.remove(0)
}

/// Create the necessary information for a read to be modified following realignment.
pub fn create_alignment_payload(map: Mapping) -> Option<AlignmentPayload> {
    if map.is_supplementary || !map.is_primary {
        return None;
    }
    if map.alignment.is_none() {
        return None;
    }

    Some(AlignmentPayload {
        pos: map.target_start,
        mapq: map.mapq,
        cigarstr: map.alignment.unwrap().cigar_str.unwrap(),
        tid: map.target_id,
    })
}

fn parse_cigar_string(cigar_str: &str) -> CigarString {
    let mut ops = Vec::new();
    let mut current_num = 0;
    for c in cigar_str.chars() {
        match c {
            '0'..='9' => current_num = current_num * 10 + (c as u32 - '0' as u32),
            _ => {
                let op = match c {
                    'M' => Cigar::Match(current_num),
                    'I' => Cigar::Ins(current_num),
                    'D' => Cigar::Del(current_num),
                    'N' => Cigar::RefSkip(current_num),
                    'S' => Cigar::SoftClip(current_num),
                    'H' => Cigar::HardClip(current_num),
                    'P' => Cigar::Pad(current_num),
                    'X' => Cigar::Diff(current_num),
                    '=' => Cigar::Equal(current_num),
                    _ => panic!("Invalid CIGAR operation: {}", c),
                };
                ops.push(op);
                current_num = 0;
            }
        }
    }

    CigarString(ops)
}

pub enum AlignerReference<'a> {
    Sequence(&'a [u8]),
    Fasta(&'a str),
}

impl Realigner {
    pub fn init_to_ref(
        &mut self,
        reference: AlignerReference,
        refname: Option<&str>,
    ) -> Result<(), Error> {
        let aligner: Aligner<Built>;
        let mut aligner_build = Aligner::builder()
            .sr()
            .with_cigar()
            .with_sam_out()
            .with_index_threads(num_cpus::get());

        aligner_build.mapopt.best_n = 1;
        aligner_build.idxopt.k = 5;
        aligner_build.idxopt.w = 5;

        match reference {
            AlignerReference::Fasta(file) => {
                aligner = aligner_build
                    .with_index(file, None)
                    .map_err(|e| Error::msg(e))?
            }
            AlignerReference::Sequence(bytes) => {
                aligner = aligner_build.with_seq(bytes).map_err(|e| Error::msg(e))?
            }
        };

        let refname = refname
            .unwrap_or("REF")
            .as_bytes()
            .iter()
            .map(|x| *x)
            .collect();

        self.aligner = aligner;
        self.refname = refname;

        Ok(())
    }

    pub fn build_empty() -> Result<Self, Error> {
        let aligner = Aligner::builder()
            .with_seq(b"ACT")
            .map_err(|e| Error::msg(e))?;

        Ok(Self {
            aligner,
            refname: "NONE".into(),
        })
    }

    pub fn realign(&mut self, seq: &[u8], outvec: &mut Vec<Mapping>) -> Result<(), Error> {
        *outvec = self
            .aligner
            .map(seq, true, false, None, None, Some(self.refname.as_slice()))
            .expect("failed to realign!");

        Ok(())
    }

    pub fn realign_region_record(&mut self, records: &mut Vec<Record>) -> Result<(), Error> {
        let mut aln: AlignmentPayload;
        let mut maps: Vec<Mapping> = vec![];

        for mut rec in records.iter_mut() {
            self.realign(rec.seq().as_bytes().as_slice(), &mut maps)?;
            if maps.is_empty() {
                // println!(
                //     "{} {:?} {:?} {:?}",
                //     self.aligner.has_index(),
                //     self.aligner.mapopt,
                //     rec,
                //     rec.cigar()
                // );

                continue;
            }

            aln = filter_maps(&mut maps);
            set_read_to_realign(&mut rec, aln);
        }

        Ok(())
    }

    pub fn realign_region_plp(&mut self, pileups: &mut Vec<PileupRef>) -> Result<(), Error> {
        let mut aln: AlignmentPayload;
        let mut maps: Vec<Mapping> = vec![];

        for plp in pileups.iter_mut() {
            let mut rec = &mut plp.borrow_mut().rec;
            self.realign(rec.seq().as_bytes().as_slice(), &mut maps)?;
            if maps.is_empty() {
                // println!(
                //     // anyhow::bail!(
                //     "{} {:?} {:?} {:?}",
                //     self.aligner.has_index(),
                //     self.aligner.mapopt,
                //     rec,
                //     rec.cigar()
                // );

                continue;
            }

            aln = filter_maps(&mut maps);
            set_read_to_realign(&mut rec, aln);
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    use std::fs::File;
    use std::io::{BufRead, BufReader};

    use rust_htslib::bam::Record;

    fn get_test_dir() -> PathBuf {
        let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        path.push("test_data");
        path
    }

    fn get_test_file(test_file: &str) -> PathBuf {
        let mut path = get_test_dir();
        path.push(test_file);
        path
    }

    fn get_fastq_sequences(test_file: &str) -> Vec<Vec<u8>> {
        let mut out = vec![];
        let f = get_test_file(test_file);
        let mut lines = BufReader::new(File::open(f).unwrap()).lines();

        while let (Some(_header), Some(seq), Some(_plus), Some(_qual)) =
            (lines.next(), lines.next(), lines.next(), lines.next())
        {
            out.push(seq.unwrap().as_bytes().into());
        }

        out
    }

    fn bam_from_fastq(test_file: &str) -> Vec<Record> {
        let reads = get_fastq_sequences(test_file);
        let records = reads
            .iter()
            .map(|r| {
                let mut record = Record::new();
                record.set(b"4", None, r, vec![255 as u8; r.len()].as_slice());
                record
            })
            .collect();

        records
    }

    #[test]
    fn test1() {
        let ref_file = get_test_file("cDNA.fasta").to_str().unwrap().to_string();
        let mut ref_records = bam_from_fastq("cDNA_reads.fq");

        let mut realigner = Realigner::build_empty().unwrap();

        realigner
            .init_to_ref(AlignerReference::Fasta(&ref_file), None)
            .unwrap();

        realigner.realign_region_record(&mut ref_records).unwrap();

        println!{"{:?}", ref_records}
        for r in &ref_records {
            assert!(!r.is_unmapped());
            println! {"{:?}, {:?}", r, r.cigar()}
        }
    }

    #[test]
    fn test2() {
        let ref_file = get_test_file("hiv.fasta").to_str().unwrap().to_string();
        let mut ref_records = bam_from_fastq("hiv_reads.fq");

        let mut realigner = Realigner::build_empty().unwrap();

        realigner
            .init_to_ref(AlignerReference::Fasta(&ref_file), None)
            .unwrap();

        realigner.realign_region_record(&mut ref_records).unwrap();


        println!{"{:?}", ref_records}
        for r in &ref_records {
            println! {"{:?}, {:?}", r, r.cigar()}
            assert!(!r.is_unmapped());
        }
    }
}
