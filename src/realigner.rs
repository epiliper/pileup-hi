#![allow(dead_code)]

use anyhow::Error;
use minimap2::{Aligner, Built, Mapping};
use rust_htslib::bam::{
    record::{Cigar, CigarString},
    Record,
};

use crate::pileup_iterator::UNINIT_POS;

pub type Remapper = Aligner<Built>;

#[derive(Debug)]
pub struct AlignmentPayload {
    pos: i32,
    mapq: u32,
    cigarstr: String,
}

pub struct Realigner {
    prev_w_start: i64,
    prev_w_end: i64,
    aligner: Remapper,
    refname: Vec<u8>,
}

/// Adjust an existing read with updated fields aligned from a realignment.
pub fn set_read_to_realign(r: &mut Record, alp: AlignmentPayload) {
    println! {"BEFORE ------------ MAPQ: {} START: {} CIGAR{}", r.mapq(), r.pos(), r.cigar()}
    r.set_mapq(alp.mapq as u8);
    r.set_pos(alp.pos as i64);
    r.set_cigar(Some(&parse_cigar_string(&alp.cigarstr)));
    println! {"AFTER ------------- MAPQ: {} START: {} CIGAR{}", r.mapq(), r.pos(), r.cigar()}
}

/// remove supplementary maps, get only one alignment per read.
/// this is meant to be called on a [Vec<Mapping>] from a single read
pub fn filter_maps(mut maps: Vec<Mapping>) -> AlignmentPayload {
    let mut aln: Vec<AlignmentPayload>;
    println! {"MAPS: {}", maps.len()}

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
    pub fn new(reference: AlignerReference, refname: Option<&str>) -> Result<Self, Error> {
        let aligner: Aligner<Built>;
        let aligner_build = Aligner::builder()
            .with_cigar()
            .sr()
            .with_sam_hit_only()
            .with_index_threads(num_cpus::get());

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

        Ok(Self {
            prev_w_start: 0,
            prev_w_end: UNINIT_POS,
            aligner,
            refname,
        })
    }

    pub fn realign_region(&mut self, mut records: Vec<Record>) -> Result<(), Error> {
        let mut aln: AlignmentPayload;
        let mut maps: Vec<Mapping>;

        for rec in records.iter_mut() {
            maps = self
                .aligner
                .map(
                    // rec.seq().as_bytes().as_slice(),
                    rec.seq().as_bytes().as_slice(),
                    true,
                    false,
                    None,
                    None,
                    Some(self.refname.as_slice()),
                )
                .expect("failed to realign");

            aln = filter_maps(maps);
            println! {"Alignment: {:?}", aln}
            set_read_to_realign(rec, aln);
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
        // let ref_reads = get_fastq_sequences("cDNA_reads.fq");
        let ref_records = bam_from_fastq("cDNA_reads.fq");

        let mut realigner = Realigner::new(AlignerReference::Fasta(&ref_file), None).unwrap();
        realigner.realign_region(ref_records).unwrap();
    }
}
