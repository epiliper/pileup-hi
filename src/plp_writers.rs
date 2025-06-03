use crate::pileup::CigarState;
use crate::refseq::RefSeq;
use anyhow::Error;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{ext::BamRecordExtensions, Record};
use std::io::Write;

const LAST_POS: u8 = b'$';
const FIRST_POS: u8 = b'^';

const F_MATCH: u8 = b'.';
const R_MATCH: u8 = b',';

pub fn get_base(mut cur_base: u8, is_reverse: bool) -> u8 {
    match is_reverse {
        false => cur_base.make_ascii_uppercase(),
        true => cur_base.make_ascii_lowercase(),
    }

    cur_base
}

pub fn get_base_to_ref(
    mut cur_base: u8,
    ref_coord: u64,
    refseq: Option<&RefSeq>,
    is_reverse: bool,
) -> Result<u8, Error> {
    if let Some(refseq) = refseq {
        let ref_base = refseq.get_base(ref_coord)?;
        if ref_base == cur_base {
            if is_reverse {
                cur_base = R_MATCH;
            } else {
                cur_base = F_MATCH;
            }
        }
        Ok(cur_base)
    } else {
        Ok(get_base(cur_base, is_reverse))
    }
}

pub fn write_match(
    cs: &CigarState,
    r: &Record,
    ipos: u32,
    pos: usize,
    seq_buf: &mut Vec<u8>,
    qual_buf: &mut Vec<u8>,
    refseq: Option<&RefSeq>,
) -> Result<(), Error> {
    let ipos = ipos as usize;
    let bam_pos = cs.bam_pos as usize;

    if pos == bam_pos {
        seq_buf.push(FIRST_POS);
        seq_buf.push(r.mapq() + 33);
    }

    let cur_base = r.seq()[ipos];

    let base = get_base_to_ref(cur_base, pos as u64, refseq, r.is_reverse())?;

    let cur_qual = r.qual()[ipos] + 33;

    seq_buf.push(base);

    if pos == r.reference_end() as usize - 1 {
        seq_buf.push(LAST_POS);
    }

    qual_buf.push(cur_qual);

    Ok(())
}

pub fn write_del(pos: usize, seq_buf: &mut Vec<u8>, del_len: usize) -> Result<(), Error> {
    write!(seq_buf, "-{}", del_len)?;
    for _ in pos..pos + del_len {
        seq_buf.push(b'N')
    }
    Ok(())
}

pub fn write_ins(
    cs: &CigarState,
    r: &Record,
    ipos: u32,
    seq_buf: &mut Vec<u8>,
    qual_buf: &mut Vec<u8>,
) -> Result<(), Error> {
    let mut k = cs.icig + 1;
    let ncig = cs.cig.len();
    let ipos = ipos + 1;
    while k < ncig {
        match cs.cig[k] {
            Cigar::Pad(l) => {
                seq_buf.extend(std::iter::repeat_n(b'*', l as usize));
            }

            Cigar::Ins(l) => {
                write!(seq_buf, "+{}", l)?;
                let (s, e) = (ipos as usize, (ipos + l) as usize);
                for i in s..e {
                    let base = get_base(r.seq()[i], r.is_reverse());
                    seq_buf.push(base);
                    qual_buf.push(r.qual()[i]);
                }
            }

            _ => break,
        }

        k += 1;
    }

    Ok(())
}
