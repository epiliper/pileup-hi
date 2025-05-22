use crate::read_buf;
use crate::read_buf::{CigarState, CIG_POS_UNINIT};
use anyhow::{Context, Error};
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{ext::BamRecordExtensions, HeaderView, Read, Reader, Record};
use std::collections::VecDeque;

const UNINIT_POS: usize = usize::MAX - 1;
const UNINIT_TID: u32 = u32::MAX - 1;

const LAST_POS: u8 = b'$';
const FIRST_POS: u8 = b'^';

const INSERTION: u8 = b'+';
const DELETION: u8 = b'-';

const F_MATCH: u8 = b'.';
const R_MATCH: u8 = b',';

pub struct PileUp {
    tid: u32,
    pos: usize,
    rbuf: read_buf::ReadBuffer,
    reader: Reader,
    header: HeaderView,
    ref_seq: Option<Vec<u8>>,
}

pub enum IterResult {
    ReferenceEnd,
    Generated,
    NoData,
}

#[derive(Debug, PartialEq, Eq)]
pub enum CigarResult {
    BeforePos(),
    Op(Cigar),
    BaseEmpty(),
}

pub fn get_base_pileup(
    cs: &CigarState,
    r: &Record,
    ipos: u32,
    pos: usize,
    kstring_buf: &mut Vec<u8>,
    ref_base: u8,
) {
    let ipos = ipos as usize;
    let iseq = cs.iseq as usize;
    let bam_pos = cs.bam_pos as usize;

    if pos == r.reference_end() as usize - 1 {
        kstring_buf.push(LAST_POS);
    }

    if pos == bam_pos {
        kstring_buf.push(FIRST_POS)
    }

    let cur_base = r.seq()[ipos];
    let base: u8;

    if ref_base != cur_base {
        base = cur_base;
    } else {
        if r.is_reverse() {
            base = R_MATCH;
        } else {
            base = F_MATCH;
        }
    };

    kstring_buf.push(base);
}

/// Get the base of a read at a specific reference coordinate, considering the cigar of the
/// alignment. This is necessary to adjust for clipped bases, insertions, and deletions that would
/// not be accounted for when blindly indexing into the read sequence.
///
/// This function will traverse the cigar string up until or after the queried reference position,
/// will adjust the read start and end according to clips, deletions, and insertions, and return
/// whether the base in the read aligned to the reference is absent (e.g. softclip) or present.
// pub fn cigar_get_pos(cs: &mut CigarState, pos: u32, ipos: &mut i32) -> CigarResult {
//     // read not yet processed
//     let mut k: usize = 0;
//     if cs.icig == CIG_POS_UNINIT {
//         cs.icig = 0;
//         match cs.cig[0] {
//             Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) => (),

//             _ => {
//                 let n_cig = cs.cig.len();
//                 for i in 0..n_cig {
//                     k = i;
//                     match cs.cig[k] {
//                         Cigar::Match(_)
//                         | Cigar::Del(_)
//                         | Cigar::RefSkip(_)
//                         | Cigar::Equal(_)
//                         | Cigar::Diff(_) => break,
//                         Cigar::Ins(l) | Cigar::SoftClip(l) => *ipos += l as i32,
//                         _ => (),
//                     }
//                     assert!(cs.icig < n_cig);
//                     cs.icig = k;
//                 }
//             }
//         }
//     } else {
//         let mut op = cs.cig[cs.icig];
//         // the queried position is not covered by this op
//         // so we jump to the next one if it exists
//         if pos - cs.bam_pos >= op.len() {
//             assert!(cs.icig < cs.cig.len());
//             op = cs.cig[cs.icig + 1];
//             // if the next op is anything but clipped, and the current op is also not clipped, then
//             // increment sequence index/len
//             //
//             // this way, if we have CLIP->NON_CLIP or NON_CLIP->CLIP, we avoid adding incrementing
//             // the seq index by the amount of softclipped bases when we jump to the next op.
//             match op {
//                 Cigar::Match(_)
//                 | Cigar::Del(_)
//                 | Cigar::RefSkip(_)
//                 | Cigar::Equal(_)
//                 | Cigar::Diff(_) => match cs.cig[cs.icig] {
//                     Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => cs.iseq += l,
//                     _ => (),
//                 },

//                 _ => (),
//             }

//             cs.icig += 1;

//             // =============================================================================
//             // The next op is a clip of some kind, so increment sequence index by length and
//             // proceed to find the next non-clip op.
//             // We do this to avoid incrementing the sequence index with non-existent "bases"
//             // spanning the clip.
//             // =============================================================================
//         } else {
//             match cs.cig[cs.icig] {
//                 Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => cs.iseq += l,
//                 _ => (),
//             };

//             // we need to find the next non-clip
//             // we also increment the ref index to account for soft clips
//             for i in cs.icig + 1..cs.cig.len() {
//                 k = i;
//                 match cs.cig[i] {
//                     Cigar::Match(_)
//                     | Cigar::Del(_)
//                     | Cigar::RefSkip(_)
//                     | Cigar::Equal(_)
//                     | Cigar::Diff(_) => break, // non-clip found
//                     Cigar::Ins(l) | Cigar::SoftClip(l) => cs.iseq += l, // no bases, add to
//                     _ => (),
//                 }
//             }

//             cs.icig = k;
//         }
//         assert!(cs.icig < cs.cig.len());
//     }

//     // now to get pileup
//     let mut op: Cigar;
//     op = cs.cig[cs.icig];

//     // queried pos is covered by next cig op
//     if cs.bam_pos + op.len() - 1 == pos && cs.icig + 1 < cs.cig.len() {
//         let mut next_op = cs.cig[cs.icig + 1];

//         // if next_op == Cigar::Ins(next_op.len()) {
//         //     // update indel
//         //     for k in cs.icig + 2..cs.cig.len() {
//         //         next_op = cs.cig[k];
//         //         if
//         //     }
//         // }
//     }

//     match op {
//         Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
//             // update pos
//         },

//         Cigar::Del(l) | Cigar::RefSkip(l) => {

//         }
//     }

//     todo!();
// }

pub fn cigar_get_pos(cs: &mut CigarState, pos: u32, ipos: &mut i32) -> CigarResult {
    let cig = &cs.cig;
    let ncig = cig.len();
    let mut op: Cigar;
    while cs.bam_pos <= pos {
        if cs.icig >= ncig {
            return CigarResult::BeforePos();
        }

        op = cig[cs.icig];
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let end_pos = cs.bam_pos + len - 1;

                if end_pos < pos {
                    cs.bam_pos += len;
                    cs.iseq += len;
                    cs.icig += 1;
                    continue;
                }

                *ipos = pos as i32 - cs.bam_pos as i32 + cs.iseq as i32;
                if end_pos == pos && cs.icig + 1 < ncig {
                    let next_op = cig[cs.icig + 1];

                    match next_op {
                        Cigar::Ins(_) => return CigarResult::Op(next_op),
                        Cigar::Del(_) => return CigarResult::Op(next_op),
                        _ => (),
                    }
                }
                return CigarResult::Op(Cigar::Match(len));
            }

            Cigar::Ins(len) | Cigar::SoftClip(len) => {
                cs.iseq += len;
                cs.icig += 1;
                continue;
            }

            Cigar::Del(len) => {
                let end_pos = cs.bam_pos + len - 1;
                if end_pos < pos {
                    cs.bam_pos += len;
                    cs.icig += 1;
                    continue;
                }

                *ipos -= 1;
                return CigarResult::Op(op);
            }

            Cigar::RefSkip(len) => {
                let end_pos = cs.bam_pos + len - 1;
                if end_pos < pos {
                    cs.bam_pos += len;
                    cs.icig += 1;
                    continue;
                }
            }
            _ => (),
        }
    }

    CigarResult::BaseEmpty()
}

impl PileUp {
    pub fn new(bam_fname: &str, tid: Option<u32>, pos: Option<usize>) -> Result<Self, Error> {
        let tid = tid.unwrap_or(UNINIT_TID);
        let pos = pos.unwrap_or(UNINIT_POS);
        let reader = Reader::from_path(bam_fname)?;
        let rbuf = read_buf::ReadBuffer::new();
        let header = reader.header().clone();
        let ref_seq = None;

        Ok(Self {
            tid,
            pos,
            rbuf,
            reader,
            header,
            ref_seq,
        })
    }

    pub fn fill_buffer(&mut self) -> Result<read_buf::BufPushResult, Error> {
        let mut ret: read_buf::BufPushResult = read_buf::BufPushResult::DifferentReference;

        // if our first read for the reference
        if self.pos == UNINIT_POS {
            if let Some(next_read) = self.reader.records().next() {
                let next_read = next_read?;
                self.pos = next_read.pos() as usize;
                self.tid = next_read.tid() as u32;
                let ret = self.rbuf.push(next_read, self.pos, self.tid);
                assert!(ret == read_buf::BufPushResult::Pushed);
            } else {
                return Ok(read_buf::BufPushResult::DifferentReference);
            }
        }

        let mut prev_pos = i64::MIN;
        for rec in self.reader.records() {
            let r = rec?;

            if r.pos() < prev_pos {
                panic!("UNSORTED BAM!")
            }

            prev_pos = r.pos();

            ret = self.rbuf.push(r, self.pos, self.tid);

            match ret {
                read_buf::BufPushResult::AfterWindow => break,
                read_buf::BufPushResult::DifferentReference => break,
                _ => (),
            }
        }
        Ok(ret)
    }

    pub fn set_pileup(&mut self) -> Result<(), Error> {
        let mut ndel @ mut nins @ mut nbases = 0;
        let mut to_remove: VecDeque<usize> = VecDeque::new();
        let mut seq: Vec<u8> = Vec::new();
        let ref_base = match &self.ref_seq {
            Some(seq) => seq[self.pos],
            None => b'N',
        };

        for (i, r) in self.rbuf.rbuf.iter_mut().enumerate() {
            let mut ipos: i32 = -1;
            let ret = cigar_get_pos(&mut r.cstate, self.pos as u32, &mut ipos);
            match ret {
                CigarResult::Op(Cigar::Match(_)) => {
                    get_base_pileup(&r.cstate, &r.rec, ipos as u32, self.pos, &mut seq, ref_base);

                    nbases += 1;
                }

                CigarResult::Op(Cigar::Ins(_)) => nins += 1,

                CigarResult::Op(Cigar::Del(_)) => {
                    if ipos != -1 {
                        ndel += 1;
                    }
                }

                CigarResult::BeforePos() => {
                    to_remove.push_back(i);
                }

                CigarResult::BaseEmpty() => (),
                _ => panic!(),
            }
        }

        print! {"{}\t{}\t{}\t{}\t", std::str::from_utf8(self.header.tid2name(self.tid))?, self.pos + 1, char::from(ref_base), nbases + nins + ndel }
        print! {"{}", std::str::from_utf8(&seq)?}

        while let Some(i) = to_remove.pop_back() {
            self.rbuf.rbuf.swap_remove(i);

            // using this option slows everything to a crawl: O(n) removal
            // self.rbuf.rbuf.remove(i);
        }

        print! {"\n"}

        Ok(())
    }

    pub fn next(&mut self) -> Result<IterResult, Error> {
        if self.pos != UNINIT_POS {
            self.pos += 1;
        }

        if self.tid == UNINIT_TID {
            self.tid = 0;
        }

        let ref_len = self.header.target_len(self.tid).context("no ref len")? as usize;

        // reached end of reference, so increment tid
        if self.pos != UNINIT_POS && self.pos >= ref_len {
            self.tid += 1;
            self.pos = UNINIT_POS;

            if self.header.target_count() > self.tid {
                // no more references
                return Ok(IterResult::NoData);
            } else {
                // another reference left
                return Ok(IterResult::ReferenceEnd);
            }
        }
        let _r = self.fill_buffer();
        self.set_pileup()?;
        Ok(IterResult::Generated)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::{record::CigarString, Record};

    #[test]
    pub fn cig_test1() {
        let cig = Vec::from([Cigar::Match(76)]);
        assert_eq!(cig[0].len(), 76)
    }

    #[test]
    pub fn cig_test2() {
        let mut record = Record::new();
        record.set(
            b"read1",
            Some(&CigarString(vec![Cigar::Match(4), Cigar::Equal(1)])),
            b"AAAAG",
            b"#####",
        );

        record.set_pos(1);

        let mut ipos = 0;

        let mut cstate = CigarState {
            cig: record.cigar(),
            icig: 0,
            iseq: 0,
            bam_pos: 1,
        };

        let mut ret = cigar_get_pos(&mut cstate, 4, &mut ipos);
        assert_eq!(ret, CigarResult::Op(Cigar::Match(4)));
        ret = cigar_get_pos(&mut cstate, 5, &mut ipos);
        assert_eq!(ret, CigarResult::Op(Cigar::Match(1)))
    }

    #[test]
    pub fn cig_test3() {
        let mut record = Record::new();
        record.set(
            b"read1",
            Some(&CigarString(vec![
                Cigar::Match(4),
                Cigar::Equal(1),
                Cigar::Ins(2),
                Cigar::Match(3),
            ])),
            b"AAAAGTTTTT",
            b"##########",
        );

        record.set_pos(104);

        let mut ipos = 0;

        let mut cstate = CigarState {
            cig: record.cigar(),
            icig: 0,
            iseq: 0,
            bam_pos: 104,
        };

        let mut ret = cigar_get_pos(&mut cstate, 107, &mut ipos);
        assert_eq!(ret, CigarResult::Op(Cigar::Match(4)));

        ret = cigar_get_pos(&mut cstate, 108, &mut ipos);
        assert_eq!(ret, CigarResult::Op(Cigar::Ins(2)));

        ret = cigar_get_pos(&mut cstate, 109, &mut ipos);
        assert_eq!(ret, CigarResult::Op(Cigar::Match(3)));
    }
}
