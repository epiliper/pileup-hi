use crate::read_buf;
use crate::read_buf::CigarState;
use anyhow::{Context, Error};
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{HeaderView, Read, Reader};
use std::collections::VecDeque;

const UNINIT_POS: usize = usize::MAX - 1;
const UNINIT_TID: u32 = u32::MAX - 1;

pub struct PileUp {
    tid: u32,
    pos: usize,
    rbuf: read_buf::ReadBuffer,
    reader: Reader,
    header: HeaderView,
}

pub enum IterResult {
    ReferenceEnd,
    Generated,
    NoData,
}

#[derive(Debug, PartialEq, Eq)]
pub enum CigarResult {
    OutOfBounds(),
    Op(Cigar),
}

pub fn cigar_get_pos(cs: &mut CigarState, pos: u32, ipos: &mut i32) -> CigarResult {
    let cig = &cs.cig;
    let ncig = cig.len();
    while cs.bam_pos <= pos {
        if cs.icig >= ncig {
            return CigarResult::OutOfBounds();
        }

        let op = cig[cs.icig];
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

    CigarResult::OutOfBounds()
}

impl PileUp {
    pub fn new(bam_fname: &str, tid: Option<u32>, pos: Option<usize>) -> Result<Self, Error> {
        let tid = tid.unwrap_or(UNINIT_TID);
        let pos = pos.unwrap_or(UNINIT_POS);
        let reader = Reader::from_path(bam_fname)?;
        let mut rbuf = read_buf::ReadBuffer::new();
        let header = reader.header().clone();
        rbuf.pos = pos;
        rbuf.tid = tid;

        Ok(Self {
            tid,
            pos,
            rbuf,
            reader,
            header,
        })
    }

    pub fn fill_buffer(&mut self) -> Result<read_buf::BufPushResult, Error> {
        let mut ret: read_buf::BufPushResult = read_buf::BufPushResult::DifferentReference;
        let mut pos: i64;

        if self.rbuf.rbuf.is_empty() || self.pos == UNINIT_POS {
            let first = self.reader.records().next().context("no read")??;
            self.pos = first.pos() as usize;
            self.tid = first.tid() as u32;
            self.rbuf.pos = first.pos() as usize;
            self.rbuf.tid = first.tid() as u32;
        }

        for rec in self.reader.records() {
            let r = rec?;

            ret = self.rbuf.push(r);

            match ret {
                read_buf::BufPushResult::AfterWindow => {
                    break;
                }
                read_buf::BufPushResult::DifferentReference => break,
                _ => (),
            }
        }
        Ok(ret)
    }

    pub fn set_pileup(&mut self) {
        let mut ndel @ mut nins @ mut nbases = 0;
        let mut to_remove: VecDeque<usize> = VecDeque::new();

        for (i, r) in self.rbuf.rbuf.iter_mut().enumerate() {
            let mut ipos: i32 = -1;
            let ret = cigar_get_pos(&mut r.cstate, self.pos as u32, &mut ipos);
            match ret {
                CigarResult::Op(Cigar::Match(_)) => {
                    let base = r.rec.seq().encoded_base(ipos as usize).to_ascii_uppercase();
                    print! {" {base}"}
                    nbases += 1;
                }

                CigarResult::Op(Cigar::Ins(_)) => nins += 1,

                CigarResult::Op(Cigar::Del(_)) => {
                    if ipos != -1 {
                        ndel += 1;
                    }
                }

                CigarResult::OutOfBounds() => {
                    // println! {"{} {} {}", self.pos, r.rec.pos(), self.rbuf.len}
                    to_remove.push_back(i);
                }

                _ => panic!(),
            }
        }

        while let Some(i) = to_remove.pop_back() {
            self.rbuf.rbuf.swap_remove(i);
        }

        print! {"\n"}
    }

    pub fn next(&mut self) -> Result<IterResult, Error> {
        self.pos += 1;
        if self.pos != UNINIT_POS + 1
            && self.pos
                >= self
                    .header
                    .target_len(self.tid)
                    .context("Unable to get ref len")? as usize
        {
            self.tid += 1;
            self.pos = UNINIT_POS;
            if self.header.target_count() <= self.tid {
                Ok(IterResult::NoData)
            } else {
                Ok(IterResult::ReferenceEnd)
            }
        } else {
            let r = self.fill_buffer();
            println! {"{:?}", r}
            println! {"{}", self.pos}
            self.set_pileup();
            Ok(IterResult::Generated)
        }
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
