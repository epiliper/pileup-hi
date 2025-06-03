use crate::overlap::{MapOverlaps, OverlapInsertResult, OverlapMap};
use crate::pileup::{cigar2rlen, CigarState, PileUp};
use rust_htslib::bam::Record;
use std::{cell::RefCell, rc::Rc};

pub struct ReadBuffer {
    pub rbuf: Vec<Rc<RefCell<PileUp>>>,
    pub overlaps: OverlapMap,
    pub len: usize,
}

#[derive(Debug)]
pub enum BufPushResult {
    AfterWindow(usize),
    Pushed(Rc<RefCell<PileUp>>),
    DifferentReference,
    Unmapped,
    MateFound(Rc<RefCell<PileUp>>),
}

impl<'a> ReadBuffer {
    pub fn c_to_next_window(&mut self, next_pos: i64, cur_pos: usize) -> usize {
        let next_pos = next_pos as usize;
        std::cmp::max(0, next_pos - (cur_pos + self.len - 1))
    }

    pub fn attempt_push(&'a mut self, r: &Record, pos: usize, tid: u32) -> BufPushResult {
        if r.is_unmapped() {
            return BufPushResult::Unmapped;
        }

        if r.tid() as u32 != tid {
            return BufPushResult::DifferentReference;
        }

        if cigar2rlen(&r) > self.len {
            self.len = cigar2rlen(&r);
        }

        if r.pos() as usize + self.len - 1 < pos {
            panic!(); // unsorted
        }

        if r.pos() as usize > pos + self.len - 1 {
            let window_start = self.c_to_next_window(r.pos(), pos);
            return BufPushResult::AfterWindow(window_start);
        }

        let cstate = CigarState {
            cig: r.cigar(),
            icig: 0,
            iseq: 0,
            bam_pos: r.pos() as u32,
        };

        let plp = PileUp {
            rec: r.clone(),
            indel: 0,
            cstate,
        };

        match self.overlaps.push(plp) {
            OverlapInsertResult::Inserted(plp_ref) => {
                return BufPushResult::MateFound(plp_ref);
            }

            OverlapInsertResult::Rejected(plp_obj) => {
                self.rbuf.push(Rc::new(RefCell::new(plp_obj)));
                let plp_ref = self.rbuf.get(self.rbuf.len() - 1).unwrap();
                return BufPushResult::Pushed(Rc::clone(plp_ref));
            }
        }
    }

    pub fn new() -> Self {
        let rbuf: Vec<Rc<RefCell<PileUp>>> = Vec::with_capacity(500);
        let overlaps = OverlapMap::new();
        let len = 0;

        Self {
            rbuf,
            overlaps,
            len,
        }
    }
}
