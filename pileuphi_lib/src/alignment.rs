use rust_htslib::bam::{
    record::{Cigar, CigarStringView},
    Record,
};

use std::cell::UnsafeCell;

pub const CIGAR_STATE_UNINIT: usize = usize::MAX - 10;

/// analogous to bam_pileup1_t
pub struct PileupAlignment {
    pub rec: Record,
    pub cstate: CigarState,
    pub indel: i32,
    pub refskip: bool,
    pub del: bool,
    pub head: bool,
    pub tail: bool,
    pub qpos: usize,
    pub cigar_index: usize,
}

impl PileupAlignment {
    pub fn new(r: Record, cstate: CigarState) -> Self {
        Self {
            rec: r,
            cstate,
            indel: 0,
            refskip: false,
            del: false,
            head: false,
            tail: false,
            qpos: 0,
            cigar_index: 0,
        }
    }
}

// pub type PileupAlignmentRef = Rc<UnsafeCell<PileupAlignment>>;
pub struct PileupAlignmentRef {
    value: UnsafeCell<PileupAlignment>,
}

impl PileupAlignmentRef {
    #[allow(clippy::mut_from_ref)]
    pub unsafe fn get(&self) -> &mut PileupAlignment {
        &mut (*self.value.get())
    }

    pub fn new(plp: PileupAlignment) -> Self {
        Self {
            value: UnsafeCell::new(plp),
        }
    }
}

pub struct CigarState {
    pub cig: CigarStringView,
    pub icig: usize,  // position in cigar string
    pub iseq: usize,  // position in read sequence that corresponds to cigar pos
    pub bam_pos: i64, // ref coord of first base
    pub read_len_from_cigar: i64,
}

pub fn cigar2rlen(r: &Record) -> i64 {
    let mut len = 0;
    for op in &r.cigar() {
        match op {
            Cigar::Match(l) | Cigar::Del(l) | Cigar::RefSkip(l) | Cigar::Equal(l) | Cigar::Diff(l) => len += l,
            _ => (),
        }
    }

    len as i64
}
