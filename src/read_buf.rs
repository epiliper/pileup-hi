use rust_htslib::bam::{
    record::{Cigar, CigarStringView},
    Record,
};

pub const CIG_POS_UNINIT: usize = usize::MAX - 1;

pub struct PileUpRecord {
    pub rec: Record,
    pub cstate: CigarState,
}

pub struct ReadBuffer {
    pub rbuf: Vec<PileUpRecord>,
    pub len: usize,
}

#[derive(Debug, Eq, PartialEq)]
pub enum BufPushResult {
    AfterWindow((Record, usize)),
    Pushed,
    DifferentReference((Record, usize)),
    Unmapped,
}

pub struct CigarState {
    pub cig: CigarStringView,
    pub icig: usize,  // position in cigar string
    pub iseq: u32,    // position in read sequence that corresponds to cigar pos
    pub bam_pos: u32, // ref coord of first base
}

pub fn cigar2rlen(r: &Record) -> usize {
    let mut len = 0;
    for op in &r.cigar() {
        match op {
            Cigar::Match(l)
            | Cigar::Del(l)
            | Cigar::RefSkip(l)
            | Cigar::Equal(l)
            | Cigar::Diff(l) => len += l,
            _ => (),
        }
    }

    len as usize
}

impl ReadBuffer {
    pub fn c_to_next_window(&mut self, next_pos: i64, cur_pos: usize) -> usize {
        let next_pos = next_pos as usize;
        std::cmp::max(0, next_pos - (cur_pos + self.len - 1))
    }

    pub fn next_c(&mut self, pos: i64) -> usize {
        let pos = pos as usize;

        std::cmp::max(0, pos - self.len - 1)
    }

    pub fn push(&mut self, r: Record, pos: usize, tid: u32) -> BufPushResult {
        if r.is_unmapped() {
            return BufPushResult::Unmapped;
        }

        if r.tid() as u32 != tid {
            return BufPushResult::DifferentReference((r, usize::MAX));
        }

        if cigar2rlen(&r) > self.len {
            self.len = cigar2rlen(&r);
        }

        if r.pos() as usize + self.len - 1 < pos {
            panic!(); // unsorted
        }

        if r.pos() as usize > pos + self.len - 1 {
            let window_start = self.c_to_next_window(r.pos(), pos);
            return BufPushResult::AfterWindow((r, window_start));
        }

        let cstate = CigarState {
            cig: r.cigar(),
            icig: 0,
            iseq: 0,
            bam_pos: r.pos() as u32,
        };

        self.rbuf.push(PileUpRecord { rec: r, cstate });
        BufPushResult::Pushed
    }

    pub fn new() -> Self {
        let rbuf: Vec<PileUpRecord> = Vec::with_capacity(500);
        let len = 0;

        Self { rbuf, len }
    }
}
