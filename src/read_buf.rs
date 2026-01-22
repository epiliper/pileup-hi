use crate::alignment::{cigar2rlen, CigarState, PileupAlignment, PileupAlignmentRef, CIGAR_STATE_UNINIT};
use crate::cigar_resolve::resolve_cigar;
use crate::overlap::{MapOverlaps, OverlapMap};
use anyhow::Error;
use log::error;
use rust_htslib::bam::Record;
use std::{cell::RefCell, collections::HashMap, rc::Rc};

/// A container storing reads to be used in pileup generation. When requested to intake a
/// read given a particular tid/coordinate, the buffer will check that the read overlaps with the
/// coordinate and upon intake, perform overlap correction if enabled.
///
/// The ReadBuffer can be capped to hold a maximum amount of reads, tracked via the 'max_depth'
/// field.
///
/// The buffer is actually a pair of alternating dynamically-resizing arrays: buffer 1 (Self::rbuf) is drained
/// during pileup generation, and the items being drained and retained are added to buffer 2
/// (Self::backup_buf). This way, we can 'remove' reads that no longer overlap with the queried
/// position from the middle of the buffer without having to do a nasty left-shift of the entire
/// array.
///
/// It is mandatory to call Self::reset() when draining the buffer. See generate_pileup() in
/// pileup_iterator.rs.
///
/// If you're wondering, using std::collections::VecDeque is signficantly slower in my experience.
pub struct ReadBuffer {
    pub rbuf: Vec<PileupAlignmentRef>,
    pub len: i64,
    pub backup_buf: Vec<PileupAlignmentRef>,
    pub overlap_map: Option<OverlapMap>,
    pub depth: usize,
    pub max_depth: usize,
    pub head: BufferBoundary,
    pub tail: BufferBoundary,
}

pub struct BufferBoundary {
    pub tid: i32,
    pub pos: i64,
}

impl BufferBoundary {
    fn new() -> Self {
        Self {
            tid: i32::MIN,
            pos: i64::MIN,
        }
    }

    pub fn set_to_rec(&mut self, r: &Record) {
        self.tid = r.tid();
        self.pos = r.pos();
    }
}

pub enum BufPushResult {
    Pushed,
    DifferentReference,
    Unmapped,
    MaxDepthMet,
    BeforePos,
}

impl ReadBuffer {
    #[inline(always)]
    /// Store the read we've decided to push in the buffer, wrapping it in an Rc pointer and
    /// performing mate overlap-correction if requested by caller.
    pub fn store(&mut self, r: &Record, read_len_from_cigar: i64, pos: i64) {
        let cstate = CigarState {
            cig: r.cigar(),
            icig: CIGAR_STATE_UNINIT,
            iseq: 0,
            bam_pos: r.pos(),
            read_len_from_cigar,
        };

        let mut plp = PileupAlignment::new(r.clone(), cstate);

        // for the love of god, don't remove this.
        //
        // if we encounter a read that starts before pos, we need to initialize its cigar state to
        // catch up to the current position, as resolve_cigar increments in one-cigar-op intervals.
        //
        // if we don't, and if the query position is at cigar block 2+, the cigar state will stay
        // stale.
        if plp.rec.pos() < pos && plp.rec.tid() == self.tail.tid {
            for i in plp.rec.pos()..pos {
                resolve_cigar(&mut plp, i);
            }
        }

        let plp_ref = Rc::new(RefCell::new(plp));

        if let Some(overlap_map) = &mut self.overlap_map {
            overlap_map.push(Rc::clone(&plp_ref));
        }

        if self.depth == 0 {
            self.head.set_to_rec(r);
            self.tail.set_to_rec(r);
        }

        self.rbuf.push(Rc::clone(&plp_ref));
        self.depth += 1;
    }

    #[inline(always)]
    /// Attempt to add a read to the buffer for pending pileup generation. This is similar to
    /// bam_plp_push() from htslib.
    ///
    /// Reads can be added if and only if:
    /// 0. they are mapped
    /// 1. they overlap with the queried position
    /// 2. have a mapped mate
    /// 3. there is room in the buffer (dictated by Self::max_depth)
    ///
    /// Notably, reads found to be in a new reference are still added.
    pub fn attempt_push(&mut self, tid: i32, pos: i64, r: &Record) -> Result<BufPushResult, Error> {
        if r.is_unmapped() {
            if let Some(ref mut ov) = self.overlap_map {
                ov.delete_read(r);
            }

            return Ok(BufPushResult::Unmapped);
        }

        // check for unsorted input
        if r.tid() < self.tail.tid {
            error!(
                "File unsorted by reference: tid {} comes after tid {}. Read name: {}",
                r.tid(),
                self.tail.tid,
                std::str::from_utf8(r.qname()).unwrap()
            );

            anyhow::bail!("Unsorted");
        }

        if r.pos() < self.tail.pos && self.tail.tid == r.tid() {
            error!(
                "File unsorted by coordinate: pos {} comes after pos {}. Read name: {}",
                r.pos(),
                self.tail.pos,
                std::str::from_utf8(r.qname()).unwrap()
            );
            anyhow::bail!("Unsorted");
        }

        let read_len_from_cigar = cigar2rlen(r);

        if read_len_from_cigar > self.len {
            self.len = read_len_from_cigar;
        }

        if r.is_mate_unmapped() {
            if let Some(ref mut ov) = self.overlap_map {
                ov.delete_read(r);
            }
            return Ok(BufPushResult::BeforePos);
        }

        if r.tid() != tid {
            self.store(r, read_len_from_cigar, pos);
            return Ok(BufPushResult::DifferentReference);
        }

        if r.pos() <= pos && self.depth >= self.max_depth {
            if let Some(ov) = &mut self.overlap_map {
                ov.delete_read(r);
            }
            return Ok(BufPushResult::MaxDepthMet);
        }

        if r.pos() + read_len_from_cigar - 1 < pos {
            return Ok(BufPushResult::BeforePos);
        }

        self.store(r, read_len_from_cigar, pos);
        Ok(BufPushResult::Pushed)
    }

    /// Create a new read buffer, with optional overlap map if disable_overlaps == false
    pub fn new(depth: usize, disable_overlaps: bool) -> Self {
        let rbuf: Vec<PileupAlignmentRef> = Vec::with_capacity(500);
        let backup_buf: Vec<PileupAlignmentRef> = Vec::with_capacity(500);

        let max_depth = if depth.cmp(&0).is_eq() { usize::MAX } else { depth };
        let len = 0;

        let head = BufferBoundary::new();
        let tail = BufferBoundary::new();

        let overlap_map = match disable_overlaps {
            false => Some(HashMap::new()),
            true => None,
        };

        Self {
            rbuf,
            backup_buf,
            overlap_map,
            len,
            depth: 0,
            max_depth,
            head,
            tail,
        }
    }

    /// Prepare the buffer for the next round of pileup generation.
    pub fn reset(&mut self) {
        assert!(self.rbuf.is_empty());
        std::mem::swap(&mut self.rbuf, &mut self.backup_buf);
        if self.depth > 0 {
            self.head.set_to_rec(&self.rbuf[0].borrow().rec);
            self.tail.set_to_rec(&self.rbuf.last().unwrap().borrow().rec);
        } else {
            self.head = BufferBoundary::new();
            self.tail = BufferBoundary::new();
        }
    }
}
