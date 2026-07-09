use crate::alignment::{cigar2rlen, CigarState, PileupAlignment, PileupAlignmentRef, CIGAR_STATE_UNINIT};
use crate::cigar_resolve::resolve_cigar;
use crate::errors::{Error, ErrorKind};
use crate::overlap::{MapOverlaps, OverlapMap};
use likely_stable::unlikely;
use rust_htslib::bam::Record;
use std::{cell::RefCell, collections::HashMap};

#[allow(clippy::large_enum_variant)]
pub enum ReadBufferEntry {
    Occupied(PileupAlignmentRef),
    Tombstone,
}

pub struct ReadBuffer {
    pub rbuf: Vec<ReadBufferEntry>,
    pub pending_del_indexes: RefCell<Vec<usize>>,
    pub n_tombstones: usize,
    pub overlap_map: Option<OverlapMap>,
    pub depth: usize,
    pub max_depth: usize,
    pub head: BufferBoundary,
    pub tail: BufferBoundary,
    pub len: i64,
}

pub struct BufferBoundary {
    pub tid: i32,
    pub pos: i64,
}

fn rbuf_get(rbuf: &mut [ReadBufferEntry], key: usize) -> Result<&PileupAlignmentRef, Error> {
    if let ReadBufferEntry::Occupied(t) = rbuf
        .get(key)
        .ok_or(Error::from(ErrorKind::ThreadingError("Invalid read buffer index")))?
    {
        Ok(t)
    } else {
        Err(Error::from(ErrorKind::ThreadingError(
            "Invalid read buffer index: hit tombstone",
        )))
    }
}

#[inline(always)]
fn rbuf_first(rbuf: &[ReadBufferEntry]) -> Option<&PileupAlignmentRef> {
    for entry in rbuf {
        match entry {
            ReadBufferEntry::Tombstone => continue,

            ReadBufferEntry::Occupied(plp) => {
                return Some(plp);
            }
        }
    }

    None
}

#[inline(always)]
fn rbuf_last(rbuf: &[ReadBufferEntry]) -> Option<&PileupAlignmentRef> {
    for entry in rbuf {
        match entry {
            ReadBufferEntry::Tombstone => continue,

            ReadBufferEntry::Occupied(plp) => {
                return Some(plp);
            }
        }
    }

    None
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

pub fn handle_overlaps(
    ov: &mut OverlapMap,
    readbuf: &mut [ReadBufferEntry],
    idx: usize,
    plp: &PileupAlignmentRef,
) -> Result<(), Error> {
    if let Some(mate) = ov.lookup(plp) {
        let other = rbuf_get(readbuf, mate)?;
        ov.nullify(plp, other);
        // otherwise store
    } else {
        ov.store(idx, plp);
    }

    Ok(())
}

impl ReadBuffer {
    #[inline(always)]
    /// Store the read we've decided to push in the buffer, wrapping it in an Rc pointer and
    /// performing mate overlap-correction if requested by caller.
    pub fn store(&mut self, r: &Record, read_len_from_cigar: i64, pos: i64) -> Result<(), Error> {
        let cstate = CigarState {
            cig: r.cigar(),
            icig: CIGAR_STATE_UNINIT,
            iseq: 0,
            bam_pos: r.pos(),
            read_len_from_cigar,
        };

        let mut plp = PileupAlignment::new(r.clone(), cstate);

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

        let plp_ref = PileupAlignmentRef::new(plp);
        let idx = self.rbuf.len();

        if let Some(overlap_map) = &mut self.overlap_map {
            handle_overlaps(overlap_map, &mut self.rbuf, idx, &plp_ref)?;
        }

        if self.depth == 0 {
            self.head.set_to_rec(r);
            self.tail.set_to_rec(r);
        }

        self.rbuf.push(ReadBufferEntry::Occupied(plp_ref));
        self.depth += 1;
        Ok(())
    }

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
    #[inline(always)]
    pub fn attempt_push(&mut self, tid: i32, pos: i64, r: &Record) -> Result<BufPushResult, Error> {
        if r.is_unmapped() {
            if let Some(ref mut ov) = self.overlap_map {
                ov.delete_read(r);
            }

            return Ok(BufPushResult::Unmapped);
        }

        // check for unsorted input
        if unlikely(r.tid() < self.tail.tid) {
            return Err(Error::from(ErrorKind::BamNotSortedByReference(self.tail.tid, r.tid())));
        }

        if unlikely(r.pos() < self.tail.pos) {
            return Err(Error::from(ErrorKind::BamNotSortedByCoordinate(self.tail.pos, r.pos())));
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
            self.store(r, read_len_from_cigar, pos)?;
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

        self.store(r, read_len_from_cigar, pos)?;
        Ok(BufPushResult::Pushed)
    }

    /// Create a new read buffer, with optional overlap map if disable_overlaps == false
    pub fn new(depth: usize, disable_overlaps: bool) -> Self {
        let rbuf: Vec<ReadBufferEntry> = Vec::with_capacity(500);

        let max_depth = if depth.cmp(&0).is_eq() { usize::MAX } else { depth };
        let len = 0;

        let head = BufferBoundary::new();
        let tail = BufferBoundary::new();

        let overlap_map = match disable_overlaps {
            false => Some(HashMap::new()),
            true => None,
        };

        let pending_del_indexes = RefCell::new(Vec::with_capacity(500));

        Self {
            rbuf,
            overlap_map,
            pending_del_indexes,
            len,
            depth: 0,
            n_tombstones: 0,
            max_depth,
            head,
            tail,
        }
    }

    fn remove_tombstones(&mut self) -> Result<(), Error> {
        let mut dest = Vec::with_capacity(self.rbuf.capacity());
        // if let Some(ref mut ov) = self.overlap_map {
        //     ov.clear();
        // }

        for entry in self.rbuf.drain(..) {
            match entry {
                ReadBufferEntry::Occupied(plp) => {
                    if let Some(ref mut ov) = self.overlap_map {
                        ov.store(dest.len(), &plp);
                        // handle_overlaps(ov, &mut self.rbuf, idx, &plp)?;
                    }
                    dest.push(ReadBufferEntry::Occupied(plp));
                }
                ReadBufferEntry::Tombstone => (),
            }
        }

        self.n_tombstones = 0;
        self.rbuf = dest;
        Ok(())
    }

    pub fn remove(&self, index: usize) {
        self.pending_del_indexes.borrow_mut().push(index);
    }

    pub fn reset(&mut self) -> Result<(), Error> {
        {
            let mut indexes = self.pending_del_indexes.borrow_mut();

            for idx in indexes.drain(..) {
                let entry = self.rbuf.get_mut(idx).expect("Bad buffer removal index");
                match entry {
                    ReadBufferEntry::Tombstone => {
                        panic!("Attempted to delete tombstone in buffer")
                    }
                    ReadBufferEntry::Occupied(plp) => {
                        self.n_tombstones += 1;
                        self.depth -= 1;

                        unsafe {
                            if let Some(ref mut overlap) = self.overlap_map {
                                overlap.delete_read(&plp.get().rec);
                            }
                        }

                        *entry = ReadBufferEntry::Tombstone;
                    }
                }
            }

            assert!(indexes.is_empty());
        }

        if self.n_tombstones >= ((self.rbuf.len() as f32) * 0.77) as usize {
            self.remove_tombstones()?;
            assert!(self.n_tombstones == 0);
        }

        unsafe {
            if self.depth > 0 {
                let first = rbuf_first(&self.rbuf).unwrap();
                self.head.tid = first.get().rec.tid();
                self.head.pos = first.get().rec.pos();

                let last = rbuf_last(&self.rbuf).unwrap();
                self.tail.tid = last.get().rec.tid();
                self.tail.pos = last.get().rec.pos();
            } else {
                self.head = BufferBoundary::new();
                self.tail = BufferBoundary::new();
            }
        }

        Ok(())
    }
}
