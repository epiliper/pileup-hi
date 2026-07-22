use crate::{
    bamio::{BamDataSource, BamReader},
    baq::realign_record,
    cigar_resolve::resolve_cigar,
    engine::MIN_BAM_READ_THREADS,
    errors::{Error, ErrorKind},
    output::{OrderedPileupOutput, PileupOutputContext},
    params::PileupParams,
    position_queue::GenomeInterval,
    read_buf::{BufPushResult, ReadBuffer, ReadBufferEntry},
    read_filter::ReadFilter,
    refseq::RefSeqHandle,
    utils::read_ends_before_pos,
};

use rust_htslib::bam::Record;

#[derive(Clone)]
enum EmitStrategy {
    /// output absolutely nothing
    #[allow(dead_code)]
    Nothing,

    /// output for a position if it has coverage
    ByPos,

    /// output for an entire reference if it has coverage anywhere
    ByRef,

    /// output all coordinates for all references regardless of coverage
    Everything,
}

/// The state machine performing the pileup generation and advancing across region coordinates.
/// The type T dictates what output it generates.
pub struct PileupIteratorCore<T: OrderedPileupOutput> {
    tid: i32,
    next_tid: i32,
    last_tid_with_cov: i32,
    pos: i64,
    next_pos: i64,
    max_pos: i64,

    emit: EmitStrategy,

    rbuf: ReadBuffer,

    dest: T,
    emit_current: bool,

    pub reader: BamReader,
    refseq: RefSeqHandle,
    read_filter: ReadFilter,
    cur_rec: Record,
    realign: bool,
    min_baseq: u8,
    min_mapq: u8,
    redo_baq: bool,

    read_len: usize,
}

impl<T: OrderedPileupOutput> PileupIteratorCore<T> {
    /// Create a new pileup iterator from a data source (e.g. bam file), a set of query regions,
    /// input params and an output type.
    pub fn new(src: &BamDataSource, refseq: RefSeqHandle, params: &PileupParams) -> Result<Self, Error> {
        let read_len = BamReader::sample_read_len(src)?;

        let reader = BamReader::new(src, MIN_BAM_READ_THREADS)?;

        let rbuf = ReadBuffer::new(params.depth, params.disable_overlaps);

        let read_filter = ReadFilter::new(params.count_orphans, params.excl_flags, params.incl_flags);

        let cur_rec = Record::new();

        let emit = if params.show_empty_regions {
            EmitStrategy::Everything
        } else if params.show_empty_coords {
            EmitStrategy::ByRef
        } else {
            EmitStrategy::ByPos
        };

        let min_baseq = params.min_baseq;
        let min_mapq = params.min_mapq;

        let pos @ next_pos @ max_pos = -1;
        let tid @ next_tid @ last_tid_with_cov = -1;

        let dest = T::new();

        Ok(Self {
            tid,
            next_tid,
            last_tid_with_cov,
            pos,
            next_pos,
            max_pos,
            rbuf,
            dest,
            emit_current: false,
            emit,
            reader,
            read_filter,
            refseq,
            cur_rec,
            min_baseq,
            min_mapq,
            realign: !params.no_baq,
            redo_baq: params.redo_baq,
            read_len,
        })
    }

    #[inline(always)]
    pub fn is_exhausted(&self) -> bool {
        self.tid == -1
    }

    #[inline(always)]
    pub fn set_exhausted(&mut self) {
        self.tid = -1;
    }

    /// Generate a pileup from all bases passing the minimum quality filter and covering the
    /// iterator's current reference position. Importantly, generate_pileup() is where stale reads no longer
    /// overlapping the query position are removed.
    #[inline(always)]
    pub fn set_pileup(&mut self) -> Result<(), Error> {
        let mut skip = false;

        // don't bother going through read buffer if it starts beyond the
        // current coordinate
        if matches!(self.emit, EmitStrategy::ByPos) {
            if self.rbuf.head.tid == self.tid && self.rbuf.head.pos > self.pos {
                skip = true;
            }

            if self.rbuf.head.tid > self.tid {
                skip = true
            }
        }

        let rbuf = &mut self.rbuf;
        let generated;

        {
            self.dest.clear();

            if !skip {
                let ctx = PileupOutputContext {
                    tid: self.tid,
                    pos: self.pos,
                    ref_name: &self.reader.cur_ref,
                    refseq: &self.refseq,
                };

                generated = generate_pileup(rbuf, &mut self.dest, self.pos, self.tid, self.min_baseq, &ctx)?;
            } else {
                generated = false;
            };
        }

        //////////////////////////
        self.pos += 1;
        /////////////////////////

        if generated || self.dest.depth() > 0 {
            self.last_tid_with_cov = self.tid;
        }

        match self.emit {
            EmitStrategy::Nothing => {
                self.dest.clear();
                self.emit_current = false;
            }

            EmitStrategy::ByPos => {
                self.emit_current = generated || self.dest.depth() > 0;
                if !self.emit_current {
                    self.dest.clear();
                }
            }

            EmitStrategy::ByRef => {
                self.emit_current = self.tid == self.last_tid_with_cov || self.rbuf.head.tid == self.tid;
            }

            EmitStrategy::Everything => self.emit_current = true,
        }

        Ok(())
    }

    /// When given a region not starting at zero, rewind by 2X read length in order to populate the
    /// overlap map to ensure read mates get nulled.
    fn preload_region(&mut self, interval: &GenomeInterval) -> Result<(), Error> {
        let rewind = (interval.start - (2 * self.read_len) as i64).max(0);

        self.pos = rewind;
        self.next_pos = self.pos;
        self.max_pos = interval.start - 1;

        self.tid = interval.tid as i32;
        self.next_tid = self.tid;

        self.reader.init_to_ref(interval.tid as u32, self.pos, interval.end)?;

        let preset = self.emit.clone();
        self.emit = EmitStrategy::Nothing;

        while self.step()?.is_some() {
            continue;
        }

        self.pos = interval.start;
        self.next_pos = self.pos;
        self.next_tid = self.tid;

        self.emit = preset;

        Ok(())
    }

    /// Update iterator state and prepare ref-specific data given a new interval.
    pub fn set_ref(&mut self, interval: &GenomeInterval) -> Result<(), Error> {
        let tid = self.reader.name_to_tid(&interval.name)?.ok_or_else(|| {
            Error::from(ErrorKind::AnomalousData(format!(
                "Interval has a reference {} not in header",
                interval.name
            )))
        })?;

        let ctx = PileupOutputContext {
            tid: self.tid,
            pos: self.pos,
            ref_name: &self.reader.cur_ref,
            refseq: &self.refseq,
        };

        // purge read buffer to remove any reads spanning the old ref to update head and tail.
        generate_pileup(&mut self.rbuf, &mut self.dest, i64::MAX, self.tid, self.min_baseq, &ctx)?;

        self.dest.clear();

        if interval.start != 0 && self.rbuf.overlap_map.is_some() {
            self.preload_region(interval)?;
        } else {
            self.reader.init_to_ref(tid as u32, interval.start, interval.end)?;
            self.pos = interval.start;
            self.next_pos = interval.start;
        }

        self.tid = tid;
        self.next_tid = self.tid;

        self.max_pos = interval.end - 1;

        Ok(())
    }

    /// Read from an input BAM until we find either 1) read starting at a coordinate beyond the end of the
    /// read buffer, 2) run out of data, 3) exceed the max position of this region, or 4) encounter
    /// a read mapping to a new reference.
    #[inline(always)]
    pub fn intake(&mut self) -> Result<IterResult, Error> {
        if self.reader.eof {
            self.next_tid = -1;
            return Ok(IterResult::ReferenceEnd);
        }

        loop {
            // we need to keep reading until we have gathered all reads overlapping a position.
            // TODO: move the IO reading logic outside
            if let Some(read) = self.reader.read_no_alloc(&mut self.cur_rec) {
                read?;
                let r = &mut self.cur_rec;

                if r.is_unmapped() {
                    continue;
                }

                // Unindexed readers start at the beginning of the file for each query.
                if r.tid() < self.tid {
                    continue;
                }

                if !self.read_filter.check_read(r) {
                    continue;
                }

                if self.realign {
                    if let Some(refseq) = self.refseq.as_ref() {
                        let flag = if self.redo_baq { 7 } else { 3 };
                        realign_record(r, refseq, refseq.len() as i64, flag)?;
                    }
                }

                // we passed queried region
                if r.pos() > self.max_pos {
                    self.rbuf.attempt_push(self.tid, self.pos, r)?;
                    self.next_pos = r.pos();
                    return Ok(IterResult::ReferenceEnd);
                }

                if r.mapq() < self.min_mapq {
                    continue;
                }

                let ret = self.rbuf.attempt_push(self.tid, self.pos, r)?;

                match ret {
                    BufPushResult::Unmapped => (),

                    BufPushResult::DifferentReference => {
                        self.next_tid = r.tid();
                        self.next_pos = r.pos();
                        return Ok(IterResult::ReferenceEnd);
                    }

                    BufPushResult::Pushed => {
                        self.next_pos = r.pos();

                        if self.next_pos > self.pos {
                            return Ok(IterResult::Generated);
                        }
                    }

                    // if we've capped our buffer to a given depth, we'll iterate over all
                    // remaining reads spanning this coordinate before stopping to generate
                    // pileups. This way we won't have to deal with them at the next position.
                    BufPushResult::MaxDepthMet | BufPushResult::BeforePos => {
                        continue;
                    }
                }
            } else {
                // we ran out of reads.
                self.reader.eof = true;
                self.next_tid = -1;
                return Ok(IterResult::ReferenceEnd);
            }
        }
    }

    pub fn step(&mut self) -> Result<Option<()>, Error> {
        loop {
            if self.pos > self.max_pos {
                self.set_exhausted();
                return Ok(None);
            }

            // we have already hit the end of the current reference.
            if self.next_tid != self.tid {
                // still have some positions left
                if self.rbuf.head.tid == self.tid
                    || matches!(self.emit, EmitStrategy::ByRef)
                    || matches!(self.emit, EmitStrategy::Everything)
                {
                    return Some(self.set_pileup()).transpose();
                } else {
                    self.set_exhausted();
                    return Ok(None); // we are done
                }
            }

            // we have reads with coordinates in the buffer that we need to process
            if self.rbuf.head.pos < self.pos
                || (matches!(self.emit, EmitStrategy::ByRef) && self.rbuf.head.tid == self.tid)
                || matches!(self.emit, EmitStrategy::Everything)
            {
                if self.pos < self.next_pos {
                    return Some(self.set_pileup()).transpose();
                }
                // or we don't, in which case we just jump ahead.
            } else {
                self.pos = self.rbuf.head.pos;

                if self.pos < self.next_pos {
                    return Some(self.set_pileup()).transpose();
                }
            }

            // we still need to sample more reads at this position
            if let Err(e) = self.intake() {
                self.set_exhausted();
                return Some(Err(e)).transpose(); // propagate error if something goes wrong here.
            }
        }
    }
}

pub enum IterResult {
    ReferenceEnd,
    Generated,
}

#[inline(always)]
/// Perform the pileup given a read buffer, optional ref sequence, an output destination, and query
/// position. Importantly, reads found to no longer overlap (pos, tid) will be removed from the
/// buffer.
fn generate_pileup<T: OrderedPileupOutput>(
    rbuf: &mut ReadBuffer,
    out: &mut T,
    pos: i64,
    tid: i32,
    min_baseq: u8,
    ctx: &PileupOutputContext,
) -> Result<bool, Error> {
    let mut generated = false;
    let mut plp;

    for (idx, raw) in rbuf.rbuf.iter().enumerate() {
        match raw {
            ReadBufferEntry::Tombstone => continue,
            ReadBufferEntry::Occupied(_plp) => plp = _plp,
        }

        unsafe {
            let r = plp.get();

            // record starts beyond position, which means that the remainder of the buffer does
            // too. Skip the rest of the records.
            if r.rec.tid() > tid || (r.rec.pos() > pos && r.rec.tid() == tid) {
                // drop(r);
                break;
            }

            // record is old and no longer overlaps the query coordinate. We discard it by not adding
            // it to the alternate buffer.
            if read_ends_before_pos(r, pos) || r.rec.tid() < tid {
                rbuf.remove(idx);
                continue;
            }

            generated = true;

            // advance to the current ref position in read and record cigar op
            resolve_cigar(r, pos);
            let qual = *r.rec.qual().get(r.qpos).unwrap_or(&0);

            if qual < min_baseq {
                continue;
            }

            out.intake(ctx, r)?;
        }
    }

    rbuf.reset()?;
    Ok(generated)
}

pub struct PileupIterator<T: OrderedPileupOutput> {
    core: Vec<PileupIteratorCore<T>>,
    tid: i32,
    pos: i64,
    output_tid: i32,
    output_pos: i64,
    pub max_pos: i64,
    refseq: RefSeqHandle,
    ref_name: String,
}

impl<T: OrderedPileupOutput> PileupIterator<T> {
    pub fn advance(&mut self) -> Result<Option<()>, Error> {
        loop {
            if self.pos > self.max_pos {
                return Ok(None);
            }

            let output_tid = self.tid;
            let output_pos = self.pos;
            let mut next_pos: i64 = i64::MAX;

            for iter in self.core.iter_mut() {
                if iter.is_exhausted() {
                    continue;
                }

                while !iter.is_exhausted() && iter.pos <= self.pos {
                    iter.step()?;
                }

                next_pos = std::cmp::min(next_pos, iter.pos);
            }

            if next_pos == i64::MAX {
                return Ok(None);
            }

            self.pos = std::cmp::min(next_pos, self.pos + 1);

            let has_output = self
                .core
                .iter()
                .any(|core| !core.is_exhausted() && core.pos == self.pos && core.emit_current);

            if has_output {
                self.output_tid = output_tid;
                self.output_pos = output_pos;
                return Ok(Some(()));
            }
        }
    }

    pub fn ctx(&self) -> PileupOutputContext<'_> {
        PileupOutputContext {
            tid: self.output_tid,
            pos: self.output_pos,
            ref_name: &self.ref_name,
            refseq: &self.refseq,
        }
    }

    pub fn current(&self) -> impl Iterator<Item = Option<&T>> {
        self.core.iter().map(|core| {
            if !core.is_exhausted() && core.pos == self.pos && core.emit_current {
                Some(&core.dest)
            } else {
                None
            }
        })
    }

    fn set_ref(&mut self, interval: &GenomeInterval) -> Result<(), Error> {
        self.tid = interval.tid as i32;
        self.pos = interval.start;
        self.output_tid = -1;
        self.output_pos = -1;
        self.max_pos = interval.end - 1;

        let mut valid = false;

        for iter in self.core.iter_mut() {
            // if we have a sample that doesn't contain the interval we want, just don't run it.
            if iter.set_ref(interval).is_err() {
                iter.set_exhausted();
            } else {
                valid = true;
            }
        }

        // if no samples have the interval we asked for, then we consider this an anomaly and return
        // an error.
        if !valid {
            return Err(Error::from(ErrorKind::AnomalousData(format!(
                "Interval {:?} is not found in any input samples.",
                interval
            ))));
        }

        Ok(())
    }

    pub fn from_query(
        src: &[BamDataSource],
        refseq: RefSeqHandle,
        interval: &GenomeInterval,
        params: &PileupParams,
    ) -> Result<Self, Error> {
        let mut iterators = Vec::with_capacity(src.len());

        for src in src {
            iterators.push(PileupIteratorCore::<T>::new(src, refseq.clone(), params)?);
        }

        let ref_name = interval.name.clone();
        let refseq = refseq.clone();

        let mut s = Self {
            core: iterators,
            ref_name,
            refseq,
            tid: -1,
            pos: -1,
            output_tid: -1,
            output_pos: -1,
            max_pos: -1,
        };

        s.set_ref(interval)?;
        Ok(s)
    }
}
