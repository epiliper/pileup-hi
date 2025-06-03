use crate::params::Params;
use crate::pileup::PileUp;
use crate::plp_writers::*;
use crate::refseq::RefSeq;
use crate::rpileup::*;
use crate::{overlap::OverlapMap, read_filter::ReadFilter};
use crate::{
    read_buf::{BufPushResult, ReadBuffer},
    rpileup::IterResult,
};
use anyhow::{Context, Error};
use rust_htslib::bam::{
    ext::BamRecordExtensions, record::Cigar, HeaderView, IndexedReader, Read, Record,
};
use std::borrow::Borrow;
use std::cell::RefCell;
use std::rc::Rc;

const UNINIT_POS: usize = usize::MAX - 1;
const UNINIT_TID: u32 = u32::MAX - 1;

pub struct PileupIterator {
    overlaps: Option<OverlapMap>,
    rbuf: ReadBuffer,
    pub iter: Vec<(u32, Option<Rc<RefCell<PileUp>>>)>, // note u32's innate limit to bam size
    reader: IndexedReader,
    cur_rec: Record,
    pos: usize,
    next_pos: usize,
    header: HeaderView,
    tid: u32,
    read_filter: ReadFilter,
    refseq: Option<RefSeq>,
    max_pos: usize,
    plp_counter: u32,
    seq_buf: Vec<u8>,
    qual_buf: Vec<u8>,
    min_baseq: u8,
    show_all: bool,
}

impl PileupIterator {
    pub fn new(params: Params) -> Result<Self, Error> {
        let tid = params.inp.tid.unwrap_or(UNINIT_TID);
        let pos @ next_pos @ max_pos = params.inp.pos.unwrap_or(UNINIT_POS);
        let mut reader = IndexedReader::from_path(params.inp.input)?;
        reader.set_threads(num_cpus::get())?;
        let rbuf = ReadBuffer::new();
        let header = reader.header().clone();
        let show_all = params.plp.show_empty_coords;
        let (seq_buf, qual_buf) = (Vec::with_capacity(500), Vec::with_capacity(500));
        let cur_rec = Record::new();
        let mut refseq = None;
        let min_baseq = params.plp.min_baseq;
        let mut iter = Vec::new();
        let mut plp_counter = 0;

        let read_filter = ReadFilter::new(
            params.plp.min_mapq,
            params.plp.count_orphans,
            params.plp.excl_flags.iter().map(|s| s.as_str()).collect(),
            params.plp.incl_flags.iter().map(|s| s.as_str()).collect(),
        )?;

        if let Some(ref_file) = params.inp.refseq {
            refseq = Some(RefSeq::from_file(ref_file)?);
        }

        let overlaps = match params.plp.disable_overlap {
            true => None,
            false => Some(OverlapMap::new()),
        };

        Ok(Self {
            tid,
            pos,
            next_pos,
            iter,
            plp_counter,
            max_pos,
            rbuf,
            reader,
            overlaps,
            header,
            min_baseq,
            read_filter,
            show_all,
            refseq,
            seq_buf,
            qual_buf,
            cur_rec,
        })
    }

    pub fn next(&mut self) -> Result<IterResult, Error> {
        if self.pos >= self.max_pos {
            return Ok(IterResult::ReferenceEnd);
        }

        let mut gen = false;

        // if we are at the next position in the bam where reads are within window range,
        // resume read intake
        if self.pos == self.next_pos {
            self.fill_buffer()?;
        }

        // if we have reads in buffer, attempt to generate plp.
        if !self.rbuf.rbuf.is_empty() {
            gen = self.set_pileup()?;
        }

        // if no reads in buffer overlapped with pos, print empty plp if enabled
        if !gen && self.show_all {
            self.write_pileup_str(b'N', 0, 0, 0)?;
        }

        // if we need to print blank plps for each col,
        // advance query coord by 1
        // else, jump to the next coord with reads in range
        if self.show_all || !self.rbuf.rbuf.is_empty() {
            self.pos += 1;
        } else {
            self.pos = self.next_pos;
        }

        return Ok(IterResult::Generated);
    }

    pub fn init_to_ref(&mut self) -> Result<IterResult, Error> {
        // todo: check if this works for bam files without refs in header
        //
        if self.tid == UNINIT_TID {
            self.tid = 0;
        } else {
            self.tid += 1;
        }

        if self.tid >= self.header.target_count() {
            Ok(IterResult::NoData)
        } else {
            if let Some(r) = self.refseq.as_mut() {
                // right now we just get the entire reference sequence.
                // Next step will be to load it in windows.
                let tidname = std::str::from_utf8(self.header.tid2name(self.tid));
                r.load_seq(
                    tidname?,
                    0,
                    self.header
                        .target_len(self.tid)
                        .context("Failed to get target length")?,
                )?
            }

            self.max_pos = self.header.target_len(self.tid).context("No ref len")? as usize;
            self.pos = 0;
            self.next_pos = 0;
            self.reader.fetch((self.tid, 0, u32::MAX))?;
            Ok(IterResult::Generated)
        }
    }

    pub fn register_plp(&mut self, plp_ref: Rc<RefCell<PileUp>>) {
        self.plp_counter += 1;
        self.iter.push((self.plp_counter, Some(plp_ref)));
    }

    pub fn prune(&mut self, v: Vec<usize>) {
        for i in v {
            self.iter[i] = (0, None)
        }
        let mut i = self.iter.len() - 1;
        loop {
            let (_idx, node) = &self.iter[i];
            if node.is_none() {
                self.iter.remove(i);
                i -= 1
            } else {
                break;
            }
        }
    }

    pub fn set_pileup(&mut self) -> Result<bool, Error> {
        let mut generated = false;
        let mut to_remove: Vec<usize> = Vec::new();

        let mut ndel @ mut nins @ mut nbases = 0;
        let ref_base = match &self.refseq {
            Some(seq) => seq.get_base(self.pos as u64)?,
            None => b'N',
        };

        for (idx, (i, raw)) in self.iter.iter_mut().enumerate() {
            if let Some(ra) = raw {
                let mut r = ra.borrow_mut();
                if r.rec.reference_end() - 1 < self.pos as i64 {
                    to_remove.push(idx);
                    // self.iter.remove(idx);

                    continue;
                }

                let mut ipos: i32 = -1;
                let ret = cigar_get_pos(&mut r.cstate, self.pos as u32, &mut ipos);

                if ipos != -1 && r.rec.qual()[ipos as usize] < self.min_baseq {
                    continue;
                }

                match ret {
                    CigarAtPos::Op(Cigar::Match(_)) => {
                        write_match(
                            &r.cstate,
                            &r.rec,
                            ipos as u32,
                            self.pos,
                            &mut self.seq_buf,
                            &mut self.qual_buf,
                            self.refseq.as_ref(),
                        )?;

                        nbases += 1;
                    }

                    CigarAtPos::Op(Cigar::Ins(l)) => {
                        nins += 1;
                        r.indel += l;
                        write_ins(
                            &r.cstate,
                            &r.rec,
                            ipos as u32,
                            &mut self.seq_buf,
                            &mut self.qual_buf,
                        )?;
                    }

                    CigarAtPos::Op(Cigar::Del(l)) => {
                        if ipos != -1 {
                            write_del(self.pos, &mut self.seq_buf, l as usize)?;
                            r.indel -= l;
                        } else {
                            self.seq_buf.push(b'*');
                        }
                        ndel += 1;
                    }

                    CigarAtPos::BeforePos() => {
                        panic!(
                            "{} {} {}",
                            r.rec.is_unmapped(),
                            self.pos,
                            r.rec.reference_end() - 1
                        );
                    }

                    CigarAtPos::BaseEmpty() => (),
                    _ => panic!(),
                }
            }
        }

        if nbases + nins + ndel > 0 {
            self.write_pileup_str(ref_base, nbases, nins, ndel)?;
            generated = true;
        }

        self.prune(to_remove);

        Ok(generated)
    }

    pub fn write_pileup_str(
        &mut self,
        ref_base: u8,
        nbases: usize,
        nins: usize,
        ndel: usize,
    ) -> Result<(), Error> {
        print! {"{}\t{}\t{}\t{}\t", std::str::from_utf8(self.header.tid2name(self.tid))?, self.pos + 1, char::from(ref_base), nbases + nins + ndel }
        if self.seq_buf.is_empty() {
            print! {"*\t"}
        } else {
            print! {"{}\t", std::str::from_utf8(&self.seq_buf)?}
            self.seq_buf.clear();
        }

        if self.qual_buf.is_empty() {
            print! {"*\t"}
        } else {
            print! {"{}\t", std::str::from_utf8(&self.qual_buf)?}
            self.qual_buf.clear();
        }

        print! {"\n"}

        Ok(())
    }

    /// Read records in to fill a read buffer spanning the current coordinate window.
    /// This will loop over records until A) a record is found that starts outside the current
    /// window, (e.g. faraway coord or different reference).
    ///
    /// When a read outside the current window is found, the [PileupIterator] will skip buffer
    /// filling / pileup generation for all coordinates between the current and the next read's
    /// start position.
    pub fn fill_buffer(&mut self) -> Result<(), Error> {
        let mut ret: BufPushResult;
        // let mut scanned = 0;

        let mut prev_pos = i64::MIN;

        if self.cur_rec.tid() == -1 {
            if let Some(rec) = self.reader.read(&mut self.cur_rec) {
                rec?;
            } else {
                // if we have no reads at all to set next pos, assume
                // we've hit the end of the reference, and set next pos to MAX
                self.next_pos = usize::MAX;
                return Ok(());
            }
        }

        loop {
            if self.cur_rec.tid() == -1 {
                break;
            }

            let r = &self.cur_rec;

            if r.is_unmapped() {
                match self.reader.read(&mut self.cur_rec) {
                    None => break,
                    Some(_) => continue,
                };
            }

            if !self.read_filter.check_read(&r) {
                match self.reader.read(&mut self.cur_rec) {
                    None => break,
                    Some(_) => continue,
                };
            }

            if r.pos() < prev_pos {
                panic!("UNSORTED BAM! {} {}", r.pos(), prev_pos)
            }

            prev_pos = r.pos();

            ret = self.rbuf.attempt_push(&r, self.pos, self.tid);

            match ret {
                BufPushResult::Unmapped => panic!(),

                BufPushResult::AfterWindow(next_pos) => {
                    self.next_pos = self.pos + next_pos;
                    break;
                }

                BufPushResult::DifferentReference => {
                    break;
                }
                BufPushResult::Pushed(plp_ref) | BufPushResult::MateFound(plp_ref) => {
                    self.cur_rec.set_tid(-1);
                    self.register_plp(plp_ref);
                    match self.reader.read(&mut self.cur_rec) {
                        Some(Ok(_)) => continue,
                        None => break,
                        Some(Err(_)) => panic!(),
                    }
                }
            };
        }

        Ok(())
    }
}
