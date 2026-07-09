use crate::alignment::PileupAlignmentRef;
use crate::errors::{Error, ErrorKind};
use rust_htslib::bam::Record;
use std::collections::HashMap;
use std::hash::{DefaultHasher, Hash, Hasher};

extern "C" {
    fn tweak_overlap_quality(a: *mut rust_htslib::htslib::bam1_t, b: *mut rust_htslib::htslib::bam1_t) -> i32;
}

pub fn tweak_overlap_qual(a: &mut Record, b: &mut Record) -> Result<(), Error> {
    unsafe {
        let ret = tweak_overlap_quality(a.inner_mut() as *mut _, b.inner_mut() as *mut _);
        if ret < 0 {
            let qname = std::str::from_utf8(a.qname())?;
            return Err(Error::from(ErrorKind::MateOverlapFailed(qname.to_string())));
        }
    }
    Ok(())
}

pub type OverlapMap = HashMap<u64, usize>;

pub trait MapOverlaps {
    fn lookup(&self, rec: &PileupAlignmentRef) -> Option<usize>;
    fn store(&mut self, val: usize, rec: &PileupAlignmentRef);
    fn nullify(&mut self, rec1: &PileupAlignmentRef, rec2: &PileupAlignmentRef);
    fn delete_hash(&mut self, r: u64);
    fn delete_read(&mut self, r: &Record);
}

pub fn hash_qname(r: &Record) -> u64 {
    let mut hasher = DefaultHasher::new();
    r.qname().hash(&mut hasher);
    hasher.finish()
}

impl MapOverlaps for OverlapMap {
    fn lookup(&self, rec: &PileupAlignmentRef) -> Option<usize> {
        unsafe {
            let _r = rec.get();
            let len = _r.cstate.read_len_from_cigar;
            let r = &mut _r.rec;

            if r.is_mate_unmapped() || !r.is_proper_pair() {
                return None;
            }

            if (r.mtid() >= 0 && (r.mtid() != r.tid()))
                || r.insert_size().abs() >= 2 * (r.seq_len() as i64) && r.mpos() >= r.pos() + len
            {
                return None;
            }

            let h = hash_qname(r);

            self.get(&h).copied()
        }
    }

    fn store(&mut self, val: usize, rec: &PileupAlignmentRef) {
        unsafe {
            let h = hash_qname(&rec.get().rec);
            self.insert(h, val);
        }
    }

    fn nullify(&mut self, cur_rec: &PileupAlignmentRef, prev_rec: &PileupAlignmentRef) {
        unsafe {
            let h = hash_qname(&cur_rec.get().rec);
            tweak_overlap_qual(&mut prev_rec.get().rec, &mut cur_rec.get().rec).unwrap();
            self.delete_hash(h)
        }
    }

    fn delete_hash(&mut self, r: u64) {
        self.remove(&r);
    }

    fn delete_read(&mut self, r: &Record) {
        let mut h = DefaultHasher::new();
        r.qname().hash(&mut h);
        self.remove(&h.finish());
    }
}
