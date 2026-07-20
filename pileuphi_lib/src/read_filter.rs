use rust_htslib::bam::Record;

pub struct ReadFilter {
    inc_flag: u16,
    exc_flag: u16,
    count_orphans: bool,
}

impl ReadFilter {
    pub fn new(count_orphans: bool, excl_flags: u16, incl_flags: u16) -> Self {
        Self {
            inc_flag: incl_flags,
            exc_flag: excl_flags,
            count_orphans,
        }
    }

    #[inline(always)]
    pub fn check_read(&mut self, read: &Record) -> bool {
        // let mut pass;

        // check if orphan pair
        if !self.count_orphans && (read.is_paired() && !read.is_proper_pair()) {
            return false;
        }

        if self.inc_flag > 0 && (self.inc_flag & read.inner.core.flag) == 0 {
            return false;
        }

        if self.exc_flag > 0 && (self.exc_flag & read.inner.core.flag) > 0 {
            return false;
        }

        true
    }
}
