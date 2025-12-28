use anyhow::Error;
use rust_htslib::{bam::Record, htslib::sam_prob_realn};

const REALN_FAIL: i32 = 4;

// Method to realign record using htslib's internal realignment tools.
pub fn realign_record(r: &mut Record, ref_seq: &[u8], ref_len: i64, flag: u32) -> Result<i32, Error> {
    unsafe {
        let ptr: *const i8 = ref_seq.as_ptr() as *const i8;
        let ret = sam_prob_realn(r.inner_mut(), ptr, ref_len, flag as std::os::raw::c_int); // temp value for testing
        if ret == REALN_FAIL {
            anyhow::bail!("realignment with sam_prob_realn failed")
        } else {
            Ok(ret)
        }
    }
}
