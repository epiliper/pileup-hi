use anyhow::Error;
use rust_htslib::{
    bam::Record,
    htslib::{htsRealnFlags, sam_prob_realn},
};

const REALN_FAIL: i32 = 4;

// Method to realign record using htslib's internal realignment tools.
pub fn realign_record(r: &mut Record, ref_seq: &[u8], ref_len: i64, _flag: htsRealnFlags) -> Result<i32, Error> {
    unsafe {
        let ptr: *const i8 = ref_seq.as_ptr() as *const i8;
        let ret = sam_prob_realn(r.inner_mut(), ptr, ref_len, 3); // temp value for testing
        if ret == REALN_FAIL {
            anyhow::bail!("realignment with sam_prob_realn failed")
        } else {
            Ok(ret)
        }
    }
}
