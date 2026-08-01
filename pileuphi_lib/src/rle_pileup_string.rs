#![allow(clippy::unused_io_amount)]
use crate::alignment::PileupAlignment;
use crate::errors::Error;
use crate::output::{OrderedPileupOutput, PileupOutputContext};
use crate::pileup_string::PileupString;

#[derive(Clone, Debug)]
pub struct RLEPileupString {
    core: PileupString,
}

#[inline(always)]
fn is_base(base: u8) -> bool {
    let b = base.to_ascii_uppercase();
    b == b'A' || b == b'G' || b == b'C' || b == b'T' || b == b'N'
}

#[inline(always)]
#[allow(dead_code)]
fn is_any(_base: u8) -> bool {
    true
}

fn format_int(buf: &[u8], ndigits: usize) -> u32 {
    unsafe {
        std::str::from_utf8_unchecked(&buf[0..ndigits])
            .parse::<u32>()
            .expect("Invalid integer parse")
    }
}

#[inline(always)]
fn rlebuf_write<W: std::io::Write>(
    inbuf: &[u8],
    writer: &mut W,
    formatbuf: &mut itoa::Buffer,
    char_should_end_substring: fn(u8) -> bool,
) -> Result<(), Error> {
    let mut last: [u8; 10] = [0; 10];
    let mut cur: [u8; 10] = [0; 10];
    let mut count = 1;
    let mut start = 0;

    let mut digits: [u8; 32] = [0; 32];
    let mut ndigits = 0;

    // tracks letters we have to iterate over to capture entire string, e.g.
    // -2AG means we have to iterate over A and G without terminating the string.
    let mut avoid_terminate: u32 = 0;

    // preload last buf
    for b in inbuf.iter().copied() {
        last[start] = b;
        start += 1;

        avoid_terminate = avoid_terminate.saturating_sub(1);

        if b.is_ascii_digit() {
            digits[ndigits] = b;
            ndigits += 1;
            continue;
        }

        if char_should_end_substring(b) && avoid_terminate == 0 {
            if ndigits > 0 {
                avoid_terminate = format_int(&digits, ndigits);
                ndigits = 0;
                continue;
            }

            break;
        }
    }

    let mut i = 0; // track position along last
    let mut lastlen = start; // track length of previous strmatch

    if start < inbuf.len() {
        for b in inbuf[start..].iter().copied() {
            cur[i] = b;

            avoid_terminate = avoid_terminate.saturating_sub(1);

            if b.is_ascii_digit() {
                digits[ndigits] = b;
                ndigits += 1;
                continue;
            }

            if b == last[i] {
                i += 1;

                if char_should_end_substring(b) && avoid_terminate == 0 {
                    count += 1;
                    i = 0;
                }
            } else {
                if char_should_end_substring(b) && avoid_terminate == 0 {
                    writer.write(formatbuf.format(count).as_bytes())?;
                    writer.write(&last[0..lastlen])?;

                    // reset unused chars, otherwise on a subsequent iteration a stale value in last (copied from cur) could result
                    // in a false match
                    for j in (i + 1)..5 {
                        cur[j] = 0;
                    }

                    last = cur;

                    count = 1;
                    lastlen = i + 1;
                    i = 0;
                } else {
                    i += 1;
                }
            }
        }
    }

    writer.write(formatbuf.format(count).as_bytes())?;
    writer.write(&last[0..lastlen])?;
    Ok(())
}

impl RLEPileupString {
    #[inline(always)]
    pub fn write_body_rle<W: std::io::Write>(&self, writer: &mut W) -> Result<(), Error> {
        let mut buf = itoa::Buffer::new();
        writer.write(buf.format(self.core.depth).as_bytes())?;
        writer.write(b"\t")?;

        if self.core.seqbuf.is_empty() {
            writer.write(b"*\t")?;
        } else {
            rlebuf_write(&self.core.seqbuf, writer, &mut buf, is_base)?;
            writer.write(b"\t")?;
        }

        if self.core.qualbuf.is_empty() {
            writer.write_all(b"*\t")?;
        } else {
            // rlebuf_write(&self.core.qualbuf, writer, &mut buf, is_any)?;
            writer.write(&self.core.qualbuf)?;
            writer.flush()?;
        }

        Ok(())
    }
}

impl OrderedPileupOutput for RLEPileupString {
    #[inline(always)]
    fn intake(&mut self, ctx: &PileupOutputContext, p: &PileupAlignment) -> Result<(), Error> {
        self.core.depth += 1;
        self.core.write_plp(p, ctx.pos, ctx.refseq)?;
        Ok(())
    }

    #[inline(always)]
    fn write_header<W: std::io::Write>(ctx: &PileupOutputContext, writer: &mut W) -> Result<(), Error> {
        PileupString::write_refinfo(ctx, writer)
    }

    #[inline(always)]
    fn write_body<W: std::io::Write>(&self, _ctx: &PileupOutputContext, writer: &mut W) -> Result<(), Error> {
        self.write_body_rle(writer)
    }

    fn write_body_empty<W: std::io::Write>(_ctx: &PileupOutputContext, writer: &mut W) -> Result<(), Error> {
        writer.write_all(b"0\t*\t*")?;
        Ok(())
    }

    #[inline(always)]
    fn depth(&self) -> u32 {
        self.core.depth()
    }

    #[inline(always)]
    fn new() -> Self {
        Self {
            core: PileupString::new(),
        }
    }

    #[inline(always)]
    fn clear(&mut self) {
        self.core.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rle_1() {
        let mut outbuf: Vec<u8> = vec![];
        let inbuf = "^[A^[A^[G".as_bytes().to_vec();
        let mut buf = itoa::Buffer::new();

        eprintln!("{}", std::str::from_utf8(&inbuf).unwrap());
        rlebuf_write(&inbuf, &mut outbuf, &mut buf, is_base).unwrap();
        assert_eq!(std::str::from_utf8(&outbuf).unwrap(), "2^[A1^[G");
    }

    #[test]
    fn test_rle_2() {
        let mut outbuf: Vec<u8> = vec![];
        let inbuf = "^[A^[A^[A[^G[^G[^G".as_bytes().to_vec();
        let mut buf = itoa::Buffer::new();

        eprintln!("{}", std::str::from_utf8(&inbuf).unwrap());
        rlebuf_write(&inbuf, &mut outbuf, &mut buf, is_base).unwrap();
        assert_eq!(std::str::from_utf8(&outbuf).unwrap(), "3^[A3[^G");
    }

    #[test]
    fn test_rle_3() {
        let mut outbuf: Vec<u8> = vec![];
        let inbuf = "AAAAGGGG".as_bytes().to_vec();
        let mut buf = itoa::Buffer::new();

        eprintln!("{}", std::str::from_utf8(&inbuf).unwrap());
        rlebuf_write(&inbuf, &mut outbuf, &mut buf, is_base).unwrap();
        assert_eq!(std::str::from_utf8(&outbuf).unwrap(), "4A4G");
    }

    #[test]
    fn test_rle_4() {
        let mut outbuf: Vec<u8> = vec![];
        let inbuf = "^[T^[T^[T^A^[T^[T".as_bytes().to_vec();
        let mut buf = itoa::Buffer::new();

        eprintln!("{}", std::str::from_utf8(&inbuf).unwrap());
        rlebuf_write(&inbuf, &mut outbuf, &mut buf, is_base).unwrap();
        assert_eq!(std::str::from_utf8(&outbuf).unwrap(), "3^[T1^A2^[T");
    }

    #[test]
    fn test_rle_5() {
        let mut outbuf: Vec<u8> = vec![];
        let inbuf = "GGGGGG".as_bytes().to_vec();
        let mut buf = itoa::Buffer::new();

        eprintln!("{}", std::str::from_utf8(&inbuf).unwrap());
        rlebuf_write(&inbuf, &mut outbuf, &mut buf, is_base).unwrap();
        assert_eq!(std::str::from_utf8(&outbuf).unwrap(), "6G");
    }

    #[test]
    fn test_rle_qual_1() {
        let mut outbuf: Vec<u8> = vec![];
        let inbuf = "FFFII!!98999".as_bytes().to_vec();
        let mut buf = itoa::Buffer::new();

        eprintln!("{}", std::str::from_utf8(&inbuf).unwrap());
        rlebuf_write(&inbuf, &mut outbuf, &mut buf, is_any).unwrap();
        assert_eq!(std::str::from_utf8(&outbuf).unwrap(), "3F2I2!191839");
    }

    #[test]
    fn test_rle_indel_1() {
        let mut outbuf: Vec<u8> = vec![];
        let inbuf = "AGGG+2GG+2GGAA".as_bytes().to_vec();
        let mut buf = itoa::Buffer::new();

        eprintln!("{}", std::str::from_utf8(&inbuf).unwrap());
        rlebuf_write(&inbuf, &mut outbuf, &mut buf, is_any).unwrap();
        assert_eq!(std::str::from_utf8(&outbuf).unwrap(), "1A3G2[+2GG]2A");
    }
}
