#![allow(clippy::unused_io_amount)]

use crate::alignment::PileupAlignment;
use crate::errors::Error;
use crate::output::{OrderedPileupOutput, PileupOutputContext};
use crate::refseq::RefSeqHandle;
use rust_htslib::bam::record::Cigar;
use std::io::Write;

const LAST_POS: u8 = b'$';
const FIRST_POS: u8 = b'^';

const F_MATCH: u8 = b'.';
const R_MATCH: u8 = b',';

const F_REFSKIP: u8 = b'>';
const R_REFSKIP: u8 = b'<';

#[derive(Clone, Debug)]
/// The standard samtools mpileup string. See the samtools mpileup docs for details.
pub struct PileupString {
    pub seqbuf: Vec<u8>,
    pub qualbuf: Vec<u8>,
    pub depth: u32,
}

unsafe impl Send for PileupString {}
unsafe impl Sync for PileupString {}

impl OrderedPileupOutput for PileupString {
    #[inline(always)]
    fn intake(&mut self, ctx: &PileupOutputContext, p: &PileupAlignment) -> Result<(), Error> {
        self.depth += 1;
        self.write_plp(p, ctx.pos, ctx.refseq)?;
        // write_plp(p, ctx.pos, &mut self.seqbuf, &mut self.qualbuf, ctx.refseq)?;
        Ok(())
    }

    #[inline(always)]
    fn write_header<W: std::io::Write>(ctx: &PileupOutputContext, writer: &mut W) -> Result<(), Error> {
        Self::write_refinfo(ctx, writer)
    }

    #[inline(always)]
    fn write_body<W: std::io::Write>(&self, _ctx: &PileupOutputContext, writer: &mut W) -> Result<(), Error> {
        self.write_body(writer)?;

        Ok(())
    }

    fn write_body_empty<W: std::io::Write>(_ctx: &PileupOutputContext, writer: &mut W) -> Result<(), Error> {
        writer.write_all(b"0\t*\t*")?;
        Ok(())
    }

    #[inline(always)]
    fn depth(&self) -> u32 {
        self.depth
    }

    #[inline(always)]
    fn new() -> Self {
        Self::new()
    }

    #[inline(always)]
    fn clear(&mut self) {
        self.depth = 0;
        self.seqbuf.clear();
        self.qualbuf.clear();
    }
}

impl PileupString {
    #[inline(always)]
    pub fn write_refinfo<W: std::io::Write>(ctx: &PileupOutputContext, writer: &mut W) -> Result<(), Error> {
        let mut buf = itoa::Buffer::new();
        writer.write(ctx.ref_name.as_bytes())?;
        writer.write(b"\t")?;

        writer.write(buf.format(ctx.pos + 1).as_bytes())?;
        writer.write(b"\t")?;

        let ref_base = if let Some(seq) = ctx.refseq.as_ref() {
            *seq.get(ctx.pos as usize).unwrap_or(&b'N')
        } else {
            b'N'
        };

        writer.write(&[ref_base])?;
        writer.write_all(b"\t")?;
        Ok(())
    }

    #[inline(always)]
    fn write_body<W: std::io::Write>(&self, writer: &mut W) -> Result<(), Error> {
        let mut buf = itoa::Buffer::new();
        writer.write(buf.format(self.depth).as_bytes())?;
        writer.write(b"\t")?;

        if self.seqbuf.is_empty() {
            writer.write(b"*\t")?;
        } else {
            writer.write(&self.seqbuf)?;
            writer.write(b"\t")?;
        }

        if self.qualbuf.is_empty() {
            writer.write_all(b"*")?;
        } else {
            writer.write_all(&self.qualbuf)?;
        }

        Ok(())
    }

    #[allow(clippy::new_without_default)]
    pub fn new() -> Self {
        Self {
            depth: 0,
            qualbuf: Vec::new(),
            seqbuf: Vec::new(),
        }
    }

    #[inline(always)]
    pub fn write_plp(&mut self, p: &PileupAlignment, pos: i64, refseq: &RefSeqHandle) -> Result<(), Error> {
        write_plp(p, pos, &mut self.seqbuf, &mut self.qualbuf, refseq)
    }
}

// cap qualitites at max of 126; this also helps avoid non-ascii output
#[inline(always)]
pub fn get_qual(qual: u8) -> u8 {
    match qual.cmp(&93).is_gt() {
        true => 126,
        false => qual + 33,
    }
}

// TODO: take arguments that determine verbosity of reported insertion, e.g. full sequence or just
// length?
#[inline(always)]
pub fn expand_insertions(
    p: &PileupAlignment,
    seq_buf: &mut Vec<u8>,
    ndel: &mut i32,
    decorate: bool,
) -> Result<(), Error> {
    let is_rev = p.rec.is_reverse();
    let mut read_pos: usize;
    let mut read_base: u8;

    *ndel = 0;
    // first measure how many insertion bases there are.
    let mut len_indel = 0;
    let ncig = p.cstate.cig.len();
    let mut k = p.cigar_index + 1;
    while k < ncig {
        match p.cstate.cig[k] {
            Cigar::Pad(l) | Cigar::Ins(l) => {
                len_indel += l;
            }
            _ => break,
        }
        k += 1;
    }

    if decorate {
        write!(seq_buf, "+{}", len_indel)?
    };

    // then produce the sequence representing the insertion
    k = p.cigar_index + 1;
    let mut offset = 1;
    while k < ncig {
        match p.cstate.cig[k] {
            Cigar::Pad(l) => seq_buf.extend(std::iter::repeat_n(b'*', l as usize)),
            Cigar::Ins(l) => {
                for _ in 0..l as usize {
                    read_pos = p.qpos + offset - p.del as usize;
                    read_base = if read_pos < p.rec.seq_len() {
                        p.rec.seq()[read_pos]
                    } else {
                        b'n'
                    };
                    offset += 1;
                    match is_rev {
                        true => seq_buf.push(read_base.to_ascii_lowercase()),
                        false => seq_buf.push(read_base.to_ascii_uppercase()),
                    }
                }
            }
            Cigar::Del(l) => {
                *ndel = l as i32;
                break;
            }
            _ => break,
        }
        k += 1;
    }

    Ok(())
}

#[inline(always)]
pub fn write_plp(
    p: &PileupAlignment,
    pos: i64,
    seq_buf: &mut Vec<u8>,
    qualbuf: &mut Vec<u8>,
    refseq: &RefSeqHandle,
) -> Result<(), Error> {
    if p.head {
        seq_buf.push(FIRST_POS);
        seq_buf.push(get_qual(p.rec.mapq()));
    }

    let is_rev = p.rec.is_reverse();
    let mut refbase: u8;

    match p.del {
        false => {
            refbase = if let Some(refseq) = refseq.as_ref() {
                refseq[pos as usize]
            } else {
                b'('
            };
            let readbase = if p.qpos < p.rec.seq_len() {
                p.rec.seq()[p.qpos]
            } else {
                b'n'
            };

            if readbase.eq_ignore_ascii_case(&refbase) {
                match is_rev {
                    true => seq_buf.push(R_MATCH),
                    false => seq_buf.push(F_MATCH),
                }
            } else {
                match is_rev {
                    true => seq_buf.push(readbase.to_ascii_lowercase()),
                    false => seq_buf.push(readbase.to_ascii_uppercase()),
                }
            }
        }

        true => {
            if p.refskip {
                match is_rev {
                    true => seq_buf.push(R_REFSKIP),
                    false => seq_buf.push(F_REFSKIP),
                };
            } else {
                seq_buf.push(b'*');
            }
        }
    }

    let mut del_len = -p.indel;
    if p.indel > 0 {
        expand_insertions(p, seq_buf, &mut del_len, true)?;
    }

    if del_len > 0 {
        write!(seq_buf, "{}", -del_len)?;
        for i in 1..=del_len as i64 {
            refbase = if let Some(refseq) = refseq.as_ref() {
                refseq[(pos + i) as usize]
            } else {
                b'N'
            };

            match is_rev {
                false => seq_buf.push(refbase.to_ascii_uppercase()),
                true => seq_buf.push(refbase.to_ascii_lowercase()),
            }
        }
    }

    if p.tail {
        seq_buf.push(LAST_POS);
    }

    // finally,we add PHRED qual
    qualbuf.push(get_qual(*p.rec.qual().get(p.qpos).unwrap_or(&0)));

    Ok(())
}
