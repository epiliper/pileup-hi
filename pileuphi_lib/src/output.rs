use crate::alignment::PileupAlignment;
use crate::errors::Error;
use crate::refseq::RefSeqHandle;

#[allow(type_alias_bounds)]
#[derive(Debug)]
pub enum PileupCoordinate<'a, T: OrderedPileupOutput> {
    NoCoverage,
    Coverage(&'a T),
}

pub enum PileupCoordinateType {
    NoCoverage,
    Coverage,
}

pub struct PileupOutputContext<'a> {
    pub tid: i32,
    pub pos: i64,
    pub ref_name: &'a str,
    pub refseq: &'a RefSeqHandle,
}

/// The interface requirements for a pileup output. It needs to give ref information,
/// intake pileup alignments, update current ref info, display depth, and write itself.
pub trait OrderedPileupOutput: Send + Sync + Clone + std::fmt::Debug {
    /// Update internal data with pileup alignment
    fn intake(&mut self, ctx: &PileupOutputContext, p: &PileupAlignment) -> Result<(), Error>;

    fn write_header<W: std::io::Write>(ctx: &PileupOutputContext, writer: &mut W) -> Result<(), Error>;

    fn write_body<W: std::io::Write>(&self, ctx: &PileupOutputContext, writer: &mut W) -> Result<(), Error>;

    fn write_body_empty<W: std::io::Write>(ctx: &PileupOutputContext, writer: &mut W) -> Result<(), Error>;

    fn depth(&self) -> u32;

    fn clear(&mut self);

    #[allow(dead_code)]
    fn new() -> Self;
}

pub fn write_multiple_outputs<'a, T: OrderedPileupOutput + 'a, W: std::io::Write>(
    ctx: &PileupOutputContext,
    iter: impl Iterator<Item = Option<&'a T>>,
    writer: &mut W,
) -> Result<(), Error> {
    T::write_header(ctx, writer)?;
    for data in iter {
        if let Some(inner) = data {
            inner.write_body(ctx, writer)?;
        } else {
            T::write_body_empty(ctx, writer)?;
        }

        writer.write_all(b"\t")?;
    }

    writeln!(writer)?;

    Ok(())
}
