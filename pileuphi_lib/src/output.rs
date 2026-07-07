use crate::alignment::PileupAlignment;
use crate::errors::Error;
use crate::refseq::RefSeqHandle;

#[allow(type_alias_bounds)]
pub enum PileupCoordinate<'a, T: OrderedPileupOutput> {
    NoCoverage,
    Coverage(&'a T),
}

/// The interface requirements for a pileup output. It needs to give ref information,
/// intake pileup alignments, update current ref info, display depth, and write itself.
pub trait OrderedPileupOutput: Send + Sync + Clone + std::fmt::Debug {
    /// Get the reference of the pileup
    #[allow(dead_code)]
    fn tid(&self) -> i32;

    /// Get the coordinate of the pileup
    #[allow(dead_code)]
    fn pos(&self) -> i64;

    /// Update internal data with pileup alignment
    fn intake(&mut self, p: &PileupAlignment, refseq: &RefSeqHandle) -> Result<(), Error>;

    /// Update reference data given ref num, pos, name, and sequence
    fn set_ref_info(&mut self, tid: i32, pos: i64, ref_name: &str, refseq: &RefSeqHandle);

    fn write<W: std::io::Write>(&self, writer: &mut W) -> Result<(), Error>;

    fn depth(&self) -> u32;

    fn clear(&mut self);

    #[allow(dead_code)]
    fn new() -> Self;
}
