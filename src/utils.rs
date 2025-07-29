use anyhow::{Context, Error};
use rust_htslib::bam::record::Cigar;
use crate::pileup::CigarState;

#[derive(Debug, PartialEq, Eq)]
pub enum CigarAtPos {
    BeforePos(),
    Op(Cigar),
    BaseEmpty(),
}

pub fn has_index(bam_file: &str) -> Result<bool, Error> {
    let potential_index = format! {"{bam_file}.bai"};

    std::fs::exists(&potential_index)
        .with_context(|| format!("Unable to check directory for file {}", &potential_index))
}

/// Get the cigar operation in a read at a given index. Intended to mimic cigar_resolver2 from
/// htslib.
///
/// If the queried index is at the end of a match operation, the function will check if the next
/// operation is a deletion or insertion, and return the corresponding operation if so.
///
/// For example:
///
/// if return == [CigarAtPos(Cigar::Del(l))], then current position is [Cigar::Match] but the very next
/// one is [Cigar::Del].
pub fn cigar_get_pos(cs: &mut CigarState, pos: u32) -> CigarAtPos {
    let cig = &cs.cig;
    let ncig = cig.len();
    let mut op: Cigar;
    while cs.bam_pos <= pos {
        if cs.icig >= ncig {
            // this should never happen, since we check cigars beforehand to at least end
            // at the queried coordinate, if not pass over it.
            return CigarAtPos::BeforePos();
        }

        op = cig[cs.icig];
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let end_pos = cs.bam_pos + len - 1;

                if end_pos < pos {
                    cs.bam_pos += len;
                    cs.iseq += len;
                    cs.icig += 1;
                    continue;
                }

                cs.del = false;
                cs.qpos = pos as usize - cs.bam_pos as usize + cs.iseq as usize;
                if end_pos == pos && cs.icig + 1 < ncig {
                    let next_op = cig[cs.icig + 1];

                    match next_op {
                        Cigar::Ins(_) => return CigarAtPos::Op(next_op),
                        Cigar::Del(_) => return CigarAtPos::Op(next_op),
                        _ => (),
                    }
                }
                return CigarAtPos::Op(Cigar::Match(len));
            }

            Cigar::Ins(len) | Cigar::SoftClip(len) => {
                cs.iseq += len;
                cs.icig += 1;
                continue;
            }

            Cigar::Del(len) => {
                let end_pos = cs.bam_pos + len - 1;
                if end_pos < pos {
                    cs.bam_pos += len;
                    cs.icig += 1;
                    continue;
                }

                // this coordinate comes after we already indicated the deletion, so
                // mark ipos to avoid repeating the deletion in this and subsequent plp cols
                cs.del = true;
                // cs.qpos = (cs.iseq + len) as usize;
                cs.qpos = cs.iseq as usize;
                return CigarAtPos::Op(op);
            }

            Cigar::RefSkip(len) => {
                let end_pos = cs.bam_pos + len - 1;
                if end_pos < pos {
                    cs.bam_pos += len;
                    cs.icig += 1;
                    continue;
                }

                return CigarAtPos::BaseEmpty();
            }
            _ => (),
        }
    }

    CigarAtPos::BaseEmpty()
}

