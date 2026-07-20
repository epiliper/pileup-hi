#[cfg(feature = "cli")]
use {crate::engine::MIN_COORDS_PER_THREAD, clap::Parser};

pub const STDOUT_ARG_STR: &str = "STDOUT";
pub const DEFAULT_MPLP_DEPTH: usize = 8000;
pub const DEFAULT_EXCL_FLAGS: u16 = 0x704;
pub const DEFAULT_MIN_BASEQ: u8 = 13;
pub const DEFAULT_MIN_MAPQ: u8 = 0;

#[cfg(feature = "cli")]
fn parse_bam_flags(value: &str) -> Result<u16, String> {
    if value.starts_with(|c: char| c.is_ascii_digit()) {
        let (digits, radix) = if let Some(digits) = value.strip_prefix("0x").or_else(|| value.strip_prefix("0X")) {
            (digits, 16)
        } else if value.len() > 1 && value.starts_with('0') {
            (&value[1..], 8)
        } else {
            (value, 10)
        };

        return u16::from_str_radix(digits, radix).map_err(|_| format!("invalid BAM flag mask: {value}"));
    }

    value.split(',').try_fold(0, |mask, name| {
        let flag = match name.trim().to_ascii_uppercase().as_str() {
            "PAIRED" => 0x1,
            "PROPER_PAIR" | "PROPER-PAIR" => 0x2,
            "UNMAP" | "UNMAPPED" => 0x4,
            "MUNMAP" | "MATE_UNMAPPED" | "MATE-UNMAPPED" => 0x8,
            "REVERSE" => 0x10,
            "MREVERSE" => 0x20,
            "READ1" | "READ-1" => 0x40,
            "READ2" | "READ-2" => 0x80,
            "SECONDARY" => 0x100,
            "QCFAIL" | "QC-FAIL" => 0x200,
            "DUP" | "DUPLICATE" => 0x400,
            "SUPPLEMENTARY" => 0x800,
            _ => return Err(format!("invalid BAM flag: {name}")),
        };
        Ok(mask | flag)
    })
}

#[cfg_attr(feature = "cli", derive(Parser, Clone))]
#[derive(Debug)]
/// Parameters specifying an input file and optional regions for processing.
pub struct InputParams {
    /// File to read (path) or stdout ("-")
    #[cfg_attr(feature = "cli", arg(value_delimiter = ','))]
    pub file: Vec<String>,

    /// only process particular bam regions (comma separated, e.g. chr1:0-8000,chr2)
    #[cfg_attr(feature = "cli", arg(short = 'r', long = "region", value_delimiter = ','))]
    pub regions: Option<Vec<String>>,
}

#[cfg_attr(feature = "cli", derive(Parser))]
#[derive(Clone, Debug)]
/// Parameters for pileup generation.
pub struct PileupParams {
    /// Emit rows for empty coords on a reference with coverage elsewhere
    #[cfg_attr(feature = "cli", arg(short = 'a'))]
    pub show_empty_coords: bool,

    /// Emit rows for all coordinates of all regions
    #[cfg_attr(feature = "cli", arg(long = "aa"))]
    pub show_empty_regions: bool,

    /// Minimum number of coordinates per worker thread. Not recommended to adjust
    #[cfg_attr(feature = "cli", arg(short = 'c', long = "thread-coords", default_value_t = MIN_COORDS_PER_THREAD))]
    pub coords_per_thread: i64,

    /// Reference fasta to use for comparison, must be indexed
    #[cfg_attr(feature = "cli", arg(short = 'f', long = "fasta-ref"))]
    pub refseq: Option<String>,

    /// The maximum number of reads to sample per position. Set to 0 to uncap
    #[cfg_attr(feature = "cli", arg(short = 'd', long = "depth", default_value_t = 8000))]
    pub depth: usize,

    /// Disable R1/R2 mate overlap correction
    #[cfg_attr(feature = "cli", arg(short = 'x'))]
    pub disable_overlaps: bool,

    /// Count reads with unmapped mates
    #[cfg_attr(feature = "cli", arg(short = 'A'))]
    pub count_orphans: bool,

    /// Allow intra-reference threading when depth is fixed. May differ from samtools mpileup output
    #[cfg_attr(feature = "cli", arg(long = "lax-depth", default_value_t = false))]
    pub lax_depth: bool,

    /// Required flags: only include reads with any mask bit set
    #[cfg_attr(
        feature = "cli",
        arg(long = "rf", alias = "incl-flags", value_name = "STR|INT", value_parser = parse_bam_flags, default_value_t = 0)
    )]
    pub incl_flags: u16,

    /// Don't consider any reads with these flags
    #[cfg_attr(
        feature = "cli",
        arg(long = "ff", alias = "excl-flags", value_name = "STR|INT", value_parser = parse_bam_flags, default_value_t = DEFAULT_EXCL_FLAGS)
    )]
    pub excl_flags: u16,

    /// Minimum mapping quality for a read's bases to be counted
    #[cfg_attr(feature = "cli", arg(short = 'q', long = "min-MQ", default_value_t = 0))]
    pub min_mapq: u8,

    /// Minimum phred score for a base to be counted
    #[cfg_attr(feature = "cli", arg(short = 'Q', long = "min-BQ", default_value_t = 13))]
    pub min_baseq: u8,

    /// Disable calcluation of base alignment quality (BAQ)
    #[cfg_attr(feature = "cli", arg(short = 'B', long = "no-BAQ", default_value_t = false))]
    pub no_baq: bool,

    /// Calculate BAQ even when BAQ already exists
    #[cfg_attr(
        feature = "cli",
        arg(short = 'E', long = "redo-BAQ", default_value_t = false, conflicts_with("no_baq"))
    )]
    pub redo_baq: bool,
}

impl InputParams {
    /// Create input parameters from a BAM file/stdin and an optional list of query region. If path is None, then input will be assumed to come from stdin. If regions is None, then all references in the BAM will be processed.
    pub fn new(path: Option<&[&str]>, regions: Option<&[&str]>) -> Self {
        Self {
            file: path
                .map(|r| r.iter().map(|r| r.to_string()).collect::<Vec<String>>())
                .unwrap_or(vec!["STDOUT_ARG_STR".to_string()]),
            regions: regions.map(|r| r.iter().map(|r| r.to_string()).collect::<Vec<String>>()),
        }
    }
}

impl Default for PileupParams {
    fn default() -> Self {
        Self {
            show_empty_coords: false,
            show_empty_regions: false,
            coords_per_thread: MIN_COORDS_PER_THREAD,
            refseq: None,
            depth: DEFAULT_MPLP_DEPTH,
            lax_depth: false,
            disable_overlaps: false,
            count_orphans: false,
            incl_flags: 0,
            excl_flags: DEFAULT_EXCL_FLAGS,
            min_mapq: DEFAULT_MIN_MAPQ,
            min_baseq: DEFAULT_MIN_BASEQ,
            no_baq: false,
            redo_baq: false,
        }
    }
}

#[cfg(all(test, feature = "cli"))]
mod tests {
    use super::parse_bam_flags;

    #[test]
    fn parses_numeric_bam_flag_masks() {
        assert_eq!(parse_bam_flags("16").unwrap(), 16);
        assert_eq!(parse_bam_flags("0x10").unwrap(), 16);
        assert_eq!(parse_bam_flags("020").unwrap(), 16);
        assert_eq!(parse_bam_flags("0").unwrap(), 0);
        assert_eq!(parse_bam_flags("65535").unwrap(), u16::MAX);
    }

    #[test]
    fn parses_named_bam_flag_masks() {
        assert_eq!(parse_bam_flags("READ2,REVERSE").unwrap(), 0x90);
        assert_eq!(parse_bam_flags("proper_pair,mreverse").unwrap(), 0x22);
        assert_eq!(parse_bam_flags("read-1,qc-fail,duplicate").unwrap(), 0x640);
    }

    #[test]
    fn rejects_invalid_bam_flag_masks() {
        assert!(parse_bam_flags("0xg").is_err());
        assert!(parse_bam_flags("65536").is_err());
        assert!(parse_bam_flags("REVERSE,UNKNOWN").is_err());
    }
}
