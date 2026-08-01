//! # pileup-hi_lib: parallel, modular API for pileup generation from SAM/BAM/CRAM
//!
//! pileup-hi_lib is an implementation of a pileup engine that can write directly to a user-defined
//! output format, and supports memory-efficient parallelization the intra-reference level.
//!
//! This engine is designed to emit identical output to samtools mpileup when used with the
//! [`outputs::PileupString`] output type.
//!
//! ## engine API examples
//! ```rust
//! use pileuphi_lib::{
//!         PileupEngine,
//!         param::{InputParams, PileupParams},
//!         outputs::PileupString,
//!         error::Error,
//!         write_multiple_outputs,
//!     };
//!
//! fn run() -> Result<(), Error> {
//!     // process a region chr1 in test.bam
//!     let inputs = InputParams::new(
//!         Some(vec!["test.bam"].as_slice()),
//!         Some(vec!["chr1:1-1000000"].as_slice())
//!         );
//!
//!     let mut params = PileupParams::default();
//!     params.depth = 0; // remove depth limit, process all bases per coordinate
//!     let engine = PileupEngine::<PileupString>::new(params)?;
//!
//!     { // example 1: streaming iterator (non-allocating)
//!         let iters = engine.iter(inputs.clone())?;
//!
//!         // one iterator per reference (just one in this example, chr1)
//!         for mut iter in iters {
//!             while iter.advance()?.is_some() {
//!                 // if you passed in multiple bams, use [`write_multiple_outputs`] for
//!                 // stacking output columns horizontally
//!                 write_multiple_outputs(&iter.ctx(), iter.current(), &mut std::io::stdout())?;
//!
//!                 // you can also iterate over each input bam's pileup like so:
//!                 for sample in iter.current() {
//!                     sample.map(|s| eprintln!("{:?}", s.depth));
//!                 }
//!             }
//!         }
//!     }
//!
//!     { // example 2: parallel write directly to file
//!         engine.write(inputs.clone(), "output.txt", /* threads */ 8)?;
//!     }
//!
//!     Ok(())
//! }
//! ```
//! ## custom output formats
//! You can implement your own output formats to use with [`PileupEngine`] via the
//! [`output::OrderedPileupOutput`] trait. See [`outputs::BaseDepthString`] and
//! [`outputs::PileupString`] for examples of implementations.

pub(crate) mod alignment;
pub(crate) mod bamio;
pub(crate) mod baq;
pub(crate) mod basedepth_string;
pub(crate) mod cigar_resolve;
pub(crate) mod engine;
pub(crate) mod errors;
pub(crate) mod jobqueue;
pub(crate) mod output;
pub(crate) mod overlap;
pub(crate) mod params;
pub(crate) mod pileup_iterator;
pub(crate) mod pileup_string;
pub(crate) mod position_queue;
pub(crate) mod read_buf;
pub(crate) mod read_filter;
pub(crate) mod refseq;
pub(crate) mod rle_pileup_string;
pub(crate) mod threading;
pub(crate) mod utils;

pub use engine::PileupEngine;
pub use jobqueue::setup_exit_handler;
pub use output::write_multiple_outputs;

pub mod error {
    pub use crate::errors::{Error, ErrorKind};
}

pub mod outputs {
    pub use crate::basedepth_string::BaseDepthString;
    pub use crate::output::OrderedPileupOutput;
    pub use crate::pileup_string::PileupString;
    pub use crate::rle_pileup_string::RLEPileupString;
}

pub mod param {
    pub use crate::params::{InputParams, PileupParams, STDOUT_ARG_STR};
}
