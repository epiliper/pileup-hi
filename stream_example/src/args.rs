use clap::Parser;
use pileuphi_lib::param::{InputParams, PileupParams, STDOUT_ARG_STR};

#[derive(Parser, Clone)]
pub struct Params {
    #[clap(flatten)]
    pub inp: InputParams,

    #[clap(flatten)]
    pub plp: PileupParams,

    /// Output to write to
    #[arg(short = 'o', long = "output", default_value_t = STDOUT_ARG_STR.to_string())]
    pub output: String,
}

pub fn parse_or_quit() -> Params {
    match Params::try_parse() {
        Ok(p) => {
            // no argument checking at the moment, leaving here for the future.
            p
        }
        Err(e) => {
            e.print().unwrap();
            std::process::exit(1)
        }
    }
}
