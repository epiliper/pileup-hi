#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

mod args;

use pileuphi_lib::{
    error::{Error, ErrorKind},
    outputs::BaseDepthString,
    outputs::PileupString,
    setup_exit_handler, PileupEngine,
};

use crate::args::{parse_or_quit, Commands};
use log::error;

#[cfg(debug_assertions)]
use log::warn;

fn _main() -> Result<(), Error> {
    let params = parse_or_quit();
    setup_exit_handler();

    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    #[cfg(debug_assertions)]
    warn!("Running the debug build of pileup-hi! This is slow.");

    match params.command {
        Commands::Plp(params) => {
            let engine = PileupEngine::<PileupString>::new(params.plp)?;
            engine.write(params.inp, &params.output, params.threads)?;
        }

        Commands::Histo(params) => {
            let engine = PileupEngine::<BaseDepthString>::new(params.plp)?;
            engine.write(params.inp, &params.output, params.threads)?;
        }
    };

    Ok(())
}

fn main() {
    if let Err(e) = _main() {
        if let ErrorKind::IOError(e) = e.kind() {
            if matches!(e.kind(), std::io::ErrorKind::BrokenPipe) {
                std::process::exit(0);
            }
        }

        error!("Error: {e}");
        std::process::exit(1);
    }
}
