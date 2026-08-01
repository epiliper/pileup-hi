#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

mod args;

use pileuphi_lib::{PileupEngine, outputs::RLEPileupString, write_multiple_outputs};

use crate::args::parse_or_quit;

fn main() {
    let params = parse_or_quit();

    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let engine = PileupEngine::<RLEPileupString>::new(params.plp).expect("failed to init pileup engine");

    let iters = engine.iter(params.inp).expect("Failed to get iterator");

    for mut iter in iters {
        while iter.advance().expect("error during pileup generation").is_some() {
            write_multiple_outputs(&iter.ctx(), iter.current(), &mut std::io::stdout()).expect("Error writing");
        }
    }
}
