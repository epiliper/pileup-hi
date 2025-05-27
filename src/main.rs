use anyhow::Error;
use clap::Parser;

mod read_buf;
mod rpileup;

#[derive(Parser)]
pub struct Args {
    pub input: String,
}

fn main() -> Result<(), Error> {
    let args = Args::parse();

    let mut pileup = rpileup::PileupIterator::new(&args.input, None, None)?;
    let mut ret: rpileup::IterResult;
    loop {
        ret = pileup.next()?;
        match ret {
            rpileup::IterResult::NoData | rpileup::IterResult::ReferenceEnd => break,
            _ => (),
        }
    }

    Ok(())
}
