use anyhow::Error;
use clap::Parser;

mod pileup;
mod read_buf;
mod rpileup;

#[derive(Parser)]
pub struct Args {
    pub input: String,

    #[arg(short = 'A', default_value_t = false)]
    pub show_all: bool,
}

fn main() -> Result<(), Error> {
    let args = Args::parse();

    let mut pileup = rpileup::PileupIterator::new(&args.input, args.show_all, None, None)?;
    let mut ret: rpileup::IterResult;

    loop {
        ret = pileup.init_to_ref()?;

        match ret {
            rpileup::IterResult::NoData => break,
            _ => loop {
                match pileup.next()? {
                    rpileup::IterResult::ReferenceEnd => break,
                    rpileup::IterResult::NoData => panic!(),
                    _ => (),
                }
            },
        }
    }

    Ok(())
}
