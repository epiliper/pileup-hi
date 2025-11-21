#![allow(dead_code)]
use crate::{alignment::PileupAlignment, threading::PileupWorkerState};
use anyhow::Error;
use crossbeam::channel::{bounded, Receiver, Sender};

const PILEUP_OUTPUT_BUF_PURGE_THRES: usize = 10_000;

/// The interface requirements for a pileup output. It needs to give ref information,
/// intake pileup alignments, update current ref info, display depth, and write itself.
pub trait OrderedPileupOutput: Send + Sync + Clone {
    /// Get the reference of the pileup
    fn tid(&self) -> i32;
    /// Get the coordinate of the pileup
    fn pos(&self) -> i64;
    /// Update internal data with pileup alignment
    fn intake(&mut self, p: &PileupAlignment, refseq: Option<&[u8]>) -> Result<(), Error>;
    /// Update reference data given ref num, pos, name, and sequence
    fn set_ref_info(&mut self, tid: i32, pos: i64, ref_name: &str, ref_seq: Option<&[u8]>);
    fn write<W: std::io::Write>(&mut self, writer: &mut W) -> Result<(), Error>;
    fn depth(&self) -> u32;
}

/// Defines how to get output data from iterators from a thread. If using a single thread, we dont'
/// have to care about queue-ing output.
pub enum OutputMethod<W: std::io::Write, T: OrderedPileupOutput> {
    WriteDirectly(W),
    QueueForOutput(Sender<T>),
}

////////////////
// Begin defs for PileupOutputAggregator
////////////////

pub enum PileupOutputState<T: OrderedPileupOutput> {
    Closed,
    Open(Sender<T>),
}

pub struct PileupOutputAggregator<T>
where
    T: OrderedPileupOutput,
{
    pub input_state: PileupOutputState<T>,
    pub worker_state: PileupWorkerState,
}

impl<T: OrderedPileupOutput + 'static> PileupOutputAggregator<T> {
    pub fn new() -> Self {
        Self {
            input_state: PileupOutputState::Closed,
            worker_state: PileupWorkerState::Off,
        }
    }

    pub fn get_output_handle(&self) -> Option<Sender<T>> {
        match &self.input_state {
            PileupOutputState::Closed => None,
            PileupOutputState::Open(s) => Some(s.clone()),
        }
    }

    pub fn terminate(self) -> Result<(), Error> {
        match (self.input_state, self.worker_state) {
            (PileupOutputState::Closed, _) | (_, PileupWorkerState::Off) => {
                anyhow::bail!("Cannot terminate an output channel that never started!")
            }

            (PileupOutputState::Open(s), PileupWorkerState::Running(j)) => {
                drop(s);
                j.join().unwrap();
                Ok(())
            }
        }
    }

    pub fn run(&mut self) {
        let (s, r): (Sender<T>, Receiver<T>) = bounded(10_000_000);
        let j = std::thread::spawn(move || {
            let mut writer = std::io::stdout().lock();
            let mut output_queue: Vec<T> = Vec::with_capacity(PILEUP_OUTPUT_BUF_PURGE_THRES);

            while let Ok(out) = r.recv() {
                output_queue.push(out);

                if output_queue.len() >= PILEUP_OUTPUT_BUF_PURGE_THRES {
                    // output_queue.sort_by(|a, b| a.tid().cmp(&b.tid()));

                    for mut out in output_queue.drain(..) {
                        out.write(&mut writer).unwrap();
                    }

                    // output_queue.shrink_to(0);
                }
            }

            // output_queue.sort_by(|a, b| a.tid().cmp(&b.tid()));
            output_queue
                .drain(..)
                .for_each(|mut item| item.write(&mut writer).unwrap())
        });

        self.worker_state = PileupWorkerState::Running(j);
        self.input_state = PileupOutputState::Open(s.clone());
    }
}
