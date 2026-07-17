use crate::{
    bamio::{BamDataSource, BamReader, OutputDataDest},
    errors::{Error, ErrorKind},
    jobqueue::IntervalJobs,
    output::{write_multiple_outputs, OrderedPileupOutput},
    params::{InputParams, PileupParams},
    pileup_iterator::PileupIterator,
    position_queue::{create_region_queue, intervals_from_header, GenomeInterval},
    refseq::{RefSeq, RefSeqHandle},
    threading::ThreadPool,
    utils::OutputWriter,
};

use log::{info, warn};
use std::marker::PhantomData;
use std::sync::Arc;

pub const BUFWRITER_CAP: usize = 2 * 1024 * 1024;
pub const MIN_BAM_READ_THREADS: usize = 2;

/// The default minimum number of coordinates to give each thread for processing.
/// This basically exists to prevent doing unnecessary work for very small regions.
/// Can be overridden if you need more horsepower for, say, high-depth regions.
pub const MIN_COORDS_PER_THREAD: i64 = 250_000;

struct PileupEngineQuery {
    intervals: Vec<GenomeInterval>,
    src: Vec<BamDataSource>,
}

impl TryFrom<InputParams> for PileupEngineQuery {
    type Error = Error;
    fn try_from(value: InputParams) -> Result<Self, Error> {
        let mut src = Vec::with_capacity(value.file.len());
        for file in value.file {
            src.push(BamDataSource::from_string(&file)?);
        }

        let tempreader = BamReader::new(&src[0], 1)?;
        let header = &tempreader.header;

        let intervals = if let Some(region) = value.regions {
            create_region_queue(&region, header)?
        } else {
            intervals_from_header(header)?
        };

        Ok(Self { intervals, src })
    }
}

pub struct PileupEngine<T: OrderedPileupOutput> {
    plp_params: PileupParams,
    _t: PhantomData<T>,
    refseq: Option<RefSeq>,
}

///////////////////////////////////////////////////////////////

impl<T: OrderedPileupOutput + 'static> PileupEngine<T> {
    /// Create a reusable pileup engine with the given generation parameters.
    pub fn new(plp_params: PileupParams) -> Result<Self, Error> {
        let refseq = if let Some(ref fasta) = plp_params.refseq {
            if !std::fs::exists(std::path::Path::new(fasta))? {
                return Err(Error::from(ErrorKind::IOError(std::io::Error::new(
                    std::io::ErrorKind::NotFound,
                    "Fasta file {fasta} doesn't exist!",
                ))));
            }
            Some(RefSeq::new(fasta.clone()))
        } else {
            None
        };

        Ok(Self {
            plp_params,
            _t: PhantomData,
            refseq,
        })
    }

    /// Return one pileup iterator per requested region.
    pub fn iter(&self, input: InputParams) -> Result<Vec<PileupIterator<T>>, Error> {
        let query: PileupEngineQuery = input.try_into()?;
        let mut ret = Vec::with_capacity(query.intervals.len());

        for interval in query.intervals.iter() {
            ret.push(PileupIterator::from_query(
                &query.src,
                self.get_refseq(&interval.name)?,
                interval,
                &self.plp_params,
            )?);
        }

        Ok(ret)
    }

    /// Generate pileups and write them to a file or stdout.
    pub fn write(&self, input: InputParams, output: &str, threads: usize) -> Result<(), Error> {
        if threads == 0 {
            return Err(Error::from(ErrorKind::IOError(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "thread count must be greater than zero",
            ))));
        }

        let query: PileupEngineQuery = input.try_into()?;
        let dest = OutputDataDest::from_string(output);

        if let OutputDataDest::File(ref f) = dest {
            if std::fs::exists(f)? {
                warn!("Output file {} already exists! Overwriting...", f);

                if let Err(e) = std::fs::remove_file(f) {
                    warn!("Failed to remove file {f}; {e}. Output will be appended...");
                };
            }
        }

        let index_found = query.src.iter().all(|f| f.has_index().unwrap_or(false));

        if index_found {
            info!("Found index for for input files");
        }

        if threads == 1 {
            self.run_all_1t(&query, &dest)
        } else if !index_found {
            warn!(
                "User asked for more than {} threads but at least one input file is unindexed. Running in single-threaded mode",
                threads
            );
            self.run_all_1t(&query, &dest)
        } else {
            info!("Running with {} threads...", threads);
            self.run_all_par(&query, &dest, threads)
        }
    }

    fn get_refseq(&self, ref_name: &str) -> Result<RefSeqHandle, Error> {
        if let Some(ref refseq) = self.refseq {
            refseq.yield_handle(ref_name)
        } else {
            Ok(Arc::new(None))
        }
    }

    /// Use a single thread for both processing and writing.
    fn run_all_1t(&self, query: &PileupEngineQuery, dest: &OutputDataDest) -> Result<(), Error> {
        for interval in query.intervals.iter() {
            let mut main_writer = OutputWriter::new(dest, BUFWRITER_CAP, true, false)?;

            let mut iterator = PileupIterator::<T>::from_query(
                &query.src,
                self.get_refseq(&interval.name)?,
                interval,
                &self.plp_params,
            )?;

            while iterator.advance()?.is_some() {
                write_multiple_outputs(&iterator.ctx(), iterator.current(), main_writer.get())?;
            }

            main_writer.flush()?;
        }
        Ok(())
    }

    /// Split up a list of input genomic intervals into smaller chunks to be processed in parallel. Chunks are first written to temporary output files before being merged into the user-specified output file.
    fn run_all_par(&self, query: &PileupEngineQuery, dest: &OutputDataDest, threads: usize) -> Result<(), Error> {
        // if we choose to limit read buffers to a given depth, then we need to process at least
        // one entire reference per thread: we can't split up references into sub-chunks per
        // thread.
        //
        // This is because the reads entering buffers at
        // coordinate x are contingent on the reads that entered in coordinate x - 1, x - 2, and
        // so on. If we wanted to get the same output with depth-limited buffers with multiple
        // threads, each thread would have to go through the entire ENTIRE reference, which
        // defeats the point of multithreading.
        //
        // I'm going to go with this constraint for now since we need to remain identical to
        // samtools mpileup output.

        let sub_ref_split_strat = if self.plp_params.depth != 0 && !self.plp_params.lax_depth {
            warn!("Depth was bound to {}, so one thread will be used per reference. For generally faster sub-reference threading, either set depth to 0 (disable) or enable lax_depth.", self.plp_params.depth);
            None
        } else {
            Some(self.plp_params.coords_per_thread)
        };

        // let sub_ref_split_strat = Some(self.plp_params.coords_per_thread);

        let mut jobs = IntervalJobs::new(&query.intervals, sub_ref_split_strat, threads as i64, dest.clone());

        // if we just get one job after splitting, just use one thread.
        if jobs.queue.len() == 1 {
            return self.run_all_1t(query, dest);
        }

        let mut pool = ThreadPool::new(threads);
        let mut n_jobs = 0;

        while !jobs.is_completed() {
            jobs.merge_completed()?;

            if let Some(worker) = pool.get_available() {
                if let Some(job) = jobs.queue.pop_front() {
                    n_jobs += 1;

                    let refseq_handle = self.get_refseq(&job.interval.name)?;

                    worker.run::<T>(n_jobs, self.plp_params.clone(), job, query.src.clone(), refseq_handle);
                }
            }
        }
        jobs.conclude()?;
        Ok(())
    }
}
