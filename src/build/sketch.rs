use std::fs;
use std::fs::File;
use std::path::PathBuf;

use crate::cli;
use finch::serialization::Sketch;
use finch::{
    errors::FinchResult, filtering::FilterParams, serialization::write_mash_file, sketch_files,
    sketch_schemes::SketchParams,
};
use rayon::iter::IntoParallelIterator;
use rayon::prelude::*;

/// Compute the value of k that minimizes the probability of
/// observing a random k-mer.
///
/// s: genome size
/// p: desired probabilty
/// Based on Fofanov et al., 2004, 10.1093/bioinformatics/bth266
pub fn k_computing(s: u32, p: f64) -> u8 {
    let x: f64 = s as f64 * (1.0f64 - p) / p;
    (x.log10() / 4.0f64.log10()).ceil() as u8
}

/// Create sketches from fasta files
pub fn create_sketches(
    filenames: &[PathBuf],
    kmer_size: u8,
    sketch_size: usize,
    seed: u64,
    outdir: &str,
) -> FinchResult<Vec<String>> {
    // Create SketchParams struct for finch
    let sketch_params = SketchParams::Mash {
        kmers_to_sketch: sketch_size * 200,
        final_size: sketch_size,
        no_strict: false,
        kmer_length: kmer_size,
        hash_seed: seed,
    };

    // Create FilterParams struct for finch
    let filter_params = FilterParams {
        filter_on: Some(false),
        abun_filter: (Some(0u32), None),
        err_filter: 1.0,
        strand_filter: 0.1,
    };

    // Process files and generate sketches
    filenames
        .iter()
        .map(|filename| {
            let filename: String = filename.display().to_string();
            let sketches = sketch_files(&[&filename], &sketch_params, &filter_params)?;
            let out_path = PathBuf::from(outdir).join(format!("{}.msh", filename));
            let mut out_file = File::create(&out_path)?;
            write_mash_file(&mut out_file, &sketches)?;
            Ok(out_path.to_string_lossy().into_owned())
        })
        .collect()
}

pub fn create_and_load_sketches(
    args: &cli::BuildArgs,
    kmer_size: u8,
) -> anyhow::Result<Vec<Sketch>> {
    let mut inputs = Vec::new();

    for file in fs::read_dir(&args.indir)? {
        inputs.push(file?.path());
    }

    // Create sketches
    let sketches_path = create_sketches(&inputs, kmer_size, args.size, args.seed, "cedar_result")?;

    // Load sketches in parallel
    let sketches: Vec<Sketch> = sketches_path
        .into_par_iter()
        .try_fold(Vec::new, |mut acc, path| -> anyhow::Result<Vec<Sketch>> {
            let mut batch = finch::open_sketch_file(path)?;
            acc.append(&mut batch);
            Ok(acc)
        })
        .try_reduce(Vec::new, |mut acc, mut batch| {
            acc.append(&mut batch);
            Ok(acc)
        })?;
    Ok(sketches)
}

#[cfg(test)]
mod tests {}
