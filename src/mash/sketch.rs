// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::fs;
use std::fs::File;
use std::path::PathBuf;

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
    fs::create_dir_all(outdir)?;

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
            let filename_str: String = filename.display().to_string();
            let sketches = sketch_files(&[&filename_str], &sketch_params, &filter_params)?;

            let basename = filename
                .file_name()
                .map(|s| s.to_string_lossy().into_owned())
                .unwrap_or_else(|| filename_str.clone());

            let out_path = PathBuf::from(outdir).join(format!("{}.msh", basename));
            let mut out_file = File::create(&out_path)?;
            write_mash_file(&mut out_file, &sketches)?;
            Ok(out_path.to_string_lossy().into_owned())
        })
        .collect()
}

pub fn create_and_load_sketches(
    indir: &str,
    sketch_size: usize,
    sketch_seed: u64,
    kmer_size: u8,
) -> anyhow::Result<Vec<Sketch>> {
    let mut inputs = Vec::new();

    for file in fs::read_dir(indir)? {
        inputs.push(file?.path());
    }

    // Create sketches
    let sketches_path =
        create_sketches(&inputs, kmer_size, sketch_size, sketch_seed, "cedar_result")?;

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
mod tests {
    use super::*;
    use std::fs;
    use std::path::Path;

    #[test]
    fn test_create_sketches() {
        let filenames = [
            PathBuf::from("test/bacam.fna"),
            PathBuf::from("test/bacsp.fna"),
        ];
        let kmer_size = 21;
        let sketch_size = 1000;
        let seed = 42;
        let outdir = "test_output";
        fs::create_dir(outdir).unwrap();

        let result = create_sketches(&filenames, kmer_size, sketch_size, seed, outdir);
        assert!(result.is_ok());

        assert!(fs::metadata(outdir).is_ok());
        for filename in &filenames {
            let output_filename = format!(
                "{}/{}.msh",
                outdir,
                Path::new(filename).file_name().unwrap().to_str().unwrap()
            );
            println!("{output_filename}");
            assert!(fs::metadata(output_filename).is_ok());
        }
        fs::remove_dir_all(outdir).unwrap();
    }

    #[test]
    fn test_create_sketches_creates_missing_outdir() {
        // The bug this guards against: create_sketches previously assumed
        // outdir already existed, silently masked by this same test
        // manually pre-creating it in test_create_sketches above. This one
        // deliberately does NOT pre-create the directory.
        let filenames = [PathBuf::from("test/bacam.fna")];
        let outdir = "test_output_no_precreate";
        assert!(
            fs::metadata(outdir).is_err(),
            "outdir must not exist before the call"
        );

        let result = create_sketches(&filenames, 21, 1000, 42, outdir);
        assert!(
            result.is_ok(),
            "create_sketches should create its own outdir: {:?}",
            result.err()
        );
        assert!(fs::metadata(outdir).is_ok());

        fs::remove_dir_all(outdir).unwrap();
    }

    #[test]
    fn test_create_sketches_output_is_basename_only_not_nested() {
        // The bug itself: an input path with a directory prefix must not
        // produce a nested output path under outdir.
        let filenames = [PathBuf::from("test/bacam.fna")];
        let outdir = "test_output_basename_check";
        fs::create_dir(outdir).unwrap();

        let result = create_sketches(&filenames, 21, 1000, 42, outdir).unwrap();

        assert_eq!(result.len(), 1);
        assert_eq!(result[0], format!("{}/bacam.fna.msh", outdir));
        assert!(
            !result[0].contains("test/bacam.fna.msh"),
            "output path incorrectly nested the input directory: {}",
            result[0]
        );

        fs::remove_dir_all(outdir).unwrap();
    }
}
