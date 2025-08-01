// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::fs::File;
use std::path::{Path, PathBuf};

use finch::{
    errors::FinchResult, filtering::FilterParams, serialization::write_mash_file, sketch_files,
    sketch_schemes::SketchParams,
};

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
    filenames: &[&str],
    kmer_size: u8,
    sketch_size: usize,
    oversketch: usize,
    seed: u64,
    outdir: &str,
) -> FinchResult<Vec<String>> {
    // Create SketchParams struct for finch
    let sketch_params = SketchParams::Mash {
        kmers_to_sketch: sketch_size * oversketch,
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
            let sketches = sketch_files(&[*filename], &sketch_params, &filter_params)?;
            let out_path = PathBuf::from(outdir).join(format!(
                "{}.msh",
                Path::new(filename).file_name().unwrap().to_string_lossy()
            ));
            let mut out_file = File::create(&out_path)?;
            write_mash_file(&mut out_file, &sketches)?;
            Ok(out_path.to_string_lossy().into_owned())
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_create_sketches() {
        // Define test parameters
        let filenames = vec!["test/bacam.fna", "test/bacsp.fna"];
        let kmer_size = 21;
        let sketch_size = 1000;
        let oversketch = 200;
        let seed = 42;
        let outdir = "test_output";
        fs::create_dir(outdir).unwrap();

        // Call the function under test
        let result = create_sketches(&filenames, kmer_size, sketch_size, oversketch, seed, outdir);
        // Verify that the function returned successfully
        assert!(result.is_ok());

        // Verify that the output directory and sketch files were created
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
}
