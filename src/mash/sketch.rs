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

use anyhow::Context;

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
        .par_iter()
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
    inputs: &[PathBuf],
    sketch_size: usize,
    sketch_seed: u64,
    kmer_size: u8,
) -> anyhow::Result<Vec<Sketch>> {
    // Create sketches
    let sketches_path =
        create_sketches(inputs, kmer_size, sketch_size, sketch_seed, "cedar_result")?;

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

/// `cedar sketch`: sketch FASTA files to .msh files in `output_dir`
/// for later reuse with --from-sketches. Mirrors `mash sketch`'s role
/// in Mash's own cli.
pub fn sketch_only(
    inputs: &[PathBuf],
    sketch_size: usize,
    seed: u64,
    kmer_size: u8,
    output_dir: &str,
) -> anyhow::Result<Vec<String>> {
    let paths = create_sketches(inputs, kmer_size, sketch_size, seed, output_dir)?;

    for p in &paths {
        println!("  sketched => {}", p);
    }

    Ok(paths)
}

/// Load pre-computed .msh sketch files from `sketch_dir`, checking that
/// every sketch's own recorder parameters (k-mer size, sketch size, seed)
/// are mutually consistent across the set. Inconsistent parameters would
/// silently produce a mixed-parameter distance matrix.
pub fn load_sketches(sketch_dir: &str) -> anyhow::Result<Vec<Sketch>> {
    let entries: Vec<PathBuf> = fs::read_dir(sketch_dir)?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().map(|e| e == "msh").unwrap_or(false))
        .collect();

    if entries.is_empty() {
        anyhow::bail!(
            "--from-sketches directory '{}' contains no .msh files",
            sketch_dir
        );
    }

    let sketches: Vec<Sketch> = entries
        .into_par_iter()
        .map(|path| {
            let batch = finch::open_sketch_file(&path)
                .with_context(|| format!("Failed to load sketch: {}", path.display()))?;
            Ok(batch)
        })
        .collect::<anyhow::Result<Vec<_>>>()?
        .into_iter()
        .flatten()
        .collect();

    // Consistency check: all sketches in a --from-sketches run
    // must have been created with the same k, sketch size, and seed.
    // Loading a mix would silently produce distance that aren't comparable
    // to each other.
    let first_params = sketches.first().map(|s| s.sketch_params.clone());
    if let Some(ref expected) = first_params {
        let mismatched: Vec<&str> = sketches
            .iter()
            .filter(|s| &s.sketch_params != expected)
            .map(|s| s.name.as_str())
            .collect();

        if !mismatched.is_empty() {
            anyhow::bail!(
                "--from-sketches: the following sketches were created with \
                different parameters (k/sketch-size/seed) than the first \
                sketch in the directory. All sketches in one run must be \
                from the same parameters to produce comparable distances:\n {}",
                mismatched.join("\n  ")
            );
        }
    }

    eprintln!(
        "Loaded {} sketch(es) from '{}' (k = {}, s = {})",
        sketches.len(),
        sketch_dir,
        sketches.first().map(|s| s.sketch_params.k()).unwrap_or(0),
        sketches
            .first()
            .map(|s| sketch_size_of(&s.sketch_params))
            .unwrap_or(0),
    );

    Ok(sketches)
}

/// Extract the final sketch size from SketchParams.
/// Needed for the consistency check log line above.
fn sketch_size_of(params: &SketchParams) -> usize {
    match params {
        SketchParams::Mash { final_size, .. } => *final_size,
        SketchParams::Scaled { .. } => 0,
        SketchParams::AllCounts { .. } => 0,
    }
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
