// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use crate::cli;
use crate::dist;
use crate::sketch;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::fs::{self, File};
use std::io::BufRead;
use std::io::BufReader;
use std::io::{self, Write};
use std::path::Path;
use std::process;
use std::sync::Once;

static INIT_RAYON: Once = Once::new();

pub fn init_rayon_pool(threads: usize) {
    INIT_RAYON.call_once(|| {
        ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .expect("Failed to initialize Rayon global thread pool");
    });
}

pub fn output_tree(output: Option<String>, newick: String) -> anyhow::Result<()> {
    match output {
        Some(path) => {
            let mut file = fs::File::create(path)?;
            file.write_all(newick.as_bytes())?;
        }
        None => {
            writeln!(io::stdout(), "{}", newick)?;
        }
    }
    Ok(())
}

pub fn manage_tempdir(
    keep: bool,
    matrix: &speedytree::DistanceMatrix,
    tempdir: &str,
    append: bool,
) -> anyhow::Result<()> {
    if keep {
        let p = Path::new(tempdir);
        dist::to_phylip(matrix, p, append)
    } else {
        fs::remove_dir_all(tempdir)?;
        Ok(())
    }
}

/// Check if provided file is in FASTA format.
///
/// **Input**: `path`: file path
///
/// **Output**: a boolean value. `true` if is FASTA formated, `false` otherwise.
///
pub fn is_fasta_format(path: &str) -> bool {
    let file = File::open(path).expect("file should exists before opening.");
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    while let Some(Ok(line)) = lines.next() {
        let trimmed_line = line.trim();
        if trimmed_line.is_empty() {
            continue;
        }
        return matches!(trimmed_line.chars().next(), Some('>'));
    }
    false
}

/// Check if file is multi FASTA or not.
///
/// **Input**: `path`: file path
///
/// **Output**: a boolean value. `true` if multi fasta, `false` otherwise.
pub fn is_multi_fasta(path: &str) -> bool {
    let file = File::open(path).expect("file should be available");
    let reader = BufReader::new(file);

    let mut header_count = 0;

    for line in reader.lines() {
        let line = line.expect("line should be parseable");
        if line.trim_start().starts_with('>') {
            header_count += 1;
            if header_count > 1 {
                return true;
            }
        }
    }
    false
}

/// Compute sequence length of each provided sequence.
///
/// **Input**: `path`: file path
///
/// **Output**: a `anyhow::Result` of tuple of String (sequence id) and usize (sequence size).
///
pub fn compute_genome_stats(filenames: &[String]) -> anyhow::Result<Vec<(String, usize)>> {
    fn helper(path: &String) -> anyhow::Result<(String, usize)> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        let mut total_len = 0;
        let mut id = String::new();
        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();

            if trimmed.is_empty() {
                continue;
            }

            if let Some(stripped) = trimmed.strip_prefix('>') {
                let tmp_id = stripped
                    .split_whitespace()
                    .next()
                    .ok_or_else(|| anyhow::anyhow!("Malformed fasta header: {}", trimmed))?;
                id = tmp_id.to_string();
            } else {
                total_len += trimmed.len();
            }
        }
        Ok((id, total_len))
    }

    let stats: Vec<(String, usize)> = filenames
        .par_iter()
        .map(helper)
        .collect::<anyhow::Result<_>>()?;

    // print stats
    for stat in &stats {
        println!("Genome: {}, size: {}", stat.0, format_genome_size(stat.1));
    }

    Ok(stats)
}

pub fn format_genome_size(size: usize) -> String {
    let actual = format!("{} bp", size);

    let approx = if size >= 1_000_000_000 {
        let gb = size as f64 / 1_000_000_000.0;
        format!("(~{:.1} Gb)", gb)
    } else if size >= 1_000_000 {
        let mb = size as f64 / 1_000_000.0;
        format!("(~{:.1} Mb)", mb)
    } else if size >= 1_000 {
        let kb = size as f64 >= 1_000.0;
        format!("(~{:.1} Kb)", kb)
    } else {
        "".to_string()
    };
    format!("{} {}", actual, approx)
}

/// Detect influential genome's size.
/// If the set of submitted genomes is less than 4, then use
/// the leave-one-out mean impact method, else use the IQR method.
/// If impactᵢ < ε × μ  → negligible.
///
/// **Inputs:**
///     - `data`: a vec of tupe containing the genome ID and size.
///     - `epsilon`: relative threshold which will be used by multiplying it to the mean (epsilon * meean)
///
/// **Output:** a vec of tuple similar to `data` containing the outliers.
///
pub fn check_genome_outliers(data: &[(String, usize)], epsilon: f64) -> anyhow::Result<()> {
    let mut outliers = Vec::with_capacity(data.len());
    if data.len() < 4 {
        // Use leave-one-out mean impact method
        let values: Vec<f64> = data.iter().map(|x| x.1 as f64).collect();
        let sum: f64 = values.iter().sum();
        let mean = sum / values.len() as f64;
        let mut max_impact = 0.0;
        let mut influential_genome = (String::new(), 0_usize);
        for x in data {
            let leave_one_sum = sum - x.1 as f64;
            let leave_one_mean = leave_one_sum / (values.len() as f64 - 1.0);
            let impact = (mean - leave_one_mean).abs();

            if impact > max_impact {
                max_impact = impact;
                influential_genome = x.clone();
            }
        }

        if max_impact > epsilon * mean {
            outliers.push(influential_genome);
        }
    } else {
        // Use IQR method
        let mut values: Vec<usize> = data.iter().map(|(_, size)| *size).collect();
        values.sort_unstable();
        let len = values.len();

        let q1 = values[len / 4];
        let q3 = values[(3 * len) / 4];
        let iqr = q3 - q1;

        let lower_bound = q1.saturating_sub(3 * iqr / 2);
        let upper_bound = q3 + 3 * iqr / 2;
        outliers = data
            .iter()
            .filter(|x| x.1 < lower_bound || x.1 > upper_bound)
            .cloned()
            .collect();
    }

    if !outliers.is_empty() {
        eprintln!("Error: outliers detected in genome sizes");
        for outlier in outliers {
            eprintln!(
                "Genome {} with size {} will negatively influence k selection",
                outlier.0,
                format_genome_size(outlier.1)
            );
        }
        process::exit(1);
    }

    Ok(())
}

pub fn validate_inputs(filenames: &[String]) -> anyhow::Result<()> {
    if filenames.len() < 3 {
        anyhow::bail!("Input validation error: At least three input FASTA files must be provided.");
    }

    let mut invalid = vec![];
    let mut multi_seq = vec![];

    for file in filenames {
        if !is_fasta_format(file) {
            invalid.push(file.clone());
        } else if is_multi_fasta(file) {
            multi_seq.push(file.clone());
        }
    }

    if !invalid.is_empty() {
        anyhow::bail!(
            "Input validation error: Only FASTA files are allowed. Invalid files: {}",
            invalid.join(", ")
        );
    }

    if !multi_seq.is_empty() {
        anyhow::bail!(
            "Input validation error: Only single-sequence FASTA files are allowed. Multi-sequence files: {}",
            multi_seq.join(", ")
        );
    }

    Ok(())
}

pub fn determine_kmer_size(args: &cli::BuildArgs, stats: &[(String, usize)]) -> u8 {
    if let Some(km) = args.kmer {
        println!("User-defined k-mer size: {}", km);
        km
    } else {
        let mean_genome_size =
            (stats.iter().map(|x| x.1 as u64).sum::<u64>() / stats.len() as u64) as u32;
        let kmer_size = sketch::k_computing(mean_genome_size, 0.01);
        println!(
            "Computed k-mer size (mean genome size: {}, probability: 0.01): {}",
            format_genome_size(mean_genome_size as usize),
            kmer_size
        );
        kmer_size
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_fasta_format_ok() {
        assert!(is_fasta_format("test/bac168.fna"));
    }

    #[test]
    fn test_is_fasta_format_not_ok() {
        assert!(!is_fasta_format("test/test.fq"));
    }
}
