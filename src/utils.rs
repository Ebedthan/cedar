use crate::cli;
use crate::dist::rescue::RescueSummary;
use crate::mash::sketch::k_computing;
use crate::mash::uncertainty::DistanceEstimate;
use crate::mash::uncertainty::Reliability;

use anyhow::{Context, Result};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

use std::fs::{self, File};
use std::io::BufRead;
use std::io::BufReader;
use std::io::{self, Write};
use std::path::PathBuf;
use std::sync::Once;

const RULE_WIDTH: usize = 64;
const LABEL_WIDTH: usize = 24;

/// Default relative-deviation threshold for flagging a genome-size
/// outlier, shared between `build` and `dist` (which otherwise has no
/// reason to depend on `build::PhyloConfig` for a single constant).
pub const DEFAULT_OUTLIER_THRESHOLD: f64 = 0.01;

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

fn fmt_opt(value: Option<f64>) -> String {
    match value {
        Some(v) => format!("{:.6}", v),
        None => "NA".to_string(),
    }
}

/// Write a `cedar dist` distance-with-uncertainty table as TSV, to a file
/// or to stdout if no output path was given.
pub fn output_distance_table(
    output: Option<String>,
    estimates: &[DistanceEstimate],
    confidence: f64,
    rescue_pvalue: f64,
) -> anyhow::Result<()> {
    let mut writer: Box<dyn Write> = match output {
        Some(path) => Box::new(fs::File::create(path)?),
        None => Box::new(io::stdout()),
    };

    // Leading metadata line: makes the file self-describing once separated
    // from the command that produced it, without repeating these two
    // constant values on every data row. Standard `#`-prefixed convention.
    writeln!(
        writer,
        "# confidence={} rescue_pvalue={}",
        confidence, rescue_pvalue
    )?;

    writeln!(
        writer,
        "genome1\tgenome2\tjaccard\tjaccard_ci_low\tjaccard_ci_high\tmash_distance\tmash_distance_ci_low\tmash_distance_ci_high\tshared_hashes\ttotal_hashes\trelative_uncertainty\tflag\tkmer_size_used\trescued"
    )?;

    for e in estimates {
        writeln!(
            writer,
            "{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}",
            e.query,
            e.reference,
            e.jaccard,
            e.jaccard_ci_low,
            e.jaccard_ci_high,
            e.mash_distance,
            e.mash_distance_ci_low,
            e.mash_distance_ci_high,
            e.shared_hashes,
            e.total_hashes,
            fmt_opt(e.relative_uncertainty),
            e.reliability,
            e.kmer_size_used,
            e.rescued,
        )?;
    }

    Ok(())
}

/// Check if provided file is in FASTA format.
///
/// **Input**: `path`: file path
///
/// **Output**: a boolean value. `true` if is FASTA formated, `false` otherwise.
///
pub fn is_fasta_format(path: &PathBuf) -> bool {
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
pub fn is_multi_fasta(path: &PathBuf) -> bool {
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
pub fn compute_genome_stats(filenames: &[PathBuf]) -> anyhow::Result<Vec<(String, usize)>> {
    fn helper(path: &PathBuf) -> anyhow::Result<(String, usize)> {
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
        let kb = size as f64 / 1_000.0;
        format!("(~{:.1} Kb)", kb)
    } else {
        "".to_string()
    };
    format!("{} {}", actual, approx)
}

/// Detect influential genome's size.
///
/// If the set of submitted genomes is less than 4, then use
/// the leave-one-out mean impact method, else use the IQR method.
/// If impactᵢ < ε × μ  => negligible.
///
/// **Inputs:**
///     - `data`: a vec of tupe containing the genome ID and size.
///     - `epsilon`: relative threshold which will be used by multiplying it to the mean (epsilon * meean)
///
/// **Output:** a vec of tuple similar to `data` containing the outliers.
///
pub fn check_genome_outliers(data: &[(String, usize)], epsilon: f64) -> Vec<(String, usize)> {
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

    outliers
}

/// What `determine_kmer_size` actually decided, kept around so the run
/// summary can report it (and, if computed, what mean genome size drove it)
/// without `determine_kmer_size` needing to print anythin itself.
#[derive(Debug, Clone)]
pub struct KmerSelection {
    pub kmer_size: u8,
    pub user_specified: bool,

    /// Only set when `kmer_size` was computed from the dataset
    /// (Fofanov's formula) rather than supplied by the user
    pub mean_genome_size: Option<u32>,
}

pub fn determine_kmer_size(
    sketch_args: &cli::SketchArgs,
    stats: &[(String, usize)],
) -> anyhow::Result<KmerSelection> {
    if let Some(km) = sketch_args.kmer {
        Ok(KmerSelection {
            kmer_size: km,
            user_specified: true,
            mean_genome_size: None,
        })
    } else {
        let mean_genome_size =
            (stats.iter().map(|x| x.1 as u64).sum::<u64>() / stats.len() as u64) as u32;
        let kmer_size = k_computing(mean_genome_size, 0.01);
        Ok(KmerSelection {
            kmer_size,
            user_specified: false,
            mean_genome_size: Some(mean_genome_size),
        })
    }
}

// Validate input files and collect them efficiently
//
// Accepts a mix of file paths and directory paths: directories ate
// recursively expanded to their contents, files are used directly.
// This is the single canonical file list computation and sketching
// consumes to avoid disagreement with the sketching process.
pub fn validate_and_collect_inputs(inputs: &[String]) -> Result<Vec<PathBuf>> {
    let mut candidates: Vec<PathBuf> = Vec::new();

    for input in inputs {
        let path = PathBuf::from(input);
        let metadata =
            fs::metadata(&path).with_context(|| format!("Failed to read path: {}", input))?;

        if metadata.is_dir() {
            let entries: Vec<_> = fs::read_dir(&path)
                .with_context(|| format!("Failed to read directory: {}", input))?
                .collect::<Result<Vec<_>, _>>()?;
            candidates.extend(entries.into_iter().map(|e| e.path()));
        } else {
            candidates.push(path);
        }
    }

    let (valid_files, validation_results): (Vec<_>, Vec<_>) = candidates
        .into_par_iter()
        .map(|path| {
            let is_valid = is_fasta_format(&path);
            let is_multi = if is_valid {
                is_multi_fasta(&path)
            } else {
                false
            };
            (path, (is_valid, is_multi))
        })
        .unzip();

    let invalid_files: Vec<_> = valid_files
        .iter()
        .zip(validation_results.iter())
        .filter_map(|(path, (is_valid, _))| if !is_valid { Some(path) } else { None })
        .collect();

    if !invalid_files.is_empty() {
        anyhow::bail!(
            "Input validation error: Only FASTA files are allowed. Invalid files: {:?}",
            invalid_files
        );
    }

    Ok(valid_files)
}

fn compact_size(size: usize) -> String {
    if size >= 1_000_000_000 {
        format!("{:.1} Gb", size as f64 / 1_000_000_000.0)
    } else if size >= 1_000_000 {
        format!("{:.1} Mb", size as f64 / 1_000_000.0)
    } else if size >= 1_000 {
        format!("{:.1} Kb", size as f64 / 1_000.0)
    } else {
        format!("{} bp", size)
    }
}

fn kv(label: &str, value: impl std::fmt::Display) -> String {
    format!("  {:<width$}{}", label, value, width = LABEL_WIDTH)
}

fn continuation(text: &str) -> String {
    format!("  {:<width$}{}", "", text, width = LABEL_WIDTH)
}

/// Print the whole run's report as one structured summary to stderr —
/// fixed-width, aligned columns under clearly separated sections.
#[allow(clippy::too_many_arguments)]
pub fn print_run_summary(
    stats: &[(String, usize)],
    outliers: &[(String, usize)],
    kmer: &KmerSelection,
    sketch_args: &cli::SketchArgs,
    uncertainty: &cli::UncertaintyArgs,
    estimates: &[DistanceEstimate],
    rescue_summary: Option<&RescueSummary>,
    output_path: &Option<String>,
    verbose: bool,
) {
    let rule = "=".repeat(RULE_WIDTH);
    let sub_rule = "-".repeat(RULE_WIDTH);

    eprintln!("{}", rule);
    eprintln!(" cedar dist - run summary");
    eprintln!("{}", rule);

    eprintln!(" INPUT");
    let sizes: Vec<usize> = stats.iter().map(|(_, s)| *s).collect();
    let (min_size, max_size) = (
        sizes.iter().min().copied().unwrap_or(0),
        sizes.iter().max().copied().unwrap_or(0),
    );
    let mean_size = if stats.is_empty() {
        0
    } else {
        sizes.iter().sum::<usize>() / stats.len()
    };
    eprintln!("{}", kv("genomes", stats.len()));
    eprintln!(
        "{}",
        kv(
            "size range",
            format!(
                "{}-{} (mean {})",
                compact_size(min_size),
                compact_size(max_size),
                compact_size(mean_size)
            )
        )
    );

    if verbose {
        for (id, size) in stats {
            eprintln!("{}", kv(&format!("  {}", id), compact_size(*size)));
        }
    }

    for (id, size) in outliers {
        eprintln!(
            "{}",
            kv("outlier", format!("{} ({})", id, compact_size(*size)))
        );
        if kmer.user_specified {
            eprintln!(
                "{}",
                continuation("=> k set manually; check its pairs below")
            );
        } else {
            eprintln!(
                "{}",
                continuation("=> may have skewed auto k; consider -k, or check its pairs below")
            );
        }
    }

    eprintln!("{}", sub_rule);

    eprintln!(" PARAMETERS");
    let k_note = if kmer.user_specified {
        if kmer.mean_genome_size.is_none() {
            "from pre-computed sketches".to_string()
        } else {
            "user-specified".to_string()
        }
    } else {
        format!(
            "auto, p = 0.01, mean {}",
            compact_size(kmer.mean_genome_size.unwrap_or(0) as usize)
        )
    };
    eprintln!(
        "{}",
        kv("k-mer size", format!("{} ({})", kmer.kmer_size, k_note))
    );
    eprintln!("{}", kv("sketch size", sketch_args.size));
    eprintln!("{}", kv("seed", sketch_args.seed));
    eprintln!(
        "{}",
        kv(
            "confidence",
            format!("{:.0}%", uncertainty.confidence * 100.0)
        )
    );
    eprintln!("{}", kv("rescue p-value", uncertainty.rescue_pvalue));

    eprintln!("{}", sub_rule);

    eprintln!(" DISTANCES");
    let reliable = estimates
        .iter()
        .filter(|e| e.reliability == Reliability::Reliable)
        .count();
    let borderline = estimates
        .iter()
        .filter(|e| e.reliability == Reliability::Borderline)
        .count();
    let unreliable = estimates
        .iter()
        .filter(|e| {
            matches!(
                e.reliability,
                Reliability::Unreliable | Reliability::NoSharedHashes
            )
        })
        .count();
    eprintln!("{}", kv("pairs computed", estimates.len()));
    eprintln!("{}", kv("reliable", reliable));
    eprintln!("{}", kv("borderline", borderline));
    eprintln!("{}", kv("unreliable", unreliable));

    eprintln!("{}", sub_rule);

    eprintln!(" RESCUE");
    match rescue_summary {
        Some(r) => {
            eprintln!("{}", kv("fixable (sketch)", r.needs_larger_sketch));
            eprintln!("{}", kv("unreliable => better", r.unreliable_rescued));
            eprintln!("{}", kv("borderline => reliable", r.borderline_rescued));
            eprintln!("{}", kv("still unresolved", r.still_unresolvable));
        }
        None => eprintln!("{}", kv("status", "disabled (--no-rescue)")),
    }

    eprintln!("{}", sub_rule);

    eprintln!(" OUTPUT");
    match output_path {
        Some(path) => eprintln!("{}", kv("TSV", path)),
        None => eprintln!("{}", kv("TSV", "not written (pass -o/--output)")),
    }
    eprintln!("{}", rule);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_fasta_format_ok() {
        assert!(is_fasta_format(&PathBuf::from("test/bac168.fna")));
    }

    #[test]
    fn test_is_fasta_format_not_ok() {
        assert!(!is_fasta_format(&PathBuf::from("test/test.fq")));
    }

    #[test]
    fn test_format_genome_size() {
        // Below 1 Kb: no approximation suffix.
        assert_eq!(format_genome_size(500), "500 bp ");
        // Kb branch: this used to print "(~true Kb)", a `bool` was being
        // formatted instead of the actual kilobase value.
        assert_eq!(format_genome_size(2_500), "2500 bp (~2.5 Kb)");
        assert_eq!(format_genome_size(3_400_000), "3400000 bp (~3.4 Mb)");
        assert_eq!(format_genome_size(5_200_000_000), "5200000000 bp (~5.2 Gb)");
    }
}
