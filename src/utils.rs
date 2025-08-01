use crate::dist;
use std::fs::{self, File};
use std::io::BufRead;
use std::io::BufReader;
use std::io::{self, Write};

pub fn compute_newick_tree(
    matrix: &speedytree::DistanceMatrix,
    is_canonical: bool,
    num_threads: usize,
) -> anyhow::Result<String> {
    if is_canonical {
        let tree =
            speedytree::NeighborJoiningSolver::<speedytree::Canonical>::default(matrix.clone())
                .solve()
                .unwrap();
        Ok(speedytree::to_newick(&tree))
    } else {
        let tree =
            speedytree::NeighborJoiningSolver::<speedytree::RapidBtrees>::default(matrix.clone())
                .set_chunk_size(std::cmp::max(matrix.size() / num_threads, 1))
                .solve()
                .unwrap();
        Ok(speedytree::to_newick(&tree))
    }
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
) -> anyhow::Result<()> {
    if keep {
        dist::to_phylip(matrix.clone(), tempdir)
    } else {
        fs::remove_dir_all(tempdir)?;
        Ok(())
    }
}

pub fn is_fasta_format(path: &str) -> bool {
    let file = File::open(path).expect("file should exists before opening.");
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    while let Some(Ok(line)) = lines.next() {
        let trimmed_line = line.trim();
        if trimmed_line.is_empty() {
            continue;
        }
        return match trimmed_line.chars().next() {
            Some('>') => true,
            _ => false,
        };
    }
    false
}

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

// Return sequence id with its length
pub fn get_seq_stats(path: &str) -> anyhow::Result<(String, usize)> {
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

        if trimmed.starts_with('>') {
            let tmp_id = trimmed[1..]
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

pub fn detect_outliers(
    data: &Vec<(String, usize)>,
    epsilon: f64,
) -> anyhow::Result<Vec<(String, usize)>> {
    let mut outliers = Vec::new();
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
            .cloned()
            .filter(|x| x.1 < lower_bound || x.1 > upper_bound)
            .collect();
    }

    Ok(outliers)
}

pub fn validate_inputs(filenames: &[String]) -> anyhow::Result<()> {
    if filenames.len() < 3 {
        anyhow::bail!("At least three input FASTA files must be provided.");
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
            "Only FASTA files are allowed. Invalid files: {}",
            invalid.join(", ")
        );
    }

    if !multi_seq.is_empty() {
        anyhow::bail!(
            "Only single-sequence FASTA files are allowed. Multi-sequence files: {}",
            multi_seq.join(", ")
        );
    }

    Ok(())
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
