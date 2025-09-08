// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use anyhow::{anyhow, Context, Result};
use rayon::prelude::*;
use regex::Regex;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use subprocess::Redirection;

use crate::cli;
use crate::utils;

mod bootstrap;
pub mod dist;
pub mod sketch;

/// Configuration for phylogenetic analysis
#[derive(Debug, Clone)]
pub struct PhyloConfig {
    pub threads: usize,
    pub bootstrap_reps: usize,
    pub outlier_threshold: f64,
    pub buffer_size: usize,
}

impl Default for PhyloConfig {
    fn default() -> Self {
        Self {
            threads: 1,
            bootstrap_reps: 1000,
            outlier_threshold: 0.01,
            buffer_size: 8192, // 8KB
        }
    }
}

pub fn build_tree_using_mash_distance(args: &cli::BuildArgs, threads: usize) -> Result<()> {
    let config = PhyloConfig {
        threads,
        bootstrap_reps: args.bootstrap.unwrap_or(1000),
        ..Default::default()
    };

    // Validate and collect input files
    let inputs = validate_and_collect_inputs(&args.indir)?;

    // Compute genome statistics in parallel and genome size outliers
    let stats = utils::compute_genome_stats(&inputs)?;
    utils::check_genome_outliers(&stats, config.outlier_threshold)?;

    // Determine k-mer size
    let kmer_size = utils::determine_kmer_size(&args, &stats);

    // Create sketches
    let sketches = sketch::create_and_load_sketches(args, kmer_size)?;

    // Build tree
    let tree_algorithm = dist::TreeAlgorithm::from_cli(args.canonical, config.threads);

    if let Some(bootstrap_reps) = args.bootstrap {
        println!("Bootstrap replicates: {}", bootstrap_reps);
        bootstrap::build_bootstrap_tree(
            sketches,
            bootstrap_reps,
            &tree_algorithm,
            args.output.clone(),
        )
    } else {
        bootstrap::build_single_tree(sketches, &tree_algorithm, &args)
    }
}

pub fn build_tree_from_orthologous_groups(
    args: &cli::BuildArgs,
    threads: usize,
) -> anyhow::Result<()> {
    let orthologous = &args.indir;
    let outdir = if let Some(p) = &args.output {
        p.to_string()
    } else {
        "cedar_result".to_string()
    };

    // Validate input early
    //utils::validate_inputs(&orthologous)?;

    // Extract species names from headers, inside [...]
    collect_species(&orthologous, &outdir)?;

    // align groups with mafft linsi
    let mut inputs = Vec::new();

    for entry in fs::read_dir(orthologous)? {
        let entry = entry?;
        inputs.push(entry.path());
    }
    let aln_outdir = format!("{}/aln", outdir);
    run_mafft_parallel(&inputs, &aln_outdir)?;

    // trimming alignments with trimal
    let mut trim_input = Vec::new();
    for entry in fs::read_dir(format!("{}/aln", orthologous).as_str())? {
        let entry = entry?;
        trim_input.push(entry.path());
    }
    let trim_outdir = format!("{}/trim", outdir);
    run_trimal(&trim_input, &trim_outdir)?;

    // concatenate alignments into supermatrix
    build_supermatrix(&outdir)?;

    // run iqtree3
    run_iqtree3(
        format!("{outdir}/supermatrix.fa").as_str(),
        if let Some(threads) = args.bootstrap {
            threads
        } else {
            1000
        },
        threads,
    )?;
    Ok(())
}

pub fn run_mafft_parallel(inputs: &[PathBuf], outdir: &str) -> Result<()> {
    println!("[2] Running MAFFT alignments...");
    fs::create_dir_all(outdir)?;

    // Process in parallel with better error collection
    let results: Vec<Result<()>> = inputs
        .par_iter()
        .map(|infile| -> Result<()> {
            let filename = infile
                .file_stem()
                .and_then(|s| s.to_str())
                .ok_or_else(|| anyhow!("Invalid filename: {}", infile.display()))?;

            let outfile = format!("{}/{}.aln", outdir, filename);

            let capture = subprocess::Exec::cmd("mafft")
                .arg("--localpair")
                .arg("--maxiterate")
                .arg("1000")
                .arg(infile)
                .stdout(Redirection::File(File::create(&outfile)?))
                .stderr(Redirection::Pipe)
                .capture()
                .with_context(|| format!("Failed to execute MAFFT for {}", infile.display()))?;

            if !capture.exit_status.success() {
                return Err(anyhow!(
                    "MAFFT failed for {}: {}",
                    infile.display(),
                    String::from_utf8_lossy(&capture.stderr)
                ));
            }

            println!("    Aligned: {}", filename);
            Ok(())
        })
        .collect();

    // Check for any failures
    let failures: Vec<_> = results.into_iter().filter_map(|r| r.err()).collect();
    if !failures.is_empty() {
        return Err(anyhow!("MAFFT failed for {} files", failures.len()));
    }

    Ok(())
}

pub fn run_trimal(inputs: &[PathBuf], outdir: &str) -> anyhow::Result<()> {
    println!("[3] Trimming alignments ({} files)...", inputs.len());

    inputs.par_iter().try_for_each(|infile| -> Result<()> {
        let stem = infile
            .file_stem()
            .and_then(|s| s.to_str())
            .ok_or_else(|| anyhow!("Invalid filename: {}", infile.display()))?;

        let outfile = format!("{}/{}.trim", outdir, stem);

        let capture = subprocess::Exec::cmd("trimal")
            .arg("-in")
            .arg(infile)
            .arg("-out")
            .arg(&outfile)
            .arg("-automated1")
            .stderr(Redirection::Pipe)
            .capture()
            .with_context(|| format!("Failed to execute trimal for {}", stem))?;

        if !capture.exit_status.success() {
            return Err(anyhow!(
                "trimal failed for {}: {}",
                stem,
                String::from_utf8_lossy(&capture.stderr)
            ));
        }

        Ok(())
    })
}

/// Concatenate alignments into a supermatrix and partitions file (parallelized with Rayon)
pub fn build_supermatrix(outdir: &str) -> anyhow::Result<()> {
    println!("[4] Building supermatrix...");

    let species_list_path = Path::new(outdir).join("species_list.txt");
    let species_list = fs::read_to_string(&species_list_path)?
        .lines()
        .map(|s| s.to_string())
        .collect::<Vec<_>>();

    let outdir_path = Path::new(outdir);
    let supermatrix_path = outdir_path.join("supermatrix.fa");
    let partitions_path = outdir_path.join("partitions.txt");
    let trim_dir = outdir_path.join("trim");

    // Collect and sort alignment files
    let mut aln_files: Vec<PathBuf> = fs::read_dir(&trim_dir)?
        .filter_map(|entry| {
            let path = entry.ok()?.path();
            if path.extension().map(|e| e == "trim").unwrap_or(false) {
                Some(path)
            } else {
                None
            }
        })
        .collect();
    aln_files.sort();

    // Process alignments in parallel with better memory management
    type AlignmentData = (String, usize, HashMap<String, String>);

    let results: Vec<AlignmentData> = aln_files
        .into_par_iter()
        .map(|aln_path| -> Result<AlignmentData> { process_single_alignment(&aln_path) })
        .collect::<Result<Vec<_>>>()?;

    // Write outputs efficiently
    write_supermatrix_and_partitions(&supermatrix_path, &partitions_path, &results, &species_list)?;

    println!("    Supermatrix built with {} genes", results.len());
    Ok(())
}

fn process_single_alignment(aln_path: &Path) -> Result<(String, usize, HashMap<String, String>)> {
    let gene_name = aln_path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown")
        .to_string();

    let file = File::open(aln_path)
        .with_context(|| format!("Failed to open alignment file: {}", aln_path.display()))?;

    let reader = BufReader::with_capacity(8192, file);
    let mut sequences: HashMap<String, String> = HashMap::new();
    let mut current_species: Option<String> = None;

    // Parse FASTA more efficiently
    for line in reader.lines() {
        let line = line?;

        if line.starts_with('>') {
            // Extract species from header [species_name]
            if let Some(species) = extract_species_from_header(&line) {
                current_species = Some(species.clone());
                sequences.entry(species).or_insert_with(String::new);
            }
        } else if let Some(species) = &current_species {
            if let Some(seq) = sequences.get_mut(species) {
                seq.push_str(line.trim());
            }
        }
    }

    // Get alignment length from first sequence
    let alignment_length = sequences.values().next().map(|s| s.len()).unwrap_or(0);

    if alignment_length == 0 {
        return Err(anyhow!("Empty alignment in {}", aln_path.display()));
    }

    // Validate alignment lengths
    for (species, seq) in &sequences {
        if seq.len() != alignment_length {
            return Err(anyhow!(
                "Inconsistent alignment length for {} in {}: expected {}, got {}",
                species,
                gene_name,
                alignment_length,
                seq.len()
            ));
        }
    }

    Ok((gene_name, alignment_length, sequences))
}

fn write_supermatrix_and_partitions(
    supermatrix_path: &Path,
    partitions_path: &Path,
    results: &[(String, usize, HashMap<String, String>)],
    species_list: &[String],
) -> Result<()> {
    // Initialize sequence storage
    let mut concatenated_seqs: HashMap<String, String> = species_list
        .iter()
        .map(|sp| (sp.clone(), String::new()))
        .collect();

    // Pre-calculate total length for efficient allocation
    let total_length: usize = results.iter().map(|(_, len, _)| len).sum();
    for seq in concatenated_seqs.values_mut() {
        seq.reserve(total_length);
    }

    // Write partitions file while building sequences
    let mut partitions_writer = BufWriter::new(File::create(partitions_path)?);
    let mut position = 1;

    for (gene_name, alignment_length, gene_sequences) in results {
        // Append sequences for each species
        for species in species_list {
            let pat = "-".repeat(*alignment_length);
            let sequence = gene_sequences
                .get(species)
                .map(|s| s.as_str())
                .unwrap_or(&pat);

            concatenated_seqs
                .get_mut(species)
                .unwrap()
                .push_str(sequence);
        }

        // Write partition entry
        let end_position = position + alignment_length - 1;
        writeln!(
            partitions_writer,
            "{}\t{}-{}",
            gene_name, position, end_position
        )?;
        position = end_position + 1;
    }
    partitions_writer.flush()?;

    // Write supermatrix
    let mut supermatrix_writer = BufWriter::new(File::create(supermatrix_path)?);
    for species in species_list {
        writeln!(supermatrix_writer, ">{}", species)?;
        if let Some(seq) = concatenated_seqs.get(species) {
            writeln!(supermatrix_writer, "{}", seq)?;
        } else {
            return Err(anyhow!("Missing sequence for species: {}", species));
        }
    }
    supermatrix_writer.flush()?;

    Ok(())
}

fn extract_species_from_header(header: &str) -> Option<String> {
    // More efficient species extraction
    if let Some(start) = header.find('[') {
        if let Some(end) = header[start..].find(']') {
            return Some(header[start + 1..start + end].to_string());
        }
    }
    None
}

fn collect_species(indir: &str, outdir: &str) -> Result<Vec<String>> {
    println!("[1] Collecting species names...");

    let species_list_path = Path::new(outdir).join("species_list.txt");
    let header_re = Arc::new(Regex::new(r"^\>.*\[(.+)\]")?);

    // Collect species names in parallel
    let input_files = collect_fasta_files(indir)?;

    let species_sets: Vec<HashSet<String>> = input_files
        .into_par_iter()
        .map(|path| -> Result<HashSet<String>> {
            let file = File::open(&path)?;
            let reader = BufReader::new(file);
            let mut local_taxa = HashSet::new();

            for line in reader.lines() {
                let line = line?;
                if let Some(caps) = header_re.captures(&line) {
                    local_taxa.insert(caps[1].to_string());
                }
            }
            Ok(local_taxa)
        })
        .collect::<Result<_>>()?;

    // Merge all species sets
    let mut all_taxa: HashSet<String> = HashSet::new();
    for species_set in species_sets {
        all_taxa.extend(species_set);
    }

    // Sort and write species list
    let mut taxa: Vec<String> = all_taxa.into_iter().collect();
    taxa.sort();

    let mut out = BufWriter::new(File::create(&species_list_path)?);
    for sp in &taxa {
        writeln!(out, "{}", sp)?;
    }
    out.flush()?;

    println!("    Found {} unique taxa.", taxa.len());
    Ok(taxa)
}

fn run_iqtree3(supermatrix: &str, bootstrap_reps: usize, threads: usize) -> anyhow::Result<()> {
    let outdir = Path::new(supermatrix)
        .parent()
        .and_then(|p| p.to_str())
        .ok_or_else(|| anyhow!("Invalid supermatrix path"))?;

    println!("[5] Running IQ-TREE3...");

    let supermatrix = format!("{}/supermatrix.fa", outdir);
    let partitions = format!("{}/partitions.txt", outdir);

    let capture = subprocess::Exec::cmd("iqtree3")
        .arg("-s")
        .arg(&supermatrix)
        .arg("-spp")
        .arg(&partitions)
        .arg("-m")
        .arg("MFP+MERGE")
        .arg("-B")
        .arg(&bootstrap_reps.to_string())
        .arg("-T")
        .arg(&threads.to_string())
        .arg("--quiet") // Reduce verbose output
        .stderr(Redirection::Pipe)
        .capture()
        .context("Failed to execute IQ-TREE3")?;

    if !capture.exit_status.success() {
        return Err(anyhow!(
            "IQ-TREE3 failed: {}",
            String::from_utf8_lossy(&capture.stderr)
        ));
    }

    println!("    IQ-TREE3 completed successfully");
    Ok(())
}

// Validate input files and collect them efficiently
fn validate_and_collect_inputs(indir: &str) -> Result<Vec<PathBuf>> {
    let entries: Vec<_> = fs::read_dir(indir)
        .with_context(|| format!("Failed to read directory: {}", indir))?
        .collect::<Result<Vec<_>, _>>()?;

    let (valid_files, validation_results): (Vec<_>, Vec<_>) = entries
        .into_par_iter()
        .map(|entry| {
            let path = entry.path();
            let is_valid = utils::is_fasta_format(&path);
            let is_multi = if is_valid {
                utils::is_multi_fasta(&path)
            } else {
                false
            };
            (path, (is_valid, is_multi))
        })
        .unzip();

    // Check for validation errors
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

    let multi_seq_files: Vec<_> = valid_files
        .iter()
        .zip(validation_results.iter())
        .filter_map(|(path, (_, is_multi))| if *is_multi { Some(path) } else { None })
        .collect();

    if !multi_seq_files.is_empty() {
        anyhow::bail!(
            "Input validation error: Multi-sequence FASTA files detected: {:?}",
            multi_seq_files
        );
    }

    Ok(valid_files)
}

fn create_output_structure(outdir: &str) -> Result<()> {
    for subdir in &["aln", "trim"] {
        fs::create_dir_all(format!("{}/{}", outdir, subdir))
            .with_context(|| format!("Failed to create directory: {}/{}", outdir, subdir))?;
    }
    Ok(())
}

fn collect_fasta_files(dir: &str) -> Result<Vec<PathBuf>> {
    let mut files: Vec<PathBuf> = fs::read_dir(dir)?
        .filter_map(|entry| {
            let path = entry.ok()?.path();
            if utils::is_fasta_format(&path) {
                Some(path)
            } else {
                None
            }
        })
        .collect();

    files.sort();
    Ok(files)
}

fn collect_alignment_files(dir: &str) -> Result<Vec<PathBuf>> {
    let mut files: Vec<PathBuf> = fs::read_dir(dir)?
        .filter_map(|entry| {
            let path = entry.ok()?.path();
            if path.extension().map(|e| e == "aln").unwrap_or(false) {
                Some(path)
            } else {
                None
            }
        })
        .collect();

    files.sort();
    Ok(files)
}
