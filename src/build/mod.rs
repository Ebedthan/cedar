// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use anyhow::anyhow;
use rayon::prelude::*;
use regex::Regex;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::path::PathBuf;
use subprocess::Redirection;

use crate::cli;
use crate::utils;

mod bootstrap;
pub mod dist;
pub mod sketch;

pub fn build_tree_from_genomes(args: &cli::BuildArgs, threads: usize) -> anyhow::Result<()> {
    let genomes = &args.indir;
    let mut inputs = Vec::new();

    // Validate inputs early
    let mut invalid = vec![];
    let mut multi_seq = vec![];

    for file in fs::read_dir(genomes)? {
        let file = file?;
        inputs.push(file.path());
        if !utils::is_fasta_format(file.path()) {
            invalid.push(file.path());
        } else if utils::is_multi_fasta(file.path()) {
            multi_seq.push(file.path());
        }
    }

    if !invalid.is_empty() {
        anyhow::bail!(
            "Input validation error: Only FASTA files are allowed. Invalid files: {:?}",
            invalid
        );
    }

    // Compute genome statistics in parallel
    let stats = utils::compute_genome_stats(&inputs)?;

    // Check for genome outliers and exit early if found
    utils::check_genome_outliers(&stats, 0.01)?;

    // Determine k-mer size
    let kmer_size = utils::determine_kmer_size(&args, &stats);

    // Set bootstrap replicates
    let reps = args.bootstrap.unwrap_or(100);
    if args.bootstrap.is_some() {
        println!("Boostrap replicated: {}", reps);
    }

    // Create temp directory
    //fs::create_dir_all(&args.tempdir)
    //    .with_context(|| format!("Could not create temp directory: {}", &args.tempdir))?;

    // Create sketches
    let sketches = sketch::create_and_load_sketches(&args, kmer_size)?;

    let tree_algorithm = dist::TreeAlgorithm::from_cli(args.canonical, threads);

    if let Some(bootstrap_reps) = args.bootstrap {
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
    run_mafft_parallel(inputs, &aln_outdir)?;

    // trimming alignments with trimal
    let mut trim_input = Vec::new();
    for entry in fs::read_dir(format!("{}/aln", orthologous).as_str())? {
        let entry = entry?;
        trim_input.push(entry.path());
    }
    let trim_outdir = format!("{}/trim", outdir);
    run_trimal(trim_input, &trim_outdir)?;

    // concatenate alignments into supermatrix
    build_supermatrix(&outdir)?;

    // run iqtree3
    run_iqtree3(
        format!("{outdir}/supermatrix.fa").as_str(),
        format!("{outdir}/partitions.txt").as_str(),
        if let Some(threads) = args.bootstrap {
            threads
        } else {
            1000
        },
        threads,
    )?;
    Ok(())
}

pub fn run_mafft_parallel(inputs: Vec<PathBuf>, outdir: &str) -> anyhow::Result<()> {
    fn helper(infile: &PathBuf, outfile: &str) -> anyhow::Result<()> {
        let capture = subprocess::Exec::cmd("linsi")
            .arg(&infile)
            .stdout(Redirection::File(File::create(outfile)?))
            .stderr(Redirection::Pipe)
            .capture()?;

        if capture.exit_status.success() {
            Ok(())
        } else {
            Err(anyhow::anyhow!(
                "linsi failed for {}: {}",
                &infile.display(),
                String::from_utf8_lossy(&capture.stderr)
            ))
        }
    }

    // create output directory if missing
    std::fs::create_dir_all(outdir)?;

    // parallel iteration over input files
    inputs.par_iter().try_for_each(|infile| {
        let outfile = format!("{}/{}.aln.fasta", outdir, infile.display());
        helper(&infile, &outfile)
    })
}

pub fn run_trimal(inputs: Vec<PathBuf>, outdir: &str) -> anyhow::Result<()> {
    fn helper(infile: &PathBuf, outdir: &str) -> anyhow::Result<()> {
        // ensure output dir exists
        fs::create_dir_all(outdir)?;

        let stem = Path::new(infile)
            .file_stem()
            .ok_or_else(|| anyhow::anyhow!("Invalid filename: {}", infile.display()))?
            .to_string_lossy();

        let outfile = format!("{}/{}.trim", outdir, stem);

        let capture = subprocess::Exec::cmd("trimal")
            .arg("-in")
            .arg(infile)
            .arg("-out")
            .arg(&outfile)
            .arg("-automated1")
            .stderr(Redirection::Pipe)
            .capture()?;

        if capture.exit_status.success() {
            Ok(())
        } else {
            Err(anyhow::anyhow!(
                "trimal failed for {}: {}",
                infile.display(),
                String::from_utf8_lossy(&capture.stderr)
            ))
        }
    }

    std::fs::create_dir_all(outdir)?;

    inputs
        .par_iter()
        .try_for_each(|infile| helper(infile, outdir))
}

/// Concatenate alignments into a supermatrix and partitions file (parallelized with Rayon)
pub fn build_supermatrix(outdir: &str) -> anyhow::Result<()> {
    let outdir = Path::new(outdir);
    let supermatrix_path = outdir.join("supermatrix.fa");
    let partitions_path = outdir.join("partitions.txt");
    let species_list_path = outdir.join("species_list.txt");
    let trim_dir = outdir.join("trim");

    // Read species list
    let species_file = File::open(&species_list_path)?;
    let reader = BufReader::new(species_file);
    let species_list: Vec<String> = reader.lines().collect::<Result<_, _>>()?;

    // Collect alignment files
    let mut aln_files: Vec<PathBuf> = fs::read_dir(&trim_dir)?
        .filter_map(|entry| {
            let p = entry.ok()?.path();
            if p.extension().map(|e| e == "trim").unwrap_or(false) {
                Some(p)
            } else {
                None
            }
        })
        .collect();
    aln_files.sort();

    // Process alignments in parallel
    let results: Vec<(String, usize, HashMap<String, String>)> = aln_files
        .into_par_iter()
        .map(|aln| -> anyhow::Result<_> {
            let gene = aln
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("unknown")
                .to_string();

            println!("  Processing {gene}...");

            let file = File::open(&aln)?;
            let reader = BufReader::new(file);

            // Parse FASTA
            let mut groupseq: HashMap<String, String> = HashMap::new();
            let mut cur_sp: Option<String> = None;

            for line in reader.lines() {
                let line = line?;
                if line.starts_with('>') {
                    // extract species inside [...]
                    if let Some(start) = line.find('[') {
                        if let Some(end) = line.find(']') {
                            let sp = line[start + 1..end].to_string();
                            cur_sp = Some(sp);
                            groupseq.entry(cur_sp.clone().unwrap()).or_default();
                        }
                    }
                } else if let Some(sp) = &cur_sp {
                    groupseq
                        .entry(sp.clone())
                        .and_modify(|seq| seq.push_str(&line));
                }
            }

            // alignment length = first sequence length
            let alen = groupseq
                .values()
                .next()
                .ok_or_else(|| anyhow!("Empty alignment in {}", aln.display()))?
                .len();

            Ok((gene, alen, groupseq))
        })
        .collect::<anyhow::Result<_>>()?;

    // Merge results sequentially (to preserve partition order)
    let mut seqs: HashMap<String, String> = species_list
        .iter()
        .map(|sp| (sp.clone(), String::new()))
        .collect();

    let mut partitions = File::create(&partitions_path)?;
    let mut pos: usize = 1;

    for (gene, alen, groupseq) in results {
        // append sequences
        for sp in &species_list {
            if let Some(s) = groupseq.get(sp) {
                seqs.entry(sp.clone()).and_modify(|seq| seq.push_str(s));
            } else {
                let gaps = "-".repeat(alen);
                seqs.entry(sp.clone()).and_modify(|seq| seq.push_str(&gaps));
            }
        }

        // partitions record
        let end = pos + alen - 1;
        writeln!(partitions, "{}\t{}-{}", gene, pos, end)?;
        pos = end + 1;
    }

    // Write supermatrix
    let mut supermatrix = File::create(&supermatrix_path)?;
    for sp in &species_list {
        writeln!(supermatrix, ">{}", sp)?;
        if let Some(seq) = seqs.get(sp) {
            writeln!(supermatrix, "{}", seq)?;
        } else {
            return Err(anyhow!("Missing sequence for species: {}", sp));
        }
    }

    Ok(())
}

fn collect_species(indir: &str, outdir: &str) -> anyhow::Result<()> {
    println!("[1] Collecting species names...");

    let indir = Path::new(indir);
    let outdir = Path::new(outdir);
    let species_list_path = outdir.join("species_list.txt");

    let header_re = Regex::new(r"^\>.*\[(.+)\]")?;

    let mut taxa: HashSet<String> = HashSet::new();

    // iterate over all .fa files
    for entry in fs::read_dir(indir)? {
        let path = entry?.path();
        if path.extension().map(|e| e == "fa").unwrap_or(false) {
            let file = File::open(&path)?;
            let reader = BufReader::new(file);

            for line in reader.lines() {
                let line = line?;
                if let Some(caps) = header_re.captures(&line) {
                    taxa.insert(caps[1].to_string());
                }
            }
        }
    }

    // sort and write species list
    let mut taxa: Vec<_> = taxa.into_iter().collect();
    taxa.sort();

    let mut out = File::create(&species_list_path)?;
    for sp in &taxa {
        writeln!(out, "{}", sp)?;
    }

    let ntaxa = taxa.len();
    println!("    Found {} unique taxa.", ntaxa);

    Ok(())
}

fn run_iqtree3(
    supermatrix: &str,
    partitions: &str,
    reps: usize,
    threads: usize,
) -> anyhow::Result<()> {
    let capture = subprocess::Exec::cmd("iqtree3")
        .arg("-s")
        .arg(supermatrix)
        .arg("-spp")
        .arg(partitions)
        .arg("-m")
        .arg("MFP+MERGE")
        .arg("-B")
        .arg(format!("{reps}").as_str())
        .arg("-T")
        .arg(format!("{threads}").as_str())
        .stderr(Redirection::Pipe)
        .capture()?;

    if capture.exit_status.success() {
        Ok(())
    } else {
        Err(anyhow::anyhow!(
            "iqtree3 failed: {}",
            String::from_utf8_lossy(&capture.stderr)
        ))
    }
}
