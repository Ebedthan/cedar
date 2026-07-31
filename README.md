# cedar

**Uncertainty-aware distance-based phylogenomics from Mash sketches.**

`cedar` computes pairwise Mash distances between genomes and propagates statistical uncertainty through every result, reporting exact Clopper–Pearson confidence intervals on each distance and flagging estimates that cannot be trusted at the requested confidence level. It is the to the best of our knowledge, the first distance-based phylogenomic tool to tell you how much to trust its own output.

## Features

- Exact 95% (or user-specified) Clopper–Pearson confidence intervals on every pairwise Jaccard and Mash distance
- Per-pair reliability tiers: **reliable**, **borderline**, **unreliable**, **no shared hashes**
- K-mer rescue mechanism: automatically searches for a smaller k-mer size that resolves borderline and unreliable pairs, subject to dual statistical checks (CI width + Mash significance test)
- Connectivity analysis: detects genomes that can only be connected via divergent pairs before building a tree, and optionally resolves them with a global k-mer search (`--include-div-pairs`)
- Bootstrap and jackknife consensus trees with branch support values
- Pre-computed sketch reuse: sketch once, run many times with `cedar sketch` + `--from-sketches`
- Parallelised throughout with rayon

## Installation

```bash
cargo install --path .
```

Minimum supported Rust version: **1.89.0**

## Quick start

```bash
# All in one pass (sketch + distances)
cedar dist genomes/*.fna -o distances.tsv

# Build a neighbour-joining tree with 100 bootstrap replicates
cedar build genomes/*.fna -B 100 -o tree.nwk

# Sketch once, reuse many times
cedar sketch genomes/*.fna -o sketches/
cedar dist  --from-sketches sketches/ -o distances.tsv
cedar build --from-sketches sketches/ -B 100 -o tree.nwk
```

## Subcommands

### `cedar sketch`

Sketch FASTA files to `.msh` files for later reuse. Accepts individual files, directories, or a mix of both.

```
cedar sketch [OPTIONS] <PATH>...
```

| Option | Default | Description |
|---|---|---|
| `-o, --output-dir <DIR>` | `cedar_sketches` | Output directory for `.msh` files |
| `-s, --size <INT>` | 1000 | Sketch size (number of hashes) |
| `-S, --seed <INT>` | 42 | Hash seed |
| `-k, --kmer <INT>` | auto | K-mer size (auto-computed from genome size if not set) |

```bash
cedar sketch genomes/ -o my_sketches/ -s 2000
cedar sketch file1.fna file2.fna dir/ -o my_sketches/
```

### `cedar dist`

Compute pairwise Mash distances with uncertainty estimates. Always prints a structured summary to stderr; the full distance table is only written when `-o` is given.

```
cedar dist [OPTIONS] <PATH>...
cedar dist --from-sketches <DIR> [OPTIONS]
```

| Option | Default | Description |
|---|---|---|
| `-o, --output <FILE>` | — | Write TSV distance table to FILE |
| `--from-sketches <DIR>` | — | Load pre-computed `.msh` files instead of re-sketching |
| `--confidence <FLOAT>` | 0.95 | Confidence level for the Jaccard/Mash-distance interval |
| `--rescue-pvalue <FLOAT>` | 0.01 | P-value threshold for the k-mer rescue significance test |
| `--no-rescue` | off | Disable k-mer rescue; report every pair at the run's single k |
| `-s, --size <INT>` | 1000 | Sketch size |
| `-S, --seed <INT>` | 42 | Hash seed |
| `-k, --kmer <INT>` | auto | K-mer size |
| `-t <INT>` | 1 | Threads |
| `-v, --verbose` | off | Print per-genome details in the run summary |

```bash
cedar dist genomes/ -o distances.tsv
cedar dist --from-sketches sketches/ --confidence 0.99 -o distances.tsv
cedar dist file1.fna file2.fna file3.fna
cedar dist genomes/*.fna --no-rescue -o raw_distances.tsv
```

Both `-v/--verbose` and `-t/--threads` are global flags and can appear before or after the subcommand.

#### TSV output format

The output file begins with a metadata comment line followed by a tab-separated table:

```
# confidence=0.95 rescue_pvalue=0.01
genome1	genome2	jaccard	jaccard_ci_low	jaccard_ci_high	mash_distance	mash_distance_ci_low	mash_distance_ci_high	shared_hashes	total_hashes	relative_uncertainty	flag	kmer_size_used	rescued
```

| Column | Description |
|---|---|
| `jaccard` | Jaccard index estimate |
| `jaccard_ci_low/high` | Exact Clopper–Pearson CI bounds on Jaccard |
| `mash_distance` | Mash distance D = −(1/k) ln(2J/(1+J)) |
| `mash_distance_ci_low/high` | CI bounds transformed through the same formula |
| `relative_uncertainty` | CI half-width / Jaccard estimate (the quantity reliability is classified on) |
| `flag` | `reliable` / `borderline` / `unreliable` / `no_shared_hashes` |
| `kmer_size_used` | The k-mer size actually used for this pair (may differ from the run k if rescued) |
| `rescued` | `true` if this pair was resolved by the k-mer rescue mechanism |

#### Stderr summary

A structured four-section report is always printed to stderr regardless of whether `-o` is given:

```
================================================================
 cedar dist - run summary
================================================================
 INPUT
  genomes                 26
  size range              4.6 Mb-5.5 Mb (mean 5.0 Mb)
----------------------------------------------------------------
 PARAMETERS
  k-mer size              21 (auto, p=0.01, mean 5.0 Mb)
  sketch size             1000
  seed                    42
  confidence              95%
  rescue p-value          0.01
----------------------------------------------------------------
 DISTANCES
  pairs computed          325
  reliable                298
  borderline              21
  unreliable              6
----------------------------------------------------------------
 RESCUE
  fixable (sketch)        4
  unreliable => better    2
  borderline => reliable  8
  still unresolved        4
----------------------------------------------------------------
 OUTPUT
  TSV                     distances.tsv
================================================================
```

### `cedar build`

Build a neighbour-joining consensus tree with optional bootstrap or jackknife support.

```
cedar build [OPTIONS] <PATH>...
cedar build --from-sketches <DIR> [OPTIONS]
```

| Option | Default | Description |
|---|---|---|
| `-o, --output <FILE>` | — | Output Newick tree |
| `--from-sketches <DIR>` | — | Load pre-computed `.msh` files instead of re-sketching. Cannot be combined with `--include-div-pairs`. |
| `-B, --boot <INT>` | — | Bootstrap replicates |
| `-J, --jack <INT>` | — | Jackknife replicates |
| `--jacknife-prop <FLOAT>` | 0.5 | Proportion of hashes to keep per jackknife replicate |
| `--canonical` | off | Compute canonical NJ tree |
| `--include-div-pairs` | off | If divergent pairs are structurally load-bearing, search for a smaller global k that resolves them |
| `--confidence <FLOAT>` | 0.95 | Confidence level |
| `--rescue-pvalue <FLOAT>` | 0.01 | Divergence significance threshold |
| `-s, --size <INT>` | 1000 | Sketch size |
| `-S, --seed <INT>` | 42 | Hash seed |
| `-k, --kmer <INT>` | auto | K-mer size |
| `-t <INT>` | 1 | Threads |
| `-v, --verbose` | off | Verbose stderr summary |

```bash
cedar build genomes/ -B 100 -o tree.nwk
cedar build --from-sketches sketches/ -B 1000 -o tree.nwk
cedar build genomes/ --include-div-pairs -o tree.nwk
cedar build genomes/ -J 100 --jacknife-prop 0.7 -o tree.nwk
```

## Statistical methodology

### K-mer size selection

When `-k` is not specified, cedar computes the smallest k that limits the probability of any k-mer appearing by chance in a random genome of the observed mean size to ≤ 0.01 (Fofanov's formula, as used in the original Mash paper). This ensures the MinHash sketch captures real sequence signal rather than background noise.

### Confidence intervals (Clopper–Pearson)

The Jaccard estimate `Ĵ = shared_hashes / total_hashes` is a sample proportion from a Binomial(s, J) model (the same approximation Mash uses). Cedar computes an exact Clopper–Pearson interval rather than a Wald approximation, using the relationship:

```
P(X ≥ x | n, p) = I_p(x, n−x+1)
```

where I_p is the regularised incomplete beta function. This interval never under-covers (it is conservative rather than anti-conservative), and is implemented via `statrs::distribution::Beta::inverse_cdf`.

The CI is computed on the Jaccard index, and both bounds are then transformed through D = −(1/k) ln(2J/(1+J)) to obtain CI bounds on the Mash distance. This avoids the approximation errors that would result from applying the delta method to D directly.

### Reliability tiers

A pair is classified by the relative CI half-width on J (half-width / Ĵ):

| Tier | Threshold | Meaning |
|---|---|---|
| `reliable` | ≤ 10% | CI is tight relative to the point estimate; trustworthy |
| `borderline` | ≤ 30% | Moderate uncertainty; usable with caution |
| `unreliable` | > 30% | Wide interval; point estimate alone is not a reliable guide |
| `no_shared_hashes` | — | No shared hashes at all; no signal to put an interval on |

These thresholds apply independently of the `--confidence` level. Requesting a 99% interval widens the interval (making the comparison more stringent), correctly increasing the fraction of pairs classified as borderline or unreliable.

### K-mer rescue mechanism

When rescue is enabled (the default; disable with `--no-rescue`), cedar attempts to improve pairs whose reliability falls below their tier's target:

- **Borderline** pairs target the *reliable* threshold (≤ 10%)
- **Unreliable / no-shared-hashes** pairs target the *borderline* threshold (≤ 30%)

Before attempting a k-mer search, cedar checks whether increasing the sketch size alone (up to 100,000 hashes) would already achieve the target — if so, it flags the pair as `fixable (sketch)` rather than attempting rescue, since sketch size governs sampling precision while k governs what sequence is captured.

For pairs where a larger sketch would not help, cedar descends through candidate k-mer sizes (step 2, down to a genome-size-derived floor from Fofanov's formula) and re-sketches only that pair in a temporary directory at each candidate k. A candidate is accepted only if it passes **both** checks:

1. Relative CI half-width ≤ target threshold
2. Mash's own significance test (Eq. 8, P-value ≤ `--rescue-pvalue`): the observed sharing must be distinguishable from what two random unrelated genomes would share by chance

Requiring both checks prevents accepting a rescue that simply reduced k until chance k-mer collisions started dominating the signal.

### Connectivity analysis (`cedar build`)

Before building a NJ tree, cedar runs Kruskal's algorithm over the complete pairwise graph to check whether any genome can only be reached via divergent (unreliable or non-significant) pairs. A divergent pair is structurally load-bearing if it ends up in the minimum spanning tree when resolvable pairs are exhausted first. This warns you before the tree is built rather than after, so you can decide whether to proceed or use `--include-div-pairs`.

`--include-div-pairs` searches for a smaller dataset-wide k that resolves connectivity without divergent edges, re-sketching the entire dataset at each candidate. This is the most expensive operation in cedar: if load-bearing divergent pairs exist, expect this to take substantially longer than a standard `cedar build`.

## Sketch reuse and `--from-sketches`

Pre-computing sketches with `cedar sketch` and reusing them with `--from-sketches` is the recommended workflow for any scenario where:

- The same genome set is compared at multiple confidence levels or k-mer sizes
- Multiple benchmark or analysis runs are performed on the same data
- Disk space is limited and raw FASTAs are deleted after sketching

All sketches in a `--from-sketches` directory must have been created with the same `-s`, `-k`, and `-S` values. Cedar validates this at load time and reports a clear error if parameters are inconsistent.

`--from-sketches` cannot be combined with `--include-div-pairs` on `cedar build`, because the global k-mer search requires re-sketching the entire dataset at candidate k values, which is impossible when the only available input is a fixed set of pre-computed sketches.

## Output notes

- All summaries and warnings are written to **stderr**. Only the TSV distance table (via `-o`) and the Newick tree (via `-o`) go to stdout or files.
- The TSV metadata comment line (`# confidence=... rescue_pvalue=...`) can be skipped in most tools: `pandas.read_csv(..., comment='#')`, `awk 'NR>2'`, `grep -v '^#'`.
- Genome size outliers (dramatically smaller or larger than the rest of the dataset) are flagged as warnings in the run summary. When k is auto-computed, this note indicates the outlier may have skewed the dataset-wide k selection; when k is user-specified, it is purely informational.

## Dependencies

| Crate | Version | Purpose |
|---|---|---|
| `clap` | 4 | CLI |
| `finch` | — | MinHash sketching and `.msh` I/O |
| `statrs` | 0.19 | Beta distribution (Clopper–Pearson CI) |
| `rayon` | 1 | Parallelism |
| `tempfile` | 3.10 | Temporary directories for rescue re-sketching |
| `rand` | 0.10 | Bootstrap/jackknife sampling |
| `anyhow` | — | Error handling |

## Citing `cedar`

Incoming `cedar` article. 

But, if you use `cedar` in your work, please also cite:

> Ondov, B.D. et al. (2016). Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biology*, 17, 132.

The distance computation and significance test are taken directly from that paper.
