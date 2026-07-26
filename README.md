# cedar

**Uncertainty-aware distance-based phylogenomics from Mash sketches.**

[![CI](https://github.com/Ebedthan/cedar/actions/workflows/ci.yml/badge.svg)](https://github.com/Ebedthan/cedar/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/Ebedthan/cedar/graph/badge.svg?token=S3OLFFHF4X)](https://codecov.io/gh/Ebedthan/cedar)
[![License](https://img.shields.io/badge/license-MIT-blue?style=flat)](https://github.com/Ebedthan/cedar/blob/main/LICENSE)

`cedar` computes MinHash (Mash) distances between genomes and, unlike existing
Mash-based tools, propagates the *statistical uncertainty* of those distances
through to every downstream result, a pairwise comparison, a tree topology,
or a decision about whether two genomes can even be reliably compared at all.

## Why cedar

Mash itself already tells you, in principle, how much to trust a distance
estimate: sketch size bounds the sampling error of the underlying Jaccard
estimate, and Mash's own significance test can tell a real signal from chance
k-mer collisions (Ondov et al., 2016, *Genome Biology*). In practice, almost
nothing downstream of Mash actually uses that information, a distance
matrix is a distance matrix, and a neighbor-joining tree treats every entry
in it as equally trustworthy.

`cedar` doesn't. It computes an **exact confidence interval** on every
pairwise Jaccard/Mash-distance estimate (Clopper–Pearson, not a normal
approximation, the latter is known to be unreliable for exactly the
divergent-genome case where it matters most), classifies each pair's
reliability, and:

- for **`cedar dist`**, attempts to *rescue* a genuinely divergent pair by
  searching for a smaller k-mer size, but only accepts a rescue if it passes
  Mash's own significance test too, not just a tighter-looking interval (a
  small enough k will always eventually produce *some* shared k-mers by
  chance alone; the significance check is what keeps a rescue honest).
- for **`cedar build`**, checks whether unreliable pairs are structurally
  necessary to keep the tree's genomes connected at all (via a minimum
  spanning tree over the resolvable/divergent distinction), flags exactly
  which ones if so, and can search for a smaller dataset-wide k-mer size to
  resolve them (`--include-div-pairs`) rather than silently building a tree
  on data it can't vouch for.

## Installation

Build from source:

```bash
git clone https://github.com/Ebedthan/cedar
cd cedar
cargo build --release
# binary at target/release/cedar
```

Requires Rust 1.89.0 or newer (MSRV, driven by the `statrs` dependency).

## Quick start

```bash
# Build a neighbor-joining tree with bootstrap support
cedar build genomes/ -B 1000 -o tree.nwk

# Compute pairwise distances with uncertainty estimates, no tree
cedar dist genomes/ -o distances.tsv
```

`genomes/` should contain one FASTA file per genome (single-sequence FASTA;
multi-FASTA input is rejected, split concatenated sequences per genome
first).

## `cedar build`

Builds a Mash-distance neighbor-joining consensus tree.

```
cedar build <DIR> [OPTIONS]

Arguments:
  <DIR>  Input directory containing FASTA files

Options:
  -o, --output <FILE>          Output tree (Newick) to FILE [default: stdout]
  -B, --boot <INT>             Bootstrap replicates (conflicts with -J)
  -J, --jack <INT>             Jackknife replicates (conflicts with -B)
      --jacknife-prop <FLOAT>  Subsampling proportion for jackknife [default: 0.5]
      --canonical              Use the canonical NJ solver instead of the
                                default parallel (RapidBTrees) one
      --include-div-pairs      If some genomes can only be connected via
                                divergent (statistically unreliable) pairs,
                                search for a smaller dataset-wide k-mer size
                                that resolves them, rather than proceeding
                                with a tree built partly on unreliable data
  -s, --size <INT>             Sketch size [default: 1000]
  -S, --seed <INT>             Hash seed [default: 42]
  -k, --kmer <INT>             K-mer size [default: computed from mean
                                genome size, targeting a 1% chance-collision
                                probability, Fofanov et al., 2004]
```

With neither `-B` nor `-J`, `cedar build` produces a single point-estimate
tree from the raw distance matrix (no consensus/support values).

Bootstrap and jackknife trees are built from a **majority-rule consensus**
across replicates (a clade must appear in more than half the replicate
trees), with each retained clade's support reported as the internal-node
label in the output Newick, following the RAxML/IQ-TREE convention
(`)95:0.0400`).

### Connectivity checking and `--include-div-pairs`

Before building the tree, `cedar build` checks whether the dataset's
genomes can be connected using only statistically reliable pairwise
distances. If not, if some genome (or cluster of genomes) can only reach
the rest of the dataset through a pair flagged Unreliable or NoSharedHashes, it prints exactly which pair(s) are load-bearing for connectivity:

```
Warning: 1 divergent pair(s) are load-bearing for connectivity and will
still be used in the tree (NJ requires a complete distance matrix, so these
pairs can't simply be dropped). Treat branches resting on them with caution:
  - genomeA vs genomeE
Rerun with --include-div-pairs to attempt resolving them with a smaller
k-mer size instead.
```

With `--include-div-pairs`, cedar searches for the smallest *dataset-wide*
k-mer size (down to a floor set by the largest genome in the dataset, below
which any apparent fix would just be k-mer-space exhaustion rather than
real homology) that resolves the connectivity gap. The whole dataset is
re-sketched at each candidate, not just the problem pairs, since Mash
distance depends directly on k, and mixing k across a single distance
matrix would break neighbor-joining's implicit assumption that every entry
is comparable.

## `cedar dist`

Computes pairwise Mash distances annotated with an uncertainty estimate on
each one, no tree is built. This is for the "is this pair a solid
species-boundary call" moment (e.g. the ~95% ANI / ~0.05 Mash-distance
species threshold), not for producing a phylogeny.

```
cedar dist <DIR> [OPTIONS]

Arguments:
  <DIR>  Input directory containing FASTA files

Options:
  -o, --output <FILE>  Output distance table (TSV) to FILE [default: stdout]
  -s, --size <INT>     Sketch size [default: 1000]
  -S, --seed <INT>     Hash seed [default: 42]
  -k, --kmer <INT>     K-mer size [default: computed from mean genome size]
```

### Output columns

| Column | Meaning |
|---|---|
| `genome1`, `genome2` | Genome identifiers (basename, extension stripped) |
| `jaccard` | Estimated Jaccard index |
| `jaccard_ci_95_low` / `_high` | Exact 95% Clopper–Pearson CI on the Jaccard estimate |
| `mash_distance` | Mash distance, `D = -(1/k)·ln(2J/(1+J))` |
| `mash_distance_ci_95_low` / `_high` | The Jaccard CI transformed directly through the same formula (not a linear approximation) |
| `shared_hashes`, `total_hashes` | Raw counts behind the Jaccard estimate |
| `relative_uncertainty` | Relative half-width of the Jaccard CI — the quantity `flag` is actually classified on (not the distance's own CI, which compresses toward 0 for close relatives and would give a misleading signal) |
| `flag` | `reliable` / `borderline` / `unreliable` / `no_shared_hashes` |
| `kmer_size_used` | K-mer size these values were actually computed at |
| `rescued` | Whether this pair's values came from the k-mer rescue mechanism rather than the run's default k |

### Reliability classification

| Flag | Relative CI half-width on Jaccard |
|---|---|
| `reliable` | ≤ 10% |
| `borderline` | 10–30% |
| `unreliable` | > 30% |
| `no_shared_hashes` | zero shared hashes — no signal to estimate uncertainty from at all |

### The rescue mechanism

For any pair flagged `unreliable` or `no_shared_hashes`, cedar first checks
whether a *larger sketch alone* (up to a ceiling) would bring it back to
`borderline`, if so, it's left alone and reported as-is, since sketch size
governs sampling precision, not whether real biological signal exists, and
resketching the whole dataset larger is a decision for the user (`-s`), not
something to do silently per-pair.

If a larger sketch wouldn't help, cedar re-sketches *only that pair's two
genomes* at progressively smaller k-mer sizes (down to a floor based on the
larger of the two genomes, below which chance k-mer matches become likely
regardless of any real relationship). A candidate k is accepted only if it
passes **both**:

- relative CI half-width at or under the borderline threshold, and
- Mash's own significance test (Ondov et al., Eq. 8): the observed sharing
  must be statistically distinguishable from two unrelated random genomes
  at that k, not merely numerically less noisy.

A rescued pair's row reports the k it was actually resolved at
(`kmer_size_used`) and `rescued=true`, so it's clear that row isn't on the
same k-mer scale as the rest of the table.

## Statistical methodology

- **Jaccard confidence intervals**: exact Clopper–Pearson, computed via the
  regularized incomplete beta function (treating the shared-hash count as
  `Binomial(s, J)`, the same hypergeometric-to-binomial approximation Mash
  itself uses), rather than a normal/Wald approximation.
- **Mash distance confidence intervals**: the Jaccard CI's bounds
  transformed *directly* through the Mash distance formula, not linearized
  via the delta method.
- **Significance testing**: Mash's own Eq. 8 (binomial P-value against an
  expected chance-collision rate derived from each genome's size and the
  k-mer size), reusing the same incomplete-beta machinery as the CI.
- **Connectivity analysis**: Kruskal's algorithm over the complete pairwise
  graph, with divergent pairs sentinel-weighted so resolvable edges are
  always preferred — any divergent edge that ends up in the minimum
  spanning tree is one the algorithm was structurally forced to use.

Reference: Ondov, B.D., Treangen, T.J., Melsted, P. et al. *Mash: fast
genome and metagenome distance estimation using MinHash.* Genome Biology
17, 132 (2016). https://doi.org/10.1186/s13059-016-0997-x

## Limitations

- `cedar build` does not perform per-pair k-mer rescue, only a
  dataset-wide search via `--include-div-pairs`, since mixing k-mer sizes
  within one distance matrix would violate neighbor-joining's assumption
  that every entry is comparable. `cedar dist`'s per-pair rescue has no
  such constraint, since it produces an unstructured table, not a tree.
- Multi-FASTA input files are rejected; each genome must be a single-record
  FASTA file.
- Designed for species/strain-level resolution. Datasets requiring
  divergent-pair rescue across most of a tree's backbone (rather than a
  handful of edges) are likely outside what Mash-based distance estimation
  can meaningfully resolve at all, regardless of k or sketch size.

## Development

```bash
cargo test          # requires committed fixtures under test/
cargo clippy --all-targets -- -D warnings
cargo fmt --all -- --check
```

### Licence
`cedar` is distributed under the terms of the MIT license.
See [LICENSE-MIT](https://github.com/Ebedthan/cedar/blob/main/LICENSE-MIT) for details.
