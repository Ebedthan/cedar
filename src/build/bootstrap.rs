use crate::{
    build::matrix::{distance_to_matrix, ComputeTree, TreeAlgorithm},
    cli,
    nwk::{build_tree_from_clades, Clade, Node, Tree},
    utils,
};
use finch::serialization::Sketch;
use rand::{rng, seq::IndexedRandom};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

use crate::mash::distance::compute_distances;

/// Per-clade bookkeeping across replicate trees: how many trees contained
/// this clade (bipartition), and the branch lengths observed for it. These
/// are independent quantities, a clade's support is its occurrence
/// frequency, its length is a separate (averaged) measurement, and must
/// never be conflated into a single number.
#[derive(Default)]
struct CladeStats {
    count: usize,
    lengths: Vec<f64>,
}

/// Majority rule: keep clades occurring in at least `threshold` of the
/// `n_trees` input trees, recording each kept clade's occurrence frequency
/// as its support (e.g. bootstrap support) separately from its branch length.
fn majority_rule(clades: Vec<Clade>, threshold: usize, n_trees: usize) -> Vec<Clade> {
    let mut clade_group: HashMap<Vec<String>, CladeStats> = HashMap::with_capacity(clades.len());

    // single pass: group clades by their leaf sets, counting occurrences and
    // collecting branch lengths independently
    for clade in clades {
        let stats = clade_group.entry(clade.leaves).or_default();
        stats.count += 1;
        if let Some(length) = clade.length {
            stats.lengths.push(length);
        }
    }

    // filter by occurrence count, average lengths, and compute support
    clade_group
        .into_iter()
        .filter(|(_, stats)| stats.count >= threshold)
        .map(|(leaves, stats)| {
            let avg_length = if stats.lengths.is_empty() {
                None
            } else {
                Some(stats.lengths.iter().sum::<f64>() / stats.lengths.len() as f64)
            };
            let support = if n_trees == 0 {
                None
            } else {
                Some(stats.count as f64 / n_trees as f64)
            };
            Clade::new_with_support(leaves, avg_length, support)
        })
        .collect()
}

/// Two clades are compatible (can coexist in the same tree) only if their
/// leaf sets are nested (one a subset of the other) or disjoint. Clades that
/// partially overlap cannot both appear in a valid tree — accepting both
/// would require duplicating a leaf under two different branches.
fn is_compatible(a: &Clade, b: &Clade) -> bool {
    let a_set: HashSet<&String> = a.leaves.iter().collect();
    let b_set: HashSet<&String> = b.leaves.iter().collect();
    let shared = a_set.intersection(&b_set).count();
    shared == 0 || shared == a_set.len() || shared == b_set.len()
}

/// Greedily reduce a set of clades to a compatible (laminar) family, required
/// for building a valid consensus tree. Clades are considered in descending
/// order of support (ties broken deterministically by leaf set) so that
/// higher-support bipartitions win over conflicting lower-support ones.
fn filter_compatible_clades(mut clades: Vec<Clade>) -> Vec<Clade> {
    clades.sort_by(|a, b| {
        b.support
            .partial_cmp(&a.support)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.leaves.cmp(&b.leaves))
    });

    let mut accepted: Vec<Clade> = Vec::with_capacity(clades.len());
    for clade in clades {
        if accepted.iter().all(|kept| is_compatible(&clade, kept)) {
            accepted.push(clade);
        }
    }
    accepted
}

/// Extract clades from a tree, returning Vec<Clade> with sorted leaves
fn get_clades(node: &Node) -> Vec<Clade> {
    let mut clades = Vec::new();

    // If this is a leaf, record a singleton clade with its branch length
    if node.children.is_empty() {
        if let Some(name) = &node.name {
            clades.push(Clade::new(vec![name.clone()], node.length));
        }
        return clades;
    }

    // Otherwise: recurse into children
    let mut leaves = Vec::new();
    for child in &node.children {
        let child_clades = get_clades(child);
        for cl in child_clades {
            clades.push(cl);
        }
        leaves.extend(get_leaf_names(child));
    }

    leaves.sort();
    clades.push(Clade::new(leaves, node.length));

    clades
}

pub fn build_consensus(trees: Vec<Tree>, threshold: usize) -> Tree {
    let n_trees = trees.len();
    let mut clades: Vec<Clade> = Vec::new();

    for tree in trees {
        clades.extend(get_clades(&tree.root));
    }

    let mut all_leaves: Vec<String> = clades
        .iter()
        .flat_map(|c| c.leaves.clone())
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();
    all_leaves.sort();

    let kept_clades = majority_rule(clades, threshold, n_trees);
    let compatible_clades = filter_compatible_clades(kept_clades);
    build_tree_from_clades(&all_leaves, &compatible_clades)
}

// Helper functions
fn get_leaf_names(node: &Node) -> Vec<String> {
    if node.children.is_empty() {
        return vec![node.name.clone().unwrap_or_default()];
    }

    let mut leaves = Vec::new();
    for child in &node.children {
        leaves.extend(get_leaf_names(child));
    }
    leaves
}

/// Sample each sketch's own hashes with replacement, resampling every
/// sketch to its own current size.
/// This is the standard nonparametric bootstrap resampling method.
pub fn sample_sketches_with_replacement(sketches: Vec<Sketch>) -> Vec<Sketch> {
    let mut rng = rng();

    sketches
        .into_iter()
        .map(|sketch| {
            let n = sketch.hashes.len();
            let sampled: Vec<_> = (0..n)
                .map(|_| sketch.hashes.choose(&mut rng).unwrap().clone())
                .collect();

            Sketch {
                hashes: sampled,
                ..sketch
            }
        })
        .collect()
}

/// Subsample each sketch's hashes *without* replacement, keeping a fixed
/// proportion of them, the jackknife counterpart to
/// `sample_sketches_with_replacement`'s bootstrap resampling. A sketch never
/// has fewer than 1 hash kept, even if `proportion` would round down to 0.
pub fn sample_sketches_without_replacement(sketches: Vec<Sketch>, proportion: f64) -> Vec<Sketch> {
    let mut rng = rng();

    sketches
        .into_iter()
        .map(|sketch| {
            let keep = ((sketch.hashes.len() as f64) * proportion).round() as usize;
            let keep = keep.clamp(1, sketch.hashes.len().max(1));
            let sampled: Vec<_> = sketch.hashes.sample(&mut rng, keep).cloned().collect();

            Sketch {
                hashes: sampled,
                ..sketch
            }
        })
        .collect()
}

pub fn build_bootstrap_tree(
    sketches: Vec<Sketch>,
    reps: usize,
    tree_algorithm: &TreeAlgorithm,
    output: Option<String>,
) -> anyhow::Result<()> {
    let trees: anyhow::Result<Vec<Tree>> = (0..reps)
        .into_par_iter()
        .map(|_| -> anyhow::Result<Tree> {
            let tmp_sketches = sample_sketches_with_replacement(sketches.clone());
            let tmp_distance = compute_distances(tmp_sketches);
            let tmp_matrix = distance_to_matrix(tmp_distance);
            let tmp_tree_newick = tmp_matrix.compute_newick_tree(tree_algorithm)?;
            Tree::from_newick(&tmp_tree_newick)
                .map_err(|e| anyhow::anyhow!("Tree parsing error: {}", e))
        })
        .collect();

    let trees = trees?;
    // Majority-rule threshold: a clade must appear in more than half the
    // replicate trees to be kept. Using `reps` itself here would require a
    // clade to appear in *every* replicate (strict/unanimous consensus),
    // which collapses to a near-unresolved tree on real data.
    let majority_threshold = reps / 2 + 1;
    let consensus_tree = build_consensus(trees, majority_threshold);
    utils::output_tree(output, consensus_tree.to_newick())
}

pub fn build_jackknife_tree(
    sketches: Vec<Sketch>,
    reps: usize,
    proportion: f64,
    tree_algorithm: &TreeAlgorithm,
    output: Option<String>,
) -> anyhow::Result<()> {
    if !(proportion > 0.0 && proportion < 1.0) {
        anyhow::bail!(
            "jackknife subsampling proportion must be strictly between 0 and 1, got {}",
            proportion
        );
    }

    let trees: anyhow::Result<Vec<Tree>> = (0..reps)
        .into_par_iter()
        .map(|_| -> anyhow::Result<Tree> {
            let tmp_sketches = sample_sketches_without_replacement(sketches.clone(), proportion);
            let tmp_distance = compute_distances(tmp_sketches);
            let tmp_matrix = distance_to_matrix(tmp_distance);
            let tmp_tree_newick = tmp_matrix.compute_newick_tree(tree_algorithm)?;
            Tree::from_newick(&tmp_tree_newick)
                .map_err(|e| anyhow::anyhow!("Tree parsing error: {}", e))
        })
        .collect();

    let trees = trees?;
    // Same majority-rule reasoning as the bootstrap path above.
    let majority_threshold = reps / 2 + 1;
    let consensus_tree = build_consensus(trees, majority_threshold);
    utils::output_tree(output, consensus_tree.to_newick())
}

pub fn build_single_tree(
    sketches: Vec<Sketch>,
    tree_algorithm: &TreeAlgorithm,
    args: &cli::BuildArgs,
) -> anyhow::Result<()> {
    let sketch_distance = compute_distances(sketches);
    let matrix = distance_to_matrix(sketch_distance);
    let newick = matrix.compute_newick_tree(tree_algorithm)?;

    utils::output_tree(args.output.clone(), newick)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    // Helper function to create a simple test tree: ((A,B),(C,D))
    fn create_test_tree() -> Tree {
        let leaf_a = Node::new(Some("A".to_string()), Some(0.1));
        let leaf_b = Node::new(Some("B".to_string()), Some(0.2));
        let leaf_c = Node::new(Some("C".to_string()), Some(0.3));
        let leaf_d = Node::new(Some("D".to_string()), Some(0.4));

        let mut internal1 = Node::new(None, Some(0.5));
        internal1.children = vec![leaf_a, leaf_b];

        let mut internal2 = Node::new(None, Some(0.6));
        internal2.children = vec![leaf_c, leaf_d];

        let mut root = Node::new(None, None);
        root.children = vec![internal1, internal2];

        Tree { root }
    }

    // Helper function to create another test tree: ((A,C),(B,D))
    fn create_alternative_tree() -> Tree {
        let leaf_a = Node::new(Some("A".to_string()), Some(0.15));
        let leaf_c = Node::new(Some("C".to_string()), Some(0.25));
        let leaf_b = Node::new(Some("B".to_string()), Some(0.35));
        let leaf_d = Node::new(Some("D".to_string()), Some(0.45));

        let mut internal1 = Node::new(None, Some(0.55));
        internal1.children = vec![leaf_a, leaf_c];

        let mut internal2 = Node::new(None, Some(0.65));
        internal2.children = vec![leaf_b, leaf_d];

        let mut root = Node::new(None, None);
        root.children = vec![internal1, internal2];

        Tree { root }
    }

    #[test]
    fn benchmark_majority_rule() {
        // Create test data
        let mut clades = Vec::new();
        for i in 0..1000 {
            clades.push(Clade::new(
                vec![format!("A{}", i % 10), format!("B{}", i % 10)],
                Some(i as f64 * 0.001),
            ));
        }

        let start = std::time::Instant::now();
        let result = majority_rule(clades, 50, 1000);
        let duration = start.elapsed();

        println!(
            "Majority rule processed {} clades in {:?}",
            result.len(),
            duration
        );
        assert!(!result.is_empty());
        assert!(result.iter().all(|c| c.support.is_some()));
    }

    #[test]
    fn test_newick_roundtrip() {
        let original = "((A:0.100000,B:0.200000):0.500000,(C:0.300000,D:0.400000):0.600000);";
        let tree = Tree::from_newick(original).expect("Failed to parse");
        let regenerated = tree.to_newick();

        // Parse both to compare structure (floating point precision may differ)
        let tree1 = Tree::from_newick(original).expect("Failed to parse original");
        let tree2 = Tree::from_newick(&regenerated).expect("Failed to parse regenerated");

        // Compare leaf names
        let leaves1 = get_leaf_names(&tree1.root);
        let leaves2 = get_leaf_names(&tree2.root);
        assert_eq!(leaves1, leaves2);
    }

    #[test]
    fn test_majority_rule_filtering() {
        let clades = vec![
            Clade::new(vec!["A".to_string(), "B".to_string()], Some(0.1)),
            Clade::new(vec!["A".to_string(), "B".to_string()], Some(0.2)), // Same clade, different length
            Clade::new(vec!["C".to_string(), "D".to_string()], Some(0.3)),
            Clade::new(vec!["A".to_string(), "C".to_string()], Some(0.4)), // Only appears once
        ];

        // Treat these 4 entries as coming from 4 replicate trees.
        let result = majority_rule(clades, 2, 4);

        // Only {A,B} should pass threshold of 2
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].leaves, vec!["A", "B"]);

        // Branch length is the average of the two observed lengths for
        // {A,B}. Compared with a tolerance since floating-point averaging
        // (0.1 + 0.2) / 2 is not bit-exact to the literal 0.15.
        let length = result[0].length.expect("length should be present");
        assert!(
            (length - 0.15).abs() < 1e-9,
            "unexpected averaged length: {}",
            length
        );
        // Support is the fraction of the 4 input trees containing this
        // clade (2/4) — tracked independently of branch length.
        assert_eq!(result[0].support, Some(0.5));
    }

    #[test]
    fn test_consensus_tree_structure() {
        let tree1 = create_test_tree(); // ((A,B),(C,D))
        let tree2 = create_alternative_tree(); // ((A,C),(B,D))

        let trees = vec![tree1, tree2];
        let consensus = build_consensus(trees, 1); // Keep clades appearing in at least 1 tree

        println!("Consensus tree: {}", consensus.to_newick());

        // No leaf should ever be duplicated under two different branches of
        // the consensus tree (this is what the incompatible-clade bug used
        // to produce: {A,B} and {A,C} both kept as siblings of the root).
        let leaf_names = get_leaf_names(&consensus.root);
        let unique_leaves: HashSet<_> = leaf_names.iter().collect();
        assert_eq!(
            leaf_names.len(),
            unique_leaves.len(),
            "a leaf was duplicated in the consensus tree: {:?}",
            leaf_names
        );
        assert_eq!(unique_leaves.len(), 4);
        assert!(unique_leaves.contains(&&"A".to_string()));
        assert!(unique_leaves.contains(&&"B".to_string()));
        assert!(unique_leaves.contains(&&"C".to_string()));
        assert!(unique_leaves.contains(&&"D".to_string()));

        // The two input trees disagree ((A,B)(C,D) vs (A,C)(B,D)), so with a
        // permissive threshold of 1 the consensus must resolve into exactly
        // one compatible pair of cherries, never a flat/star topology with
        // all four leaves as direct children of the root.
        assert_eq!(consensus.root.children.len(), 2);
        for child in &consensus.root.children {
            assert_eq!(
                child.children.len(),
                2,
                "expected a fully-resolved cherry, got: {:?}",
                child
            );
        }
    }

    #[test]
    fn test_branch_length_averaging_issue() {
        // Branch length and clade support must be tracked independently.
        let clades = vec![
            Clade::new(vec!["A".to_string(), "B".to_string()], Some(0.1)),
            Clade::new(vec!["A".to_string(), "B".to_string()], Some(0.3)),
        ];

        let result = majority_rule(clades, 2, 2);
        let kept_clade = result.first().expect("clade should pass threshold");

        // Length is the average of the two observed values (tolerance since
        // floating-point averaging isn't guaranteed bit-exact)...
        let length = kept_clade.length.expect("length should be present");
        assert!(
            (length - 0.2).abs() < 1e-9,
            "unexpected averaged length: {}",
            length
        );
        // ...while support reflects that the clade occurred in both of the
        // 2 input trees, and is not conflated with the length above.
        assert_eq!(kept_clade.support, Some(1.0));
    }

    #[test]
    fn test_jackknife_proportion_validation() {
        // proportion must be strictly between 0 and 1; out-of-range values
        // should fail fast rather than silently keep 0 or all hashes.
        assert!(build_jackknife_tree(vec![], 5, 1.5, &TreeAlgorithm::Canonical, None).is_err());
        assert!(build_jackknife_tree(vec![], 5, 0.0, &TreeAlgorithm::Canonical, None).is_err());
        assert!(build_jackknife_tree(vec![], 5, 1.0, &TreeAlgorithm::Canonical, None).is_err());
        assert!(build_jackknife_tree(vec![], 5, -0.1, &TreeAlgorithm::Canonical, None).is_err());
    }

    fn make_test_sketch(name: &str, n_hashes: u64) -> Sketch {
        Sketch {
            name: name.to_string(),
            seq_length: n_hashes * 21,
            num_valid_kmers: n_hashes,
            comment: String::new(),
            hashes: (0..n_hashes)
                .map(|h| finch::sketch_schemes::KmerCount {
                    hash: h,
                    kmer: Vec::new(),
                    count: 1,
                    extra_count: 0,
                    label: None,
                })
                .collect(),
            filter_params: finch::filtering::FilterParams::default(),
            sketch_params: finch::sketch_schemes::SketchParams::default(),
        }
    }

    #[test]
    fn test_sample_sketches_without_replacement_keeps_proportion() {
        let sketches = vec![make_test_sketch("genome1", 1000)];
        let sampled = sample_sketches_without_replacement(sketches, 0.5);

        assert_eq!(sampled.len(), 1);
        assert_eq!(sampled[0].hashes.len(), 500);
        assert_eq!(sampled[0].name, "genome1");

        // sampling is without replacement: no duplicate hash values
        let unique: HashSet<_> = sampled[0].hashes.iter().map(|k| k.hash).collect();
        assert_eq!(unique.len(), 500);
    }

    #[test]
    fn test_sample_sketches_without_replacement_never_keeps_zero() {
        // A tiny sketch with a very small proportion should still keep at
        // least 1 hash rather than rounding down to an empty sketch.
        let sketches = vec![make_test_sketch("genome1", 3)];
        let sampled = sample_sketches_without_replacement(sketches, 0.1);

        assert_eq!(sampled[0].hashes.len(), 1);
    }

    #[test]
    fn test_newick_parsing_edge_cases() {
        // Test names without branch lengths
        let result = Tree::from_newick("(A,B);");
        println!("{:?}", result);
        assert!(result.is_ok());

        // Test malformed input
        let result = Tree::from_newick("(A,B;"); // Missing closing parenthesis
        println!("{:?}", result);
        assert!(result.is_err());

        // Test single leaf
        let result = Tree::from_newick("A;");
        println!("{:?}", result);
        assert!(result.is_ok());
    }

    #[test]
    fn test_building_trees() {
        let inputs = vec![
            "((A:0.1,B:0.2):0.3,C:0.6);",
            "((A:0.15,B:0.25):0.35,C:0.55);",
            "((A:0.12,B:0.22):0.33,C:0.58);",
        ];
        let mut trees = Vec::new();
        for input_str in inputs {
            match Tree::from_newick(input_str) {
                Ok(tree) => trees.push(tree),
                Err(e) => eprintln!("Error: {}", e),
            }
        }
        let consensus = build_consensus(trees, 3);
        println!("{}", consensus.to_newick());
    }

    #[test]
    fn test_sample_sketches_with_replacement_preserves_own_size() {
        let sketches = vec![make_test_sketch("genome1", 1000)];
        let sampled = sample_sketches_with_replacement(sketches);

        assert_eq!(sampled.len(), 1);
        assert_eq!(sampled[0].hashes.len(), 1000);
        assert_eq!(sampled[0].name, "genome1");
    }

    #[test]
    fn test_sample_sketches_with_replacement_respects_each_sketchs_own_size() {
        // Sketches of different lengths must each be resampled to their own
        // length, not forced to one shared/global size.
        let sketches = vec![make_test_sketch("small", 10), make_test_sketch("big", 500)];
        let sampled = sample_sketches_with_replacement(sketches);

        let small = sampled.iter().find(|s| s.name == "small").unwrap();
        let big = sampled.iter().find(|s| s.name == "big").unwrap();
        assert_eq!(small.hashes.len(), 10);
        assert_eq!(big.hashes.len(), 500);
    }
}
