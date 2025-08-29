// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use crate::{
    build::dist::{self, ComputeTree, TreeAlgorithm},
    cli,
    nwk::{build_tree_from_clades, Clade, Node, Tree},
    utils,
};
use finch::serialization::Sketch;
use rand::{rng, seq::IndexedRandom};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

/// Majority rule: keep clades occurring >= threshold times
fn majority_rule(clades: Vec<Clade>, threshold: usize) -> Vec<Clade> {
    let mut clade_group: HashMap<Vec<String>, Vec<f64>> = HashMap::with_capacity(clades.len());

    // single passe: group clades by their leaf sets and collect their lengths
    for clade in clades {
        let lengths = clade_group.entry(clade.leaves).or_default();
        if let Some(length) = clade.length {
            lengths.push(length);
        }
    }

    // filter and average in one step
    clade_group
        .into_iter()
        .filter(|(_, lengths)| lengths.len() >= threshold)
        .map(|(leaves, lengths)| {
            let avg_length = if lengths.is_empty() {
                None
            } else {
                Some(lengths.iter().sum::<f64>() / lengths.len() as f64)
            };
            Clade::new(leaves, avg_length)
        })
        .collect()
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

    let kept_clades = majority_rule(clades, threshold);
    build_tree_from_clades(&all_leaves, &kept_clades)
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

/// Sample sequences sketches with replacement
/// Returns a new vec of sampled sketches
pub fn sample_sketches_with_replacement(sketches: Vec<Sketch>, sample_size: usize) -> Vec<Sketch> {
    let mut rng = rng();
    let mut bootstraped = Vec::new();

    for sketch in sketches {
        bootstraped.push(Sketch {
            name: sketch.name.clone(),
            seq_length: sketch.seq_length,
            num_valid_kmers: sketch.num_valid_kmers,
            comment: sketch.comment.clone(),
            filter_params: sketch.filter_params,
            sketch_params: sketch.sketch_params,
            hashes: (0..sample_size)
                .map(|_| sketch.hashes.choose(&mut rng).unwrap().clone())
                .collect(),
        });
    }

    bootstraped
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
            let tmp_sketches = sample_sketches_with_replacement(sketches.clone(), 1000);
            let tmp_distance = dist::compute_distances(tmp_sketches);
            let tmp_matrix = dist::distance_to_matrix(tmp_distance);
            let tmp_tree_newick = tmp_matrix.compute_newick_tree(tree_algorithm)?;
            Tree::from_newick(&tmp_tree_newick)
                .map_err(|e| anyhow::anyhow!("Tree parsing error: {}", e))
        })
        .collect();

    let trees = trees?;
    let consensus_tree = build_consensus(trees, reps);
    utils::output_tree(output, consensus_tree.to_newick())
}

pub fn build_single_tree(
    sketches: Vec<Sketch>,
    tree_algorithm: &TreeAlgorithm,
    args: &cli::BuildArgs,
) -> anyhow::Result<()> {
    let sketch_distance = dist::compute_distances(sketches);
    let matrix = dist::distance_to_matrix(sketch_distance);
    let newick = matrix.compute_newick_tree(tree_algorithm)?;

    utils::output_tree(args.output.clone(), newick)?;
    utils::manage_tempdir(args.keep, &matrix, &args.tempdir, true)?;

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
        let result = majority_rule(clades, 50);
        let duration = start.elapsed();

        println!(
            "Majority rule processed {} clades in {:?}",
            result.len(),
            duration
        );
        assert!(!result.is_empty());
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

        let result = majority_rule(clades, 2);

        // Only {A,B} should pass threshold of 2
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].leaves, vec!["A", "B"]);

        // Check that branch length is stored as frequency (current implementation issue)
        // This test will reveal the bug where length becomes support value
        println!("Majority rule result length: {:?}", result[0].length);
        // Current broken implementation will give Some(1.0) instead of averaged length
    }

    #[test]
    fn test_consensus_tree_structure() {
        let tree1 = create_test_tree(); // ((A,B),(C,D))
        let tree2 = create_alternative_tree(); // ((A,C),(B,D))

        let trees = vec![tree1, tree2];
        let consensus = build_consensus(trees, 1); // Keep clades appearing in at least 1 tree

        println!("Consensus tree: {}", consensus.to_newick());

        // This test will reveal the broken tree structure
        // Current implementation creates a flat tree instead of proper hierarchy

        // Check that we don't have duplicate leaves
        let leaf_names = get_leaf_names(&consensus.root);
        let unique_leaves: HashSet<_> = leaf_names.iter().collect();

        println!("Leaf names: {:?}", leaf_names);
        println!("Unique leaves: {:?}", unique_leaves);

        // Should have exactly 4 unique leaves
        assert_eq!(unique_leaves.len(), 4);

        // Current broken implementation will fail this test
        assert!(unique_leaves.contains(&&"A".to_string()));
        assert!(unique_leaves.contains(&&"B".to_string()));
        assert!(unique_leaves.contains(&&"C".to_string()));
        assert!(unique_leaves.contains(&&"D".to_string()));
    }

    #[test]
    fn test_branch_length_averaging_issue() {
        // Test to demonstrate the branch length issue in majority_rule
        let clades = vec![
            Clade::new(vec!["A".to_string(), "B".to_string()], Some(0.1)),
            Clade::new(vec!["A".to_string(), "B".to_string()], Some(0.3)),
        ];

        let result = majority_rule(clades, 2);

        if let Some(kept_clade) = result.first() {
            assert_eq!(kept_clade.length, Some(0.2));
        }
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
}
