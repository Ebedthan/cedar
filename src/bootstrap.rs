// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use crate::newick::{Clade, Node, Tree};
use finch::serialization::Sketch;
use rand::{rng, seq::IndexedRandom};
use std::collections::{HashMap, HashSet};

/// Recursively build a tree from clades
fn build_tree_from_clades(all_leaves: &[String], clades: &[Clade]) -> Tree {
    // Root = clade containing all leaves
    let mut sorted_all = all_leaves.to_vec();
    sorted_all.sort();
    let root_clade = clades
        .iter()
        .find(|c| c.leaves == sorted_all)
        .cloned()
        .unwrap_or_else(|| Clade::new(sorted_all.clone(), None));

    let root = build_node(&root_clade, clades);
    Tree { root }
}

fn build_node(clade: &Clade, all_clades: &[Clade]) -> Node {
    if clade.leaves.len() == 1 {
        // Leaf
        return Node {
            name: Some(clade.leaves[0].clone()),
            length: clade.length,
            children: vec![],
        };
    }

    let clade_set: HashSet<String> = clade.leaf_set();
    let mut children = Vec::new();
    let mut covered = HashSet::new();

    // Find maximal proper subclades of this clade
    for sub in all_clades {
        if sub.leaves.len() < clade.leaves.len() {
            let sub_set = sub.leaf_set();
            if sub_set.is_subset(&clade_set) {
                // ensure it's maximal
                let contained_in_other = all_clades.iter().any(|other| {
                    other.leaves.len() < clade.leaves.len()
                        && other.leaves.len() > sub.leaves.len()
                        && sub_set.is_subset(&other.leaf_set())
                        && other.leaf_set().is_subset(&clade_set)
                });
                if !contained_in_other {
                    children.push(build_node(sub, all_clades));
                    covered.extend(sub_set);
                }
            }
        }
    }

    // Add missing leaves not covered by any subclade
    for taxon in clade_set.difference(&covered) {
        children.push(Node {
            name: Some(taxon.clone()),
            length: all_clades
                .iter()
                .find(|c| c.leaves.len() == 1 && c.leaves[0] == *taxon)
                .and_then(|c| c.length),
            children: vec![],
        });
    }

    children.sort_by(|a, b| {
        let aname = a.name.clone().unwrap_or_default();
        let bname = b.name.clone().unwrap_or_default();
        aname.cmp(&bname)
    });

    Node {
        name: None,
        length: clade.length,
        children,
    }
}

/// Majority rule: keep clades occurring >= threshold times
fn majority_rule(clades: Vec<Clade>, threshold: usize) -> Vec<Clade> {
    let mut clade_map: HashMap<Vec<String>, (usize, Vec<f64>)> = HashMap::new();

    // group clades by their leaf sets and collect their lengths
    for clade in clades {
        let entry = clade_map.entry(clade.leaves).or_insert((0, Vec::new()));
        entry.0 += 1;
        if let Some(length) = clade.length {
            entry.1.push(length);
        }
    }

    let mut kept = Vec::new();
    for (leaves, (count, lengths)) in clade_map {
        if count >= threshold {
            // average the branch lengths
            let avg_length = if !lengths.is_empty() {
                Some(lengths.iter().sum::<f64>() / lengths.len() as f64)
            } else {
                None
            };
            // always keep singleton clades with length
            if leaves.len() == 1 && avg_length.is_none() {
                continue;
            }

            kept.push(Clade::new(leaves, avg_length));
        }
    }
    kept
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

/*
#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

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
    fn test_clade_extraction() {
        let tree = create_test_tree();
        let clades = get_clades(&tree.root);

        // Should have 3 clades: {A,B}, {C,D}, {A,B,C,D}
        assert_eq!(clades.len(), 3);

        let clade_sets: HashSet<Vec<String>> = clades.iter().map(|c| c.leaves.clone()).collect();

        assert!(clade_sets.contains(&vec!["A".to_string(), "B".to_string()]));
        assert!(clade_sets.contains(&vec!["C".to_string(), "D".to_string()]));
        assert!(clade_sets.contains(&vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string()
        ]));
    }

    #[test]
    fn test_majority_rule_filtering() {
        let clades = vec![
            Clade::new(vec!["A".to_string(), "B".to_string()], Some(0.1)),
            Clade::new(vec!["A".to_string(), "B".to_string()], Some(0.2)), // Same clade, different length
            Clade::new(vec!["C".to_string(), "D".to_string()], Some(0.3)),
            Clade::new(vec!["A".to_string(), "C".to_string()], Some(0.4)), // Only appears once
        ];

        let result = majority_rule(clades, 2, 2);

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
    fn test_is_subset_sorted() {
        let a = vec!["A".to_string(), "B".to_string()];
        let b = vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string(),
        ];
        let c = vec!["B".to_string(), "E".to_string()];

        assert!(is_subset_sorted(&a, &b));
        assert!(!is_subset_sorted(&b, &a));
        assert!(!is_subset_sorted(&c, &b));
        assert!(is_subset_sorted(&a, &a));
    }

    #[test]
    fn test_maximal_subsets() {
        let clade1 = Clade::new(vec!["A".to_string(), "B".to_string()], Some(0.1));
        let clade2 = Clade::new(
            vec!["A".to_string(), "B".to_string(), "C".to_string()],
            Some(0.2),
        );
        let clade3 = Clade::new(vec!["C".to_string(), "D".to_string()], Some(0.3));

        let clades = vec![&clade1, &clade2, &clade3];
        let maximal = maximal_subsets(&clades);

        // clade1 should be filtered out as it's a subset of clade2
        assert_eq!(maximal.len(), 2);
        assert!(maximal.contains(&&clade2));
        assert!(maximal.contains(&&clade3));
        assert!(!maximal.contains(&&clade1));
    }

    #[test]
    fn test_branch_length_averaging_issue() {
        // Test to demonstrate the branch length issue in majority_rule
        let clades = vec![
            Clade::new(vec!["A".to_string(), "B".to_string()], Some(0.1)),
            Clade::new(vec!["A".to_string(), "B".to_string()], Some(0.3)),
        ];

        let result = majority_rule(clades, 2, 2);

        if let Some(kept_clade) = result.first() {
            // Current implementation gives Some(1.0) instead of Some(0.2)
            // This test will fail with current implementation
            println!("Expected average: 0.2, Got: {:?}", kept_clade.length);

            // This assertion will fail, revealing the bug
            assert_eq!(kept_clade.length, Some(0.2));

            // Current broken implementation produces this:
            // assert_eq!(kept_clade.length, Some(1.0));
        }
    }

    #[test]
    fn test_tree_building_creates_invalid_structure() {
        // Test to demonstrate the tree building issue
        let clades = vec![
            Clade::new(vec!["A".to_string(), "B".to_string()], Some(0.1)),
            Clade::new(vec!["C".to_string(), "D".to_string()], Some(0.2)),
        ];

        let all_leaves = vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string(),
        ];
        let tree = build_tree_from_clades(&all_leaves, &clades);

        println!("Built tree structure:");
        print_tree_structure(&tree.root, 0);

        // Current implementation creates invalid structure:
        // - Each clade becomes a direct child of root
        // - Leaves are duplicated in multiple places
        // - Branch lengths become node names

        let leaf_count = count_total_leaves(&tree.root);
        println!("Total leaf count: {}", leaf_count);

        // This will show the problem - more leaves than expected due to duplication
        // Should be 4, but current implementation creates more
        assert!(leaf_count > 4); // This assertion demonstrates the bug
    }

    fn count_total_leaves(node: &Node) -> usize {
        if node.children.is_empty() {
            return 1;
        }

        node.children
            .iter()
            .map(|child| count_total_leaves(child))
            .sum()
    }

    fn print_tree_structure(node: &Node, depth: usize) {
        let indent = "  ".repeat(depth);
        if node.children.is_empty() {
            println!("{}Leaf: {:?} (len: {:?})", indent, node.name, node.length);
        } else {
            println!(
                "{}Internal: {:?} (len: {:?}) - {} children",
                indent,
                node.name,
                node.length,
                node.children.len()
            );
            for child in &node.children {
                print_tree_structure(child, depth + 1);
            }
        }
    }

    #[test]
    fn test_newick_parsing_edge_cases() {
        // Test empty names
        let result = Tree::from_newick("(,);");
        assert!(result.is_ok());

        // Test names without branch lengths
        let result = Tree::from_newick("(A,B);");
        assert!(result.is_ok());

        // Test malformed input
        let result = Tree::from_newick("(A,B;"); // Missing closing parenthesis
        assert!(result.is_err());

        // Test single leaf
        let result = Tree::from_newick("A;");
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
}*/
