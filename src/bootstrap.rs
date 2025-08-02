// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::collections::{HashMap, HashSet};

use finch::serialization::Sketch;
use phylo::{
    prelude::{Newick, RootedMetaTree, RootedTree, RootedTreeNode},
    tree::PhyloTree,
};
use rand::{rng, seq::IndexedRandom};

type Clade = HashSet<String>;

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

fn extract_clades(newick: String) -> anyhow::Result<Vec<Clade>> {
    let mut clades = Vec::new();

    let tree = PhyloTree::from_newick(newick.as_bytes())?;

    for node in tree.get_node_ids() {
        if tree.is_leaf(node) {
            continue;
        }

        let mut leaf_names = HashSet::new();
        for leaf in tree.get_node_children(node) {
            if let Some(label) = tree.get_node_taxa(leaf.get_id()) {
                leaf_names.insert(label.clone());
            }
        }

        // To keep only non-trivial clades (size < total leaf set)
        // Internal node are does not have children so it result in an
        // empty set which i discard
        if !leaf_names.is_empty() {
            clades.push(leaf_names);
        }
    }

    for clade in &clades {
        println!("{:?}", clade);
    }

    Ok(clades)
}

/// Takes a list of clades and returns the majority-rule clades.
/// A clade is a HashSet<String> reprensenting a group of taxa.
/// Only clades appearint in more than `threshold` trees are kept.
fn majority_rule_clades(
    clades: Vec<HashSet<String>>,
    total_trees: usize,
    threshold: f64,
) -> anyhow::Result<Vec<HashSet<String>>> {
    let mut clade_counts: HashMap<Vec<String>, usize> = HashMap::new();

    for clade in clades {
        // convert hahset to sorted vec<string> for hashing
        let mut sorted_clade: Vec<String> = clade.iter().cloned().collect();
        sorted_clade.sort();
        *clade_counts.entry(sorted_clade).or_insert(0) += 1;
    }

    let required_count = (threshold * total_trees as f64).ceil() as usize;
    let majority_clades: Vec<HashSet<String>> = clade_counts
        .into_iter()
        .filter_map(|(sorted, count)| {
            if count >= required_count {
                Some(sorted.into_iter().collect::<HashSet<String>>())
            } else {
                None
            }
        })
        .collect();

    Ok(majority_clades)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_clades() {
        assert!(extract_clades("((A,B),(C,D));".to_string()).is_ok());
        //assert!(extract_clades("((A,C),(B,D));".to_string()).is_ok());
    }
}
