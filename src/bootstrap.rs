// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::collections::{HashMap, HashSet};

use finch::serialization::Sketch;
use phylo::{
    node::Node,
    prelude::{
        Newick, RootedMetaNode, RootedMetaTree, RootedTree, RootedTreeNode, RootedWeightedNode,
        RootedWeightedTree, TreeNodeID, DFS,
    },
    tree::{PhyloTree, SimpleRootedTree},
};
use rand::{rng, seq::IndexedRandom};

type Clade = HashSet<String>;
type Taxon = String;
type CladeKey = Vec<Taxon>;

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

fn build_parent_map(tree: &PhyloTree) -> HashMap<usize, Option<usize>> {
    let mut parent_map: HashMap<usize, Option<usize>> = HashMap::new();
    for node_id in tree.get_nodes().map(|n| n.get_id()) {
        parent_map.insert(node_id, None);
    }

    for parent in tree.get_nodes() {
        for child_id in parent.get_children() {
            parent_map.insert(child_id, Some(parent.get_id()));
        }
    }
    parent_map
}

/// For a single tree, return:
/// - the full taxon set,
/// - a map from internal clade (sorted taxa vec, size >= 2, <all) -> edge length to parent,
/// - a map from leaf taxon -> terminal edge length.
fn collect_clade_and_leaf_lengths(tree: &PhyloTree) -> Vec<(Vec<String>, f64)> {
    let parent_map = build_parent_map(tree);
    let mut results: Vec<(Vec<String>, f64)> = Vec::new();

    for node in tree.get_nodes() {
        let branch_length: f64 = parent_map[&node.get_id()]
            .and_then(|p_id| tree.get_edge_weight(p_id, node.get_id()))
            .map(|w| w as f64)
            .unwrap_or(0.0);

        if let Some(child_id) = node.get_children().next() {
            let taxa = collect_descendant_taxa(tree, child_id);
            results.push((taxa, branch_length));
        } else {
            results.push((vec![node.get_taxa().unwrap().to_string()], branch_length));
        }
    }

    results
}

fn collect_descendant_taxa(tree: &PhyloTree, node_id: usize) -> Vec<String> {
    let mut taxa = Vec::new();
    let mut stack = vec![node_id];

    while let Some(current) = stack.pop() {
        let current_node = tree.get_node(current).unwrap();
        if let Some(_) = current_node.get_children().next() {
            stack.extend(current_node.get_children());
        } else {
            taxa.push(current_node.get_taxa().unwrap().to_string());
        }
    }

    taxa.sort();
    taxa
}

/* ----------- topology construction ----------- */
fn leaves_of(tree: &PhyloTree, node: &Node<String, f32, f32>, acc: &mut HashSet<Taxon>) {
    if let Some(n) = &node.get_taxa() {
        acc.insert(n.to_string());
    } else {
        for c in node.get_children() {
            leaves_of(tree, tree.get_node(c).unwrap(), acc);
        }
    }
}

fn node_leaf_set(tree: &PhyloTree, node: &Node<String, f32, f32>) -> HashSet<Taxon> {
    let mut s = HashSet::new();
    leaves_of(tree, node, &mut s);
    s
}

fn clade_subset(clade: &CladeKey, set: &HashSet<Taxon>) -> bool {
    clade.iter().all(|t| set.contains(t))
}

fn find_minimal_container_id(tree: &PhyloTree, node_id: usize, clade: &CladeKey) -> Option<usize> {
    let node_ref = tree.get_node(node_id).unwrap();

    if node_ref.get_taxa().is_some() {
        return None;
    }

    let here = node_leaf_set(tree, node_ref);
    if !clade_subset(clade, &here) {
        return None;
    }

    for child_id in node_ref.get_children() {
        let child_leaves = {
            let child_ref = tree.get_node(child_id).unwrap();
            node_leaf_set(tree, child_ref)
        };
        if clade_subset(clade, &child_leaves) {
            if let Some(found) = find_minimal_container_id(tree, child_id, clade) {
                return Some(found);
            }
            return Some(child_id);
        }
    }

    Some(node_id)
}

fn find_minimal_container_mut<'a>(
    tree: &'a mut PhyloTree,
    start_node_id: usize,
    clade: &CladeKey,
) -> Option<&'a mut Node<String, f32, f32>> {
    if let Some(id) = find_minimal_container_id(tree, start_node_id, clade) {
        // Only *now* do we take the mutable borrow
        Some(tree.get_node_mut(id).unwrap())
    } else {
        None
    }
}

/*fn insert_clade(root: &mut Node, clade: &CladeKey, clade_len: f64) {

}*/

/*
/// Takes a list of clades and returns the majority-rule clades.
/// A clade is a HashSet<String> representing a group of taxa.
/// Only clades appearing in more than `threshold` trees are kept.
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
}*/

/*
fn build_consensus_tree(clades: Vec<Clade>) -> anyhow::Result<SimpleRootedTree<String, f32, f32>> {
    let all_taxa: HashSet<String> = clades.iter().flat_map(|c| c.clone()).collect();

    let mut tree = SimpleRootedTree::new(0);

    fn insert_clade(node: &mut Node<String, f32, f32>, clade: &Clade) -> bool {
        let node_taxa: HashSet<String> = get_leaves(node);

        if clade.is_subset(&node_taxa) {
            for child in &mut node.get_children() {
            let mut child_node = Node::new(child);
            child_node.set_taxa()
                if insert_clade(node.)
            }
        }
        false
    }

    Ok(tree)
}*/

fn get_leaves(node: Node<String, f32, f32>) -> HashSet<String> {
    let mut leaves = HashSet::new();
    if let Some(name) = &node.get_taxa() {
        leaves.insert(name.to_string());
    }

    for child in node.get_children() {
        leaves.extend(get_leaves(Node::new(child)));
    }
    leaves
}

#[cfg(test)]
mod tests {
    use super::*;

    /*#[test]
    fn test_extract_clades() {
        let input_str = String::from("((A:0.1,B:0.2),C:0.6);");
        let tree = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
        println!("{:?}", tree);
        assert!(extract_clades("((A,B),(C,D));".to_string()).is_ok());
        //assert!(extract_clades("((A,C),(B,D));".to_string()).is_ok());
    }*/
}
