// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use finch::serialization::Sketch;
use std::collections::{HashMap, HashSet};

use rand::{rng, seq::IndexedRandom};

#[derive(Debug, Clone, PartialEq)]
pub struct Clade {
    pub leaves: Vec<String>, // sorted leaves
    pub length: Option<f64>, // not used for Hash/Eq
}

#[derive(Debug, Clone)]
pub struct Node {
    pub name: Option<String>,
    pub length: Option<f64>,
    pub children: Vec<Node>,
}

impl Node {
    pub fn new(name: Option<String>, length: Option<f64>) -> Self {
        Node {
            name,
            length,
            children: Vec::new(),
        }
    }
}

#[derive(Debug)]
pub struct Tree {
    pub root: Node,
}

impl Tree {
    pub fn from_newick(s: &str) -> Result<Self, String> {
        let chars: Vec<char> = s.chars().collect();
        let mut pos = 0;
        let root = parse_subtree(&chars, &mut pos)?;
        Ok(Self { root })
    }
}

fn parse_subtree(chars: &[char], pos: &mut usize) -> Result<Node, String> {
    skip_whitespace(chars, pos);

    if *pos >= chars.len() {
        return Err("Unexpected end of string".into());
    }

    let mut node = Node::new(None, None);

    if chars[*pos] == '(' {
        // internal node
        *pos += 1; // skip '('
        loop {
            let child = parse_subtree(chars, pos)?;
            node.children.push(child);

            skip_whitespace(chars, pos);
            if chars.get(*pos) == Some(&',') {
                *pos += 1; // skip comma
            } else {
                break;
            }
        }
        if chars.get(*pos) != Some(&')') {
            return Err("Expected ')'".into());
        }
        *pos += 1; // skip ')'
    } else {
        // leaf node
        let name = parse_name(chars, pos);
        node.name = Some(name);
    }

    // optional branch length
    skip_whitespace(chars, pos);
    if chars.get(*pos) == Some(&':') {
        *pos += 1; // skip ':'
        let len_str = parse_branch_length(chars, pos);
        node.length = len_str.parse::<f64>().ok();
    }

    skip_whitespace(chars, pos);

    Ok(node)
}

fn parse_name(chars: &[char], pos: &mut usize) -> String {
    let mut name = String::new();
    while *pos < chars.len() {
        match chars[*pos] {
            ':' | ',' | ')' | '(' | ';' => break,
            c => {
                name.push(c);
                *pos += 1;
            }
        }
    }
    name.trim().to_string()
}

fn parse_branch_length(chars: &[char], pos: &mut usize) -> String {
    let mut num = String::new();
    while *pos < chars.len() {
        match chars[*pos] {
            ',' | ')' | ';' => break,
            c => {
                num.push(c);
                *pos += 1;
            }
        }
    }
    num.trim().to_string()
}

fn skip_whitespace(chars: &[char], pos: &mut usize) {
    while *pos < chars.len() && chars[*pos].is_whitespace() {
        *pos += 1;
    }
}

/// Recursively build a tree from clades
fn build_tree_from_clades(all_leaves: &[String], clades: &[Clade]) -> Tree {
    let root_node = build_node(all_leaves, clades);
    Tree { root: root_node }
}

fn build_node(current_leaves: &[String], clades: &[Clade]) -> Node {
    // Find child clades that are proper subsets of current leaves
    let child_clades: Vec<&Clade> = clades
        .iter()
        .filter(|c| {
            c.leaves.len() < current_leaves.len() && is_subset_sorted(&c.leaves, current_leaves)
        })
        .collect();

    if child_clades.is_empty() {
        // Leaf node
        if current_leaves.len() != 1 {
            panic!("Expected single leaf, got {:?}", current_leaves);
        }
        return Node::new(Some(current_leaves[0].clone()), None);
    }

    // Keep only maximal subsets to avoid duplicates
    let maximal_clades = maximal_subsets(&child_clades);

    let mut children_nodes = Vec::new();
    for clade in maximal_clades {
        let child_node = build_node(&clade.leaves, clades);
        let mut node_with_len = child_node.clone();
        node_with_len.length = clade.length;
        children_nodes.push(node_with_len);
    }

    Node {
        name: None,
        length: None,
        children: children_nodes,
    }
}

/// Returns only clades that are not contained in any other clade in the set
fn maximal_subsets<'a>(clades: &[&'a Clade]) -> Vec<&'a Clade> {
    clades
        .iter()
        .filter(|c1| {
            !clades
                .iter()
                .any(|c2| *c1 != c2 && is_subset_sorted(&c1.leaves, &c2.leaves))
        })
        .copied()
        .collect()
}

/// Check if a sorted slice `a` is a subset of sorted slice `b`
fn is_subset_sorted(a: &[String], b: &[String]) -> bool {
    let mut i = 0;
    let mut j = 0;

    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Equal => {
                i += 1;
                j += 1;
            }
            std::cmp::Ordering::Greater => {
                return false;
            }
            std::cmp::Ordering::Less => {
                j += 1;
            }
        }
    }
    i == a.len()
}

/// Majority rule: keep clades occurring >= threshold times
fn majority_rule(clades: &[Clade], threshold: usize) -> Vec<Clade> {
    let mut counts: HashMap<Vec<String>, Vec<&Clade>> = HashMap::new();

    for clade in clades {
        let mut leaves_sorted = clade.leaves.clone();
        leaves_sorted.sort(); // ensure deterministic order
        counts.entry(leaves_sorted).or_default().push(clade);
    }

    let mut kept = Vec::new();
    for (_leaves, clade_refs) in counts.into_iter() {
        if clade_refs.len() >= threshold {
            // Average branch length over occurrences
            let avg_length = if clade_refs.iter().any(|c| c.length.is_some()) {
                Some(
                    clade_refs.iter().filter_map(|c| c.length).sum::<f64>()
                        / clade_refs.len() as f64,
                )
            } else {
                None
            };

            let mut clade_copy = (*clade_refs[0]).clone();
            clade_copy.length = avg_length;
            kept.push(clade_copy);
        }
    }

    kept
}

/// Extract clades from a tree, returning Vec<Clade> with sorted leaves
fn get_clades(node: &Node) -> Vec<Clade> {
    let mut clades = Vec::new();

    fn helper(node: &Node, clades: &mut Vec<Clade>) -> Vec<String> {
        if node.children.is_empty() {
            // Tip node: return its leaf name
            return vec![node.name.clone().unwrap()];
        }

        let mut taxa = Vec::new();
        for child in &node.children {
            let child_taxa = helper(child, clades);
            taxa.extend(child_taxa);
        }

        taxa.sort(); // ensure deterministic order
        clades.push(Clade {
            leaves: taxa.clone(),
            length: node.length,
        });

        taxa
    }

    helper(node, &mut clades);
    clades
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

fn build_consensus(trees: Vec<Tree>, threshold: usize) -> Tree {
    let mut clades: Vec<Clade> = Vec::new();

    for tree in trees {
        clades.extend(get_clades(&tree.root));
    }

    let all_leaves: Vec<String> = clades
        .iter()
        .flat_map(|c| c.leaves.clone())
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();
    let kept_clades = majority_rule(&clades, threshold);
    let consensus = build_tree_from_clades(&all_leaves, &kept_clades);

    consensus
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_clades() {
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
        let consensus = build_consensus(trees, 2);
        println!("{:?}", consensus);
    }
}
