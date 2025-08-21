// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

/// Newick Format Minimal Tree Parser Module
///
/// This module contains function to reads a tree from a newick file
/// and store the result as a tree struct. It also provide convenient
/// functions to output newick format from a tree struct and simple functions
/// to read node of a tree (into `Node` struct) and clade of a tree (into `Clade` struct).
///
///
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, PartialEq)]
pub struct Clade {
    pub leaves: Vec<String>, // sorted leaves for deterministic behavior
    pub length: Option<f64>, // branch length, excluded from equality
}

impl Clade {
    pub fn new(mut leaves: Vec<String>, length: Option<f64>) -> Self {
        leaves.sort(); // enforce deterministic ordering in O(n log n)
        leaves.dedup(); // remove duplicates in O(n)
        Clade { leaves, length }
    }

    /// Return the set of leaf of a clade
    /// Use cached hash set for frequent subset operations
    /// Time: O(n) vs O(n) for each conversion, but avoids repeated allocations
    pub fn leaf_set(&self) -> HashSet<String> {
        self.leaves.iter().cloned().collect()
    }

    /// Subset checking using sorted vectors
    /// More efficient than hashset operations with O(min(a,b)) time instead of O(a * b)
    pub fn is_subset_of_sorted(&self, other: &[String]) -> bool {
        if self.leaves.len() > other.len() {
            return false;
        }

        let mut i = 0;
        let mut j = 0;

        while i < self.leaves.len() && j < other.len() {
            match self.leaves[i].cmp(&other[j]) {
                std::cmp::Ordering::Equal => {
                    i += 1;
                    j += 1;
                }
                // self[i] not in other
                std::cmp::Ordering::Less => return false,
                std::cmp::Ordering::Greater => j += 1, // advance other
            }
        }

        i == self.leaves.len() // all elements of self were found
    }
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

    fn to_newick(&self) -> String {
        if self.children.is_empty() {
            // leaf node
            let name = self.name.clone().unwrap_or_default();
            match self.length {
                Some(len) => format!("{}:{:.4}", name, len),
                None => name,
            }
        } else {
            // internal node
            let children_str: Vec<String> = self.children.iter().map(|c| c.to_newick()).collect();
            let joined = children_str.join(",");

            let label = self.name.clone().unwrap_or_default();
            match self.length {
                Some(len) if !label.is_empty() => format!("({}){}:{:.4}", joined, label, len),
                Some(len) => format!("({}):{:.4}", joined, len),
                None if !label.is_empty() => format!("({}){}", joined, label),
                None => format!("({})", joined),
            }
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

    pub fn to_newick(&self) -> String {
        format!("{};", self.root.to_newick())
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
                   // internal node label
        skip_whitespace(chars, pos);
        if *pos < chars.len() {
            if let Some(c) = chars.get(*pos) {
                if !matches!(c, ':' | ',' | ')' | ';') {
                    let label = parse_name(chars, pos);
                    if !label.is_empty() {
                        node.name = Some(label);
                    }
                }
            }
        }
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
pub fn build_tree_from_clades(all_leaves: &[String], clades: &[Clade]) -> Tree {
    // Create lookup structures once
    let mut clade_by_leaves: HashMap<Vec<String>, &Clade> = HashMap::new();
    let mut clades_by_size: HashMap<usize, Vec<&Clade>> = HashMap::new();

    // Pre-process clades for faster lookup
    for clade in clades {
        clade_by_leaves.insert(clade.leaves.clone(), clade);
        clades_by_size
            .entry(clade.leaves.len())
            .or_default()
            .push(clade);
    }

    // sort all leaves
    let mut sorted_all = all_leaves.to_vec();
    sorted_all.sort();

    // own a fallback root clade if not present in the input
    let fallback_root;
    let root_clade: &Clade = if let Some(c) = clade_by_leaves.get(&sorted_all) {
        c
    } else {
        fallback_root = Clade::new(sorted_all.clone(), None);
        &fallback_root
    };

    let root = build_node(root_clade, &clade_by_leaves, &clades_by_size);
    Tree { root }
}

fn build_node(
    clade: &Clade,
    clade_by_leaves: &HashMap<Vec<String>, &Clade>,
    clades_by_size: &HashMap<usize, Vec<&Clade>>,
) -> Node {
    if clade.leaves.len() == 1 {
        // Leaf
        return Node {
            name: Some(clade.leaves[0].clone()),
            length: clade.length,
            children: vec![],
        };
    }

    let mut children = Vec::new();
    let mut covered = HashSet::new();

    // Find maximal proper subclades of this clade.
    // Optimized by only considering clades smaller than current clade
    for size in 1..clade.leaves.len() {
        if let Some(candidates) = clades_by_size.get(&size) {
            for &sub_clade in candidates {
                if sub_clade.is_subset_of_sorted(&clade.leaves) {
                    let is_maximal =
                        !is_contained_in_larger_subclade(sub_clade, clade, clades_by_size);

                    if is_maximal {
                        children.push(build_node(sub_clade, clade_by_leaves, clades_by_size));
                        covered.extend(sub_clade.leaves.iter().cloned());
                    }
                }
            }
        }
    }

    // Add uncovered leaves as direct children
    for taxon in &clade.leaves {
        if !covered.contains(taxon) {
            // direct lookup
            let leaf_length = if let Some(leaf_clade) = clade_by_leaves.get(&vec![taxon.clone()]) {
                leaf_clade.length
            } else {
                None
            };

            children.push(Node {
                name: Some(taxon.clone()),
                length: leaf_length,
                children: vec![],
            });
        }
    }

    // sort children for deterministic output
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

/// Check if clade is contained in a larger subclade
/// Uses size-based filtering to reduce search space
fn is_contained_in_larger_subclade(
    sub_clade: &Clade,
    parent_clade: &Clade,
    clades_by_size: &HashMap<usize, Vec<&Clade>>,
) -> bool {
    // only check clades larger than sub_clade but smaller than parent
    for size in (sub_clade.leaves.len() + 1)..parent_clade.leaves.len() {
        if let Some(candidates) = clades_by_size.get(&size) {
            for &larger_clade in candidates {
                if larger_clade.is_subset_of_sorted(&parent_clade.leaves)
                    && sub_clade.is_subset_of_sorted(&larger_clade.leaves)
                {
                    return true;
                }
            }
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_clade_creation() {
        let clade = Clade::new(
            vec!["C".to_string(), "A".to_string(), "B".to_string()],
            Some(0.5),
        );

        // Should be sorted and deduplicated
        assert_eq!(clade.leaves, vec!["A", "B", "C"]);
        assert_eq!(clade.length, Some(0.5));

        // Test deduplication
        let clade_with_dups = Clade::new(
            vec!["A".to_string(), "B".to_string(), "A".to_string()],
            Some(0.3),
        );
        assert_eq!(clade_with_dups.leaves, vec!["A", "B"]);
    }

    #[test]
    fn test_newick_parsing_simple() {
        let newick = "(A:0.1,B:0.2);";
        let tree = Tree::from_newick(newick).expect("Failed to parse simple tree");

        assert!(tree.root.name.is_none());
        assert_eq!(tree.root.children.len(), 2);

        let child1 = &tree.root.children[0];
        let child2 = &tree.root.children[1];

        assert_eq!(child1.name, Some("A".to_string()));
        assert_eq!(child1.length, Some(0.1));
        assert_eq!(child2.name, Some("B".to_string()));
        assert_eq!(child2.length, Some(0.2));
    }

    #[test]
    fn test_newick_parsing_complex() {
        let newick = "((A:0.1,B:0.2):0.5,(C:0.3,D:0.4):0.6);";
        let tree = Tree::from_newick(newick).expect("Failed to parse complex tree");

        assert_eq!(tree.root.children.len(), 2);

        let left_internal = &tree.root.children[0];
        let right_internal = &tree.root.children[1];

        assert_eq!(left_internal.length, Some(0.5));
        assert_eq!(left_internal.children.len(), 2);
        assert_eq!(right_internal.length, Some(0.6));
        assert_eq!(right_internal.children.len(), 2);
    }

    #[test]
    fn test_subset_checking_performance() {
        let small_clade = Clade::new(vec!["A".to_string(), "B".to_string()], None);
        let large_set = (0..1000)
            .map(|i| format!("taxon_{}", i))
            .collect::<Vec<_>>();
        let mut large_set_with_ab = large_set.clone();
        large_set_with_ab.extend(vec!["A".to_string(), "B".to_string()]);
        large_set_with_ab.sort();

        // This should be much faster than HashSet-based checking
        assert!(small_clade.is_subset_of_sorted(&large_set_with_ab));
        assert!(!small_clade.is_subset_of_sorted(&large_set));
    }

    #[test]
    fn test_optimized_tree_building() {
        let clades = vec![
            Clade::new(vec!["A".to_string()], Some(0.1)),
            Clade::new(vec!["B".to_string()], Some(0.2)),
            Clade::new(vec!["C".to_string()], Some(0.3)),
            Clade::new(vec!["A".to_string(), "B".to_string()], Some(0.5)),
            Clade::new(
                vec!["A".to_string(), "B".to_string(), "C".to_string()],
                Some(0.8),
            ),
        ];

        let all_leaves = vec!["A".to_string(), "B".to_string(), "C".to_string()];
        let tree = build_tree_from_clades(&all_leaves, &clades);

        println!("Optimized tree: {}", tree.to_newick());

        // Should produce: ((A:0.1,B:0.2):0.5,C:0.3):0.8;
        assert_eq!(tree.root.children.len(), 2); // AB subtree and C
    }
}
