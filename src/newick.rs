// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::collections::HashSet;

#[derive(Debug, Clone, PartialEq)]
pub struct Clade {
    pub leaves: Vec<String>, // sorted leaves
    pub length: Option<f64>, // not used for Hash/Eq
}

impl Clade {
    pub fn new(mut leaves: Vec<String>, length: Option<f64>) -> Self {
        leaves.sort(); // enforce deterministic ordering
        leaves.dedup(); // remove duplicates
        Clade { leaves, length }
    }

    pub fn leaf_set(&self) -> HashSet<String> {
        self.leaves.iter().cloned().collect()
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
            // leaf
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

#[cfg(test)]
mod tests {
    use super::*;

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
}
