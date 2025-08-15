// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use finch::serialization::Sketch;
use std::collections::{HashMap, HashSet};

use rand::{rng, seq::IndexedRandom};

#[derive(Debug, Clone)]
pub struct Tree {
    pub root: Node,
}

#[derive(Debug, Clone)]
pub struct Node {
    pub name: Option<String>,
    pub length: Option<f64>,
    pub children: Vec<Node>,
}

impl Tree {
    pub fn from_newick(s: &str) -> Result<Self, String> {
        let chars: Vec<char> = s.chars().collect();
        let mut pos = 0;
        let root = parse_subtree(&chars, &mut pos)?;
        Ok(Self { root })
    }
}

impl Node {
    pub fn new(name: Option<String>, length: Option<f64>) -> Self {
        Self {
            name,
            length,
            children: Vec::new(),
        }
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

fn get_clades(node: &Node) -> (HashSet<String>, Vec<(HashSet<String>, f64)>) {
    if node.children.is_empty() {
        let mut set = HashSet::new();
        set.insert(node.name.clone().unwrap());
        return (set, vec![]);
    }
    let mut taxa = HashSet::new();
    let mut clades = vec![];
    for child in &node.children {
        let (child_taxa, mut child_clades) = get_clades(child);
        taxa.extend(child_taxa.clone());
        clades.append(&mut child_clades);
    }
    clades.push((taxa.clone(), node.length.unwrap_or(0.0)));
    (taxa, clades)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_clades() {
        let input_str = "((A:0.1,B:0.2),C:0.6);";
        match Tree::from_newick(input_str) {
            Ok(tree) => println!("{:#?}", tree),
            Err(e) => eprintln!("Error: {}", e),
        }
    }
}
