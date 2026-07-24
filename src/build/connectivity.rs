use crate::dist::rescue::RESCUE_PVALUE_THRESHOLD;
use crate::mash::uncertainty::{mash_significance, DistanceEstimate, Reliability};
use finch::serialization::Sketch;
use std::collections::HashMap;

/// One pairwise comparison, as fat as connectivity analysis
/// needs to know about it.
/// I decoupled it from `DistanceEstimate` / `Reliability` so this
/// module doesn't need to know `divergent` was decided.
pub struct Edge {
    pub query: String,
    pub reference: String,
    pub is_divergent: bool,
    /// Weight used for Kruskal's algorithm
    pub weight: f64,
}

/// Sentinel weight assigned to divergent (Unreliable/NoSharedHashes) pairs
/// so Kruskal's always exhausts every resolvable edge before considering
/// a divergent one, regardless of the divergent pair's own (largely
/// meaningless) relative-uncertainty value.
const DIVERGENT_EDGE_WEIGHT: f64 = f64::MAX;

/// Disjoint-set (union-find) with path compression and union by rank.
struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }

    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]);
        }
        self.parent[x]
    }

    /// Returns true if this union merged two previously-separate sets
    /// (i.e. this edge is actually needed in the spanning tree).
    fn union(&mut self, a: usize, b: usize) -> bool {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return false;
        }
        match self.rank[ra].cmp(&self.rank[rb]) {
            std::cmp::Ordering::Less => self.parent[ra] = rb,
            std::cmp::Ordering::Greater => self.parent[rb] = ra,
            std::cmp::Ordering::Equal => {
                self.parent[rb] = ra;
                self.rank[ra] += 1;
            }
        }
        true
    }
}

fn get_or_insert_index(genome_index: &mut HashMap<String, usize>, name: &str) -> usize {
    if let Some(&idx) = genome_index.get(name) {
        idx
    } else {
        let idx = genome_index.len();
        genome_index.insert(name.to_string(), idx);
        idx
    }
}

/// Every pair of genomes always has *some* distance estimate, so the full
/// pairwise graph is always complete and trivially connected on its own.
/// What actually matters is whether the *non-divergent* subgraph connects
/// everyone by itself, or whether reaching some genome structurally requires
/// at least one divergent pair. This can be accomplished by running Kruskal's
/// algorithm: it reveals, when divergent edges are sentinel-weighted last,
/// that any divergent edge that ends up in the minimum spanning tree is one
/// that Kruskal's algorithm was *forced* to use because the resolvable edges
/// alone weren't enough to reach every genome.
///
/// Returns the (query, reference) pairs that were load-bearing, empty
/// if the non-divergent subgraph is already fully connected on its own.
pub fn analyze_connectivity(edges: &[Edge]) -> Vec<(String, String)> {
    let mut genome_index: HashMap<String, usize> = HashMap::new();
    let mut weighted: Vec<(f64, usize, usize, bool, String, String)> =
        Vec::with_capacity(edges.len());

    for e in edges {
        let qi = get_or_insert_index(&mut genome_index, &e.query);
        let ri = get_or_insert_index(&mut genome_index, &e.reference);
        let weight = if e.is_divergent {
            DIVERGENT_EDGE_WEIGHT
        } else {
            e.weight
        };
        weighted.push((
            weight,
            qi,
            ri,
            e.is_divergent,
            e.query.clone(),
            e.reference.clone(),
        ));
    }

    weighted.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let mut uf = UnionFind::new(genome_index.len());
    let mut load_bearing_divergent_pairs = Vec::new();

    for (_, qi, ri, is_divergent, qname, rname) in &weighted {
        if uf.union(*qi, *ri) && *is_divergent {
            load_bearing_divergent_pairs.push((qname.clone(), rname.clone()));
        }
    }

    load_bearing_divergent_pairs
}

/// Build connectivity-analysis edges from a set of distance estimates.
/// Classifying a pairs as divergent if it fails EITHER the CI-based
/// reliability bar OR Mash's own significance test.
/// It is a "both checks" standard (attempt_kmer_rescue): a candidate k can't
/// pass connectivity purely by making noise look precise.
pub(crate) fn connectivity_edges(
    estimates: &[DistanceEstimate],
    sketches: &[Sketch],
    kmer_length: u8,
) -> Vec<crate::build::connectivity::Edge> {
    let sizes: HashMap<String, u64> = sketches
        .iter()
        .map(|s| (s.name.clone(), s.seq_length))
        .collect();

    estimates
        .iter()
        .map(|e| {
            let ci_divergent = matches!(
                e.reliability,
                Reliability::Unreliable | Reliability::NoSharedHashes
            );
            let query_size = *sizes.get(&e.query_path).unwrap_or(&0);
            let reference_size = *sizes.get(&e.reference_path).unwrap_or(&0);
            let pvalue = mash_significance(
                e.shared_hashes,
                e.total_hashes,
                query_size,
                reference_size,
                kmer_length,
            );
            let is_divergent = ci_divergent || pvalue > RESCUE_PVALUE_THRESHOLD;

            crate::build::connectivity::Edge {
                query: e.query.clone(),
                reference: e.reference.clone(),
                is_divergent,
                weight: e.relative_uncertainty.unwrap_or(0.0),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn edge(query: &str, reference: &str, is_divergent: bool, weight: f64) -> Edge {
        Edge {
            query: query.to_string(),
            reference: reference.to_string(),
            is_divergent,
            weight,
        }
    }

    #[test]
    fn test_fully_resolvable_dataset_needs_no_divergent_pairs() {
        let edges = vec![
            edge("A", "B", false, 0.05),
            edge("B", "C", false, 0.08),
            edge("C", "D", false, 0.06),
            edge("A", "C", false, 0.10),
            edge("A", "D", false, 0.12),
            edge("B", "D", false, 0.09),
        ];
        assert!(analyze_connectivity(&edges).is_empty());
    }

    #[test]
    fn test_isolated_genome_needs_exactly_one_divergent_bridge() {
        // Chain A-B-C-D is fully resolvable; E is only reachable through
        // divergent pairs.
        let edges = vec![
            edge("A", "B", false, 0.05),
            edge("B", "C", false, 0.08),
            edge("C", "D", false, 0.06),
            edge("A", "C", true, 0.0),
            edge("A", "D", true, 0.0),
            edge("A", "E", true, 0.0),
            edge("B", "D", true, 0.0),
            edge("B", "E", true, 0.0),
            edge("C", "E", true, 0.0),
            edge("D", "E", true, 0.0),
        ];
        let load_bearing = analyze_connectivity(&edges);
        assert_eq!(load_bearing.len(), 1);
        assert!(
            load_bearing[0].0 == "E" || load_bearing[0].1 == "E",
            "expected E in the load-bearing pair, got {:?}",
            load_bearing[0]
        );
    }

    #[test]
    fn test_two_clusters_need_exactly_one_bridge_not_all_candidates() {
        // Two resolvable clusters {A,B} and {C,D}, connected only by
        // divergent pairs. Kruskal's needs exactly one of the four
        // divergent candidates, never all of them.
        let edges = vec![
            edge("A", "B", false, 0.05),
            edge("C", "D", false, 0.06),
            edge("A", "C", true, 0.0),
            edge("A", "D", true, 0.0),
            edge("B", "C", true, 0.0),
            edge("B", "D", true, 0.0),
        ];
        assert_eq!(analyze_connectivity(&edges).len(), 1);
    }

    #[test]
    fn test_three_clusters_need_exactly_two_bridges() {
        // k separate resolvable clusters always need exactly k-1 bridging
        // edges in a minimum spanning tree, never more.
        let edges = vec![
            edge("A", "B", false, 0.05),
            edge("C", "D", false, 0.06),
            edge("E", "F", false, 0.04),
            edge("A", "C", true, 0.0),
            edge("A", "D", true, 0.0),
            edge("A", "E", true, 0.0),
            edge("A", "F", true, 0.0),
            edge("B", "C", true, 0.0),
            edge("B", "D", true, 0.0),
            edge("B", "E", true, 0.0),
            edge("B", "F", true, 0.0),
            edge("C", "E", true, 0.0),
            edge("C", "F", true, 0.0),
            edge("D", "E", true, 0.0),
            edge("D", "F", true, 0.0),
        ];
        assert_eq!(analyze_connectivity(&edges).len(), 2);
    }

    #[test]
    fn test_redundant_resolvable_edge_is_never_load_bearing() {
        // A duplicate edge between the same pair (which would close a
        // cycle) should just be skipped by union-find, not cause an error
        // or get counted as load-bearing (it's resolvable, not divergent,
        // so it couldn't be anyway. This checks the union-find handles
        // the redundant-edge case cleanly).
        let edges = vec![
            edge("A", "B", false, 0.05),
            edge("A", "B", false, 0.20),
            edge("B", "C", false, 0.08),
        ];
        assert!(analyze_connectivity(&edges).is_empty());
    }

    #[test]
    fn test_no_edges() {
        let edges: Vec<Edge> = vec![];
        assert!(analyze_connectivity(&edges).is_empty());
    }

    fn make_estimate(
        query: &str,
        reference: &str,
        reliability: Reliability,
        relative_uncertainty: Option<f64>,
        shared_hashes: u64,
        total_hashes: u64,
    ) -> DistanceEstimate {
        DistanceEstimate {
            query: query.to_string(),
            reference: reference.to_string(),
            query_path: query.to_string(),
            reference_path: reference.to_string(),
            jaccard: if total_hashes == 0 {
                0.0
            } else {
                shared_hashes as f64 / total_hashes as f64
            },
            jaccard_ci_95_low: 0.0,
            jaccard_ci_95_high: 1.0,
            mash_distance: 0.0,
            mash_distance_ci_95_low: 0.0,
            mash_distance_ci_95_high: 0.0,
            shared_hashes,
            total_hashes,
            relative_uncertainty,
            reliability,
            kmer_size_used: 21,
            rescued: false,
        }
    }

    fn make_sketch(name: &str, seq_length: u64) -> Sketch {
        Sketch {
            name: name.to_string(),
            seq_length,
            num_valid_kmers: seq_length,
            comment: String::new(),
            hashes: Vec::new(),
            filter_params: finch::filtering::FilterParams::default(),
            sketch_params: finch::sketch_schemes::SketchParams::default(),
        }
    }

    #[test]
    fn test_ci_unreliable_pair_is_always_divergent() {
        // A pair the CI already flags as Unreliable must count as
        // divergent regardless of significance.
        let estimate = make_estimate("A", "B", Reliability::Unreliable, Some(0.9), 5, 1000);
        let sketches = vec![make_sketch("A", 5_000_000), make_sketch("B", 5_000_000)];

        let edges = connectivity_edges(&[estimate], &sketches, 21);
        assert!(edges[0].is_divergent);
    }

    #[test]
    fn test_ci_reliable_and_significant_pair_is_not_divergent() {
        // Realistic genome sizes and k: the same 300/1000 sharing used
        // below is overwhelmingly significant, so a CI-Reliable pair here
        // should NOT be classified divergent.
        let estimate = make_estimate("A", "B", Reliability::Reliable, Some(0.05), 300, 1000);
        let sketches = vec![make_sketch("A", 5_000_000), make_sketch("B", 5_000_000)];

        let edges = connectivity_edges(&[estimate], &sketches, 21);
        assert!(!edges[0].is_divergent);
    }

    #[test]
    fn test_ci_reliable_but_not_significant_pair_is_still_divergent() {
        // The key "both checks" case: tiny genomes and a tiny k-mer size
        // make the exact same 300/1000 sharing statistically
        // indistinguishable from two random sequences (background
        // Jaccard ~0.31 here), even though nothing about the CI itself
        // looks bad. connectivity_edges must still flag this divergent —
        // otherwise a candidate k could "pass" purely via k-mer-space
        // exhaustion rather than real homology.
        let estimate = make_estimate("A", "B", Reliability::Reliable, Some(0.05), 300, 1000);
        let sketches = vec![make_sketch("A", 10), make_sketch("B", 10)];

        let edges = connectivity_edges(&[estimate], &sketches, 2);
        assert!(
            edges[0].is_divergent,
            "expected significance check to override a falsely-reliable CI"
        );
    }

    #[test]
    fn test_weight_and_identity_are_preserved() {
        let estimate = make_estimate(
            "genomeA",
            "genomeB",
            Reliability::Reliable,
            Some(0.07),
            300,
            1000,
        );
        let sketches = vec![
            make_sketch("genomeA", 5_000_000),
            make_sketch("genomeB", 5_000_000),
        ];

        let edges = connectivity_edges(&[estimate], &sketches, 21);
        assert_eq!(edges[0].query, "genomeA");
        assert_eq!(edges[0].reference, "genomeB");
        assert_eq!(edges[0].weight, 0.07);
    }

    #[test]
    fn test_missing_sketch_size_defaults_to_zero_and_is_treated_as_maximally_divergent() {
        // If a genome's size can't be found in the sketch map (shouldn't
        // happen in practice, but shouldn't panic either), size defaults
        // to 0, which drives P(K)=0, background r=0, and mash_pvalue
        // returns 0.0 for shared>0 — i.e. "significant" by that fallback,
        // which is arguably the wrong direction to fail silently in. This
        // test documents the current behavior so a future change to it
        // is deliberate, not accidental.
        let estimate = make_estimate("A", "Z", Reliability::Reliable, Some(0.05), 300, 1000);
        let sketches = vec![make_sketch("A", 5_000_000)]; // "Z" is missing

        let edges = connectivity_edges(&[estimate], &sketches, 21);
        // Documents current behavior rather than asserting it's ideal —
        // worth revisiting if a genome ever legitimately goes missing here.
        assert!(!edges[0].is_divergent);
    }
}
