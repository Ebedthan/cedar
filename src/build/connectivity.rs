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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nwk::Clade; // unused placeholder import removed below if not needed
}
