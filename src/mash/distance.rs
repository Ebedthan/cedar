use finch::distance::distance;
use finch::serialization::Sketch;
use finch::serialization::SketchDistance;
use itertools::Itertools;
use std::path::Path;

/// Cached basename extraction
pub fn extract_basename(path: &str) -> String {
    Path::new(path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_string()
}

/// Compute distance between sketches
/// Uses rayon for parallel processing
pub fn compute_distances(sketches: Vec<Sketch>) -> Vec<SketchDistance> {
    use rayon::prelude::*;
    sketches
        .into_iter()
        .combinations_with_replacement(2)
        .collect::<Vec<_>>() // collect combinations first
        .into_par_iter() // then parallelize distance computation
        .filter_map(|pair| {
            let dist = distance(&pair[0], &pair[1], false).ok()?;
            // early validation to avoid storing invalid distance
            (dist.mash_distance <= 1.0 && dist.mash_distance >= 0.0).then_some(dist)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use finch::filtering::FilterParams;
    use finch::sketch_schemes::{KmerCount, SketchParams};
    use std::fs;

    #[test]
    fn test_extract_basename() {
        assert_eq!(extract_basename("/data/genomes/A.fna"), "A");
        assert_eq!(extract_basename("A.fna"), "A");
        assert_eq!(extract_basename("/data/genomes/A.fasta.gz"), "A.fasta"); // only the last extension is stripped
        assert_eq!(extract_basename("relative/path/B"), "B"); // no extension at all
        assert_eq!(extract_basename(""), "");
        assert_eq!(extract_basename("/trailing/slash/"), "slash");
    }

    /// Same construction verified against real `finch::distance::distance()`
    /// in earlier tests: with matching sketch sizes, `shared_range` produces
    /// exactly `|shared_range| / (|shared_range| + unique_count)` Jaccard
    /// between two genomes sharing the same range, as long as each genome's
    /// `unique_offset` range never collides with another's.
    fn make_sketch(
        name: &str,
        shared_range: std::ops::Range<u64>,
        unique_offset: u64,
        unique_count: u64,
    ) -> Sketch {
        let mut hashes: Vec<KmerCount> = shared_range
            .map(|h| KmerCount {
                hash: h,
                kmer: vec![],
                count: 1,
                extra_count: 0,
                label: None,
            })
            .collect();
        hashes.extend((0..unique_count).map(|i| KmerCount {
            hash: unique_offset + i,
            kmer: vec![],
            count: 1,
            extra_count: 0,
            label: None,
        }));
        hashes.sort_by_key(|k| k.hash);
        Sketch {
            name: name.to_string(),
            seq_length: 5_000_000,
            num_valid_kmers: hashes.len() as u64,
            comment: String::new(),
            hashes,
            filter_params: FilterParams::default(),
            sketch_params: SketchParams::Mash {
                kmers_to_sketch: 1000,
                final_size: 1000,
                no_strict: false,
                kmer_length: 21,
                hash_seed: 42,
            },
        }
    }

    #[test]
    fn test_compute_distances_synthetic() {
        // Doesn't depend on the (currently missing) test/sketches fixtures.
        let a = make_sketch("A", 0..950, 1_000_000, 50);
        let b = make_sketch("B", 0..950, 2_000_000, 50);

        let distances = compute_distances(vec![a, b]);

        // combinations_with_replacement(2) over 2 sketches: self-A, self-B,
        // and A-vs-B -- 3 entries total, matching the fixture-based test's
        // expected count for a 2-genome directory.
        assert_eq!(distances.len(), 3);

        let cross = distances
            .iter()
            .find(|d| d.query != d.reference)
            .expect("expected exactly one cross-pair");
        assert!((cross.jaccard - 0.95).abs() < 1e-9);
        assert_eq!(cross.common_hashes, 950);
        assert_eq!(cross.total_hashes, 1000);

        // Self-pairs are always identical (jaccard=1.0, distance=0.0).
        for d in distances.iter().filter(|d| d.query == d.reference) {
            assert_eq!(d.jaccard, 1.0);
            assert_eq!(d.mash_distance, 0.0);
        }
    }

    // Test compute_distances function
    #[test]
    fn test_compute_distances() {
        let mut sketches = Vec::new();
        for file in fs::read_dir("test/sketches").unwrap() {
            sketches.push(finch::open_sketch_file(file.unwrap().path()).unwrap());
        }
        let distances = compute_distances(sketches.into_iter().flatten().collect_vec());
        assert_eq!(distances.len(), 3);
        for distance in &distances {
            assert!(distance.mash_distance <= 1.0);
        }
    }
}
