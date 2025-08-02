// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use finch::serialization::Sketch;
use rand::{rng, seq::IndexedRandom};

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
