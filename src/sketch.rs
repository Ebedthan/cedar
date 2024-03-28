// Copyright 2024 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::fs::File;
use std::path::{Path, PathBuf};

use finch::{
    filtering::FilterParams, serialization::write_mash_file, sketch_files,
    sketch_schemes::SketchParams,
};

use finch::errors::FinchResult;

/// Create sketches from fasta/fastq files
pub fn create_sketches(
    filenames: Vec<&str>,
    kmer_size: u8,
    sketch_size: usize,
    oversketch: usize,
    seed: u64,
    outdir: &str,
) -> FinchResult<Vec<String>> {
    // Create SketchParams struct for finch
    let sketch_params = SketchParams::Mash {
        kmers_to_sketch: sketch_size * oversketch,
        final_size: sketch_size,
        no_strict: false,
        kmer_length: kmer_size,
        hash_seed: seed,
    };

    // Create FilterParams struct for finch
    let filter_params: FilterParams = FilterParams {
        filter_on: Some(false),
        abun_filter: (Some(0u32), None),
        err_filter: 1f64,
        strand_filter: 0.1f64,
    };

    // List of created sketch
    let mut sketches_list: Vec<String> = Vec::new();

    // Create directory
    for filename in filenames {
        let sketches = sketch_files(&[filename], &sketch_params, &filter_params)?;
        let mut out = PathBuf::from(outdir);
        out.push(
            Path::new(filename)
                .file_name()
                .unwrap()
                .to_str()
                .unwrap()
                .to_owned()
                + ".msh",
        );
        sketches_list.push(out.clone().into_os_string().into_string().unwrap());
        let mut out = File::create(&out)?;
        write_mash_file(&mut out, &sketches)?;
    }

    Ok(sketches_list)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_create_sketches() {
        // Define test parameters
        let filenames = vec!["test/GCF_000009265.fna", "test/GCA_000005845.fna.xz"];
        let kmer_size = 21;
        let sketch_size = 1000;
        let oversketch = 200;
        let seed = 42;
        let outdir = "test_output";
        fs::create_dir(outdir).unwrap();

        // Call the function under test
        let result = create_sketches(
            filenames.clone(),
            kmer_size,
            sketch_size,
            oversketch,
            seed,
            outdir,
        );
        // Verify that the function returned successfully
        assert!(result.is_ok());

        // Verify that the output directory and sketch files were created
        assert!(fs::metadata(outdir).is_ok());
        for filename in &filenames {
            let output_filename = format!(
                "{}/{}.msh",
                outdir,
                Path::new(filename).file_name().unwrap().to_str().unwrap()
            );
            println!("{output_filename}");
            assert!(fs::metadata(output_filename).is_ok());
        }

        fs::remove_dir_all(outdir).unwrap();
    }
}
