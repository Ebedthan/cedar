use std::fs::{self, File};
use std::path::{Path, PathBuf};

use anyhow::format_err;
use clap::ArgMatches;
use finch::{
    filtering::FilterParams, serialization::write_mash_file, sketch_files,
    sketch_schemes::SketchParams,
};

pub fn generate_sketch_files(matches: &ArgMatches) -> Vec<String> {
    // Read file names
    let filenames = matches
        .get_many::<String>("INPUT")
        .unwrap_or_default()
        .map(|v| v.as_str())
        .collect::<Vec<_>>();

    // Read input parameters
    let kmer_size: u8 = matches.get_one::<String>("kmer").unwrap().parse().unwrap();
    let sketch_size: usize = matches.get_one::<String>("size").unwrap().parse().unwrap();
    let oversketch: usize = matches
        .get_one::<String>("oversketch")
        .unwrap()
        .parse()
        .unwrap();
    let seed: u64 = matches.get_one::<String>("seed").unwrap().parse().unwrap();

    // Sketch parameters
    let sketch_params = SketchParams::Mash {
        kmers_to_sketch: sketch_size * oversketch,
        final_size: sketch_size,
        no_strict: false,
        kmer_length: kmer_size,
        hash_seed: seed,
    };

    // Filter parameters
    let filter_params: FilterParams = FilterParams {
        filter_on: Some(false),
        abun_filter: (Some(0u32), None),
        err_filter: 1f64,
        strand_filter: 0.1f64,
    };

    // Create output directory
    if !Path::new("treetrust").exists() {
        fs::create_dir("treetrust")
            .map_err(|_| format_err!("Could not create directory treetrust"))
            .expect("Directory creation failed");
    }

    // List of created sketch
    let mut sketches_list: Vec<String> = Vec::new();

    // Create directory
    for filename in filenames {
        let sketches = sketch_files(&[filename], &sketch_params, &filter_params)
            .expect("Files sketching failed");
        let mut output = PathBuf::new();
        output.push("treetrust");
        output.push(
            Path::new(filename)
                .file_name()
                .unwrap()
                .to_str()
                .unwrap()
                .to_owned()
                + ".msh",
        );
        sketches_list.push(output.clone().into_os_string().into_string().unwrap());
        let mut out = File::create(&output)
            .map_err(|_| format_err!("Could not open {}", output.display()))
            .expect("Create output file failed");
        write_mash_file(&mut out, &sketches).expect("Write mash file failed");
    }

    sketches_list
}
