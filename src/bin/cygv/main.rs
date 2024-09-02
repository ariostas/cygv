use clap::Parser;
use std::fs;
use std::io::{self, Read};
use yaml_rust2::YamlLoader;

/// Compute GV and GW invariants of CY manifold specified by the input file
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Input file
    #[arg(short, long, default_value_t = String::from("stdin"))]
    file: String,
    /// Output file
    #[arg(short, long, default_value_t = String::from("stdout"))]
    output: String,
}

fn main() {
    let args = Args::parse();

    let input_str = if args.file == "stdin" {
        let mut input = String::new();
        io::stdin()
            .read_to_string(&mut input)
            .expect("Failed to read from stdin");
        input
    } else {
        fs::read_to_string(&args.file)
            .unwrap_or_else(|_| panic!("Failed to read from file {}", args.file))
    };

    let input_data = YamlLoader::load_from_str(&input_str).expect("Failed to parse YAML");

    for cy in input_data {
        println!("{:?}", cy);

        // get cy data

        // save or print out invariants
    }
}
