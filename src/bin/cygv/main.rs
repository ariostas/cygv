//! Command line interface to compute GV and GW invariants of CY manifolds.

use clap::Parser;
use cygv::io::Input;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

const LONG_ABOUT: &str = "\
Compute GV and GW invariants of the CY manifolds specified in the input.

The input is a stream of YAML documents, each of which specifies a single CY
manifold, and the results are written back as one YAML document per input
document, sorted by the degree of the curve classes under the grading vector.
A minimal input document looks as follows.

  generators: [[0, -1], [1, 2]]
  grading_vector: [3, -1]
  q: [[1, 1, 1, 0, 1, 2], [0, 0, -1, 1, 1, -1]]
  intnums: [[0, 0, 0, 2], [0, 0, 1, 1], [0, 1, 1, -1], [1, 1, 1, 5]]
  min_points: 100

The required fields are the generators of the Mori cone, a grading vector that
has a positive inner product with every generator, the GLSM charge matrix q, and
the triple intersection numbers, given as a list of [i, j, k, value] entries.
The optional fields are the following.

  name           A label that is echoed back in the output.
  nefpart        The nef partition, given as a list of lists of indices of the
                 rays. It defaults to an empty list, i.e. a hypersurface.
  invariants     Either gv or gw. It defaults to gv.
  max_deg        The maximum degree of the curve classes to compute.
  min_points     The minimum number of curve classes to compute.
  target_points  The curve classes that must be included in the computation.
  is_threefold   Whether the CY is a threefold. It is deduced from q and the nef
                 partition when it is not given.
  prec           The precision, in bits, of the floating-point arithmetic. Exact
                 rational arithmetic is used when it is not given.
  n_threads      The number of threads to use. It defaults to the number of
                 threads that are available.
  pool_size      The number of coefficients kept in the number pools.

At most one of max_deg, min_points, and target_points may be given. When none of
them is given, the generators are taken to be the complete list of curve classes
to use, which must then include the origin.";

/// Compute GV and GW invariants of the CY manifolds specified in the input.
#[derive(Parser, Debug)]
#[command(version, about, long_about = LONG_ABOUT)]
struct Args {
    /// Input file with the YAML data. Reads from stdin when it is not given or is "-"
    #[arg(short, long, value_name = "FILE")]
    file: Option<PathBuf>,
    /// Output file. Writes to stdout when it is not given or is "-"
    #[arg(short, long, value_name = "FILE")]
    output: Option<PathBuf>,
    /// Number of threads to use. It overrides the value given in the input
    #[arg(short = 'j', long, value_name = "N")]
    threads: Option<u32>,
}

/// Whether a path refers to the standard input or output.
fn is_stdio(path: Option<&PathBuf>) -> bool {
    match path {
        None => true,
        Some(p) => p == Path::new("-"),
    }
}

/// Read the input data from a file or from the standard input.
fn read_input(path: Option<&PathBuf>) -> Result<String, Box<dyn Error>> {
    if is_stdio(path) {
        let mut data = String::new();
        io::stdin()
            .read_to_string(&mut data)
            .map_err(|e| format!("failed to read from stdin: {e}"))?;
        return Ok(data);
    }
    let path = path.expect("the path is present when it is not stdin");
    std::fs::read_to_string(path)
        .map_err(|e| format!("failed to read from {}: {e}", path.display()).into())
}

/// Open the output file, or the standard output.
fn open_output(path: Option<&PathBuf>) -> Result<Box<dyn Write>, Box<dyn Error>> {
    if is_stdio(path) {
        return Ok(Box::new(BufWriter::new(io::stdout())));
    }
    let path = path.expect("the path is present when it is not stdout");
    let file =
        File::create(path).map_err(|e| format!("failed to write to {}: {e}", path.display()))?;
    Ok(Box::new(BufWriter::new(file)))
}

fn run(args: &Args) -> Result<(), Box<dyn Error>> {
    let data = read_input(args.file.as_ref())?;
    let mut inputs = Input::load_all(&data)?;
    if let Some(n_threads) = args.threads {
        for input in inputs.iter_mut() {
            input.n_threads = Some(n_threads);
        }
    }

    let mut writer = open_output(args.output.as_ref())?;
    for input in inputs.iter() {
        let results = input.compute();
        input.write_results(&mut writer, &results)?;
        // The results of each document are flushed as they are computed, so that they
        // can be inspected while the remaining ones are still being computed.
        writer.flush()?;
    }
    Ok(())
}

fn main() -> ExitCode {
    let args = Args::parse();
    match run(&args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("Error: {e}");
            ExitCode::FAILURE
        }
    }
}
