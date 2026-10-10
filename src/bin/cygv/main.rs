//! Command line interface to compute GV and GW invariants of CY manifolds.

use clap::Parser;
use cygv::io::Input;
use cygv::Method;
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
  target_points  Curve classes P whose backward lightcones are computed: the
                 classes C such that P - C is in the Mori cone.
  is_threefold   Whether the CY is a threefold. It is deduced from q and the nef
                 partition when it is not given.
  n_threads      The number of threads to use. It defaults to the number of
                 threads that are available.
  strata         The strata of a vex phase of the toric variety, a list of
                 mappings with the three divisors of a stratum and its class,
                 e.g. {divisors: [1, 2, 4], class: [1, \"-1/2\"]}. Only for
                 threefold hypersurfaces.
  method         Which pipeline to use: auto, modular or reference. It defaults
                 to auto. The modular pipeline is much faster, but does not
                 support min_points or strata; auto uses it when it can, and the
                 reference pipeline otherwise.

The modular pipeline is much faster with cone data, which programs such as
Normaliz compute (see examples/cones.yaml). With it the generators can be left
out, but then only the modular pipeline can run.

  mori_hyperplanes  The inequalities h . C >= 0 of the Mori cone.
  mori_equations    The equations e . C = 0 of its span, when it is a face of the
                    Mori cone.
  cones             For every set T of one or two divisors, the Hilbert basis of
                    the cone K_T of the curve classes in the Mori cone that only
                    meet the divisors in T negatively, as a list of mappings,
                    e.g. {divisors: [0, 5], hilbert_basis: [[1, 0], [0, 1]]}.
                    The elements above max_deg can be left out.

min_points cannot be combined with max_deg or target_points. max_deg can be
combined with target_points, and then lowers their largest degree, which is the
maximum degree otherwise. When none of them is given, the generators are taken
to be the complete list of curve classes to use, which must then include the
origin.";

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
    /// Pipeline to use. It overrides the value given in the input
    #[arg(short, long, value_parser = parse_method)]
    method: Option<Method>,
}

/// Parse the name of a method.
fn parse_method(name: &str) -> Result<Method, String> {
    Method::from_name(name).ok_or_else(|| "expected auto, modular or reference".to_owned())
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
    for input in inputs.iter_mut() {
        if let Some(n_threads) = args.threads {
            input.n_threads = Some(n_threads);
        }
        if let Some(method) = args.method {
            input.method = method;
        }
    }

    let mut writer = open_output(args.output.as_ref())?;
    for input in inputs.iter() {
        let results = input.compute()?;
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
