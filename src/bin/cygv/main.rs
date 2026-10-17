//! Command line interface to compute GV and GW invariants of CY manifolds.

use clap::Parser;
use cygv::io::Input;
use cygv::{InvariantKind, Method};
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

  q: [[1, 1, 1, 0, 1, 2], [0, 0, -1, 1, 1, -1]]
  intnums: [[0, 0, 0, 2], [0, 0, 1, 1], [0, 1, 1, -1], [1, 1, 1, 5]]
  generators: [[0, -1], [1, 2]]
  max_deg: 10

Every document gives the Calabi-Yau, its Mori cone, and which curve classes to
compute.

The Calabi-Yau:
  q              The GLSM charges, one charge vector per Kähler parameter
                 (required). A divisor index is a position in these vectors.
  intnums        The intersection numbers, as [i, j, k, value] entries or a
                 mapping from \"i,j,k\" to values (required). For an n-fold, i is
                 a reference surface.
  nefpart        The nef partition of a complete intersection, as lists of
                 divisor indices. Left out for a hypersurface.
  fan            The maximal cones of the fan of the ambient toric variety, as
                 lists of divisor indices. It gives the strata of a vex phase,
                 the Mori cone and the default grading vector, and the
                 intersection numbers are checked against it. Threefold
                 hypersurfaces only.
  strata         Instead of the fan, the strata of a vex phase, as mappings such
                 as {divisors: [1, 2, 4], class: [1, \"-1/2\"]}.

The Mori cone, one of (or neither, with a fan):
  generators        Vectors that span it. Every lattice point of the cone they
                    span is used.
  mori_hyperplanes  The inequalities h . C >= 0 that cut it out, with
  mori_equations    the equations e . C = 0 of its span if it is a face.

  grading_vector    The degree is the inner product with it, which must be
                    positive on the Mori cone. It defaults to the primitive sum
                    of the inequalities of the Mori cone.

Which curve classes, one of:
  max_deg        Those of degree at most this.
  target_points  Those below the given curve classes P (the C with P - C in the
                 Mori cone), up to their largest degree, or max_deg if it is
                 also given and lower.
  min_points     The lowest degrees that hold at least this many classes.
  curve_classes  Exactly the classes listed (zero is added if it is missing).

Options:
  name           A label that is echoed back in the output.
  invariants     Either gv or gw. It defaults to gv.
  method         Which pipeline to use: auto, modular or reference. It defaults
                 to auto, which uses the much faster modular pipeline for
                 max_deg and target_points, and the reference one otherwise.
  certify        Whether to prove the results exact: the modular pipeline then
                 bounds every invariant, and takes more primes until the bounds
                 are below half their product. Not in vex phases. It defaults to
                 false; the reference pipeline is exact anyway.
  n_threads      The number of threads to use. It defaults to one per core.
  hilbert_bases  For every set T of one or two divisors, the Hilbert basis of
                 the cone K_T of the curve classes in the Mori cone that only
                 meet the divisors in T negatively, as computed by programs such
                 as Normaliz, e.g. {divisors: [0, 5], hilbert_basis: [[1, 0],
                 [0, 1]]}, and in a vex phase also for the three divisors of
                 every stratum. The modular pipeline computes them otherwise
                 (see examples/cones.yaml).

The output gives the dimension of the CY, the grading vector, and a report of
the pipeline that ran and whether the results are certified.";

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
    /// Maximum degree, for every input. It overrides the value given in the input, and
    /// can only be used with max_deg or target_points
    #[arg(short = 'd', long, value_name = "N")]
    max_deg: Option<u32>,
    /// Which invariants to compute, gv or gw, for every input
    #[arg(short, long, value_parser = parse_invariants)]
    invariants: Option<InvariantKind>,
    /// Prove the results exact, for every input
    #[arg(long, conflicts_with = "no_certify")]
    certify: bool,
    /// Do not prove the results exact, for every input
    #[arg(long)]
    no_certify: bool,
}

/// Parse the name of a kind of invariants.
fn parse_invariants(name: &str) -> Result<InvariantKind, String> {
    match name.to_lowercase().as_str() {
        "gv" => Ok(InvariantKind::GV),
        "gw" => Ok(InvariantKind::GW),
        _ => Err("expected gv or gw".to_owned()),
    }
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
            input.options.n_threads = Some(n_threads as usize);
        }
        if let Some(method) = args.method {
            input.options.method = method;
        }
        if let Some(max_deg) = args.max_deg {
            input.set_max_deg(max_deg)?;
        }
        if let Some(kind) = args.invariants {
            input.kind = kind;
        }
        if args.certify || args.no_certify {
            input.options.certify = args.certify;
        }
        // Every document is checked before any is computed.
        input.validate()?;
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
