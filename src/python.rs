//! The Python bindings: one function, `_compute`, that takes the keywords of
//! `cygv.Problem` and of its `gv` and `gw` methods as a dictionary, which the Python wrapper
//! (`python/cygv/hkty.py`) has already turned into lists, integers and strings.
//!
//! The computation runs on a thread of its own with the GIL released, so that other Python
//! threads go on meanwhile, while the calling thread waits for it, polling for signals: on
//! ctrl+c it interrupts the computation ([`Interrupt`]) and raises the `KeyboardInterrupt`
//! once it has stopped. The progress goes to Python's `sys.stderr`, which Jupyter shows,
//! through the same polls. The results come back as numpy arrays that own the memory of
//! the Rust vectors they were built in, without a copy.

use crate::cone::ConeError;
use crate::modular::Error as ModularError;
use crate::reference::{self, series_inversion::error::SeriesInversionError};
use crate::{
    compute, Ambient, CandidateCone, Computed, Error, Geometry, Interrupt, InvariantKind,
    Invariants, Method, MoriCone, Options, Problem, Selection, Stratum,
};
use malachite::{Integer, Rational};
use nalgebra::DMatrix;
use numpy::{IntoPyArray, PyArrayMethods};
use pyo3::conversion::FromPyObjectOwned;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use std::collections::HashMap;
use std::sync::mpsc::{self, RecvTimeoutError};
use std::sync::Mutex;
use std::time::{Duration, Instant};

pyo3::create_exception!(
    cygv,
    ComputationError,
    PyRuntimeError,
    "A valid request whose computation failed, e.g. because the invariants did not stabilize \
     as integers or the cones were too large."
);

/// How often the wait for a computation looks for signals and forwards the progress.
const POLL: Duration = Duration::from_millis(100);

/// The stack of the thread of a computation, as large as that of a main thread.
const STACK: usize = 8 << 20;

fn value_error(msg: impl Into<String>) -> PyErr {
    PyValueError::new_err(msg.into())
}

/// The value of a keyword, treating `None` as left out.
fn get<'py, T: FromPyObjectOwned<'py>>(
    args: &Bound<'py, PyDict>,
    key: &str,
) -> PyResult<Option<T>> {
    match args.get_item(key)? {
        Some(v) if !v.is_none() => v.extract::<T>().map(Some).map_err(Into::into),
        _ => Ok(None),
    }
}

/// Turn a list of columns into a matrix, requiring the columns to have equal lengths. An
/// empty list gives `n_rows` rows and no columns.
fn to_matrix(m: Vec<Vec<i32>>, name: &str, n_rows: usize) -> PyResult<DMatrix<i32>> {
    let n_rows = m.first().map_or(n_rows, Vec::len);
    if m.iter().any(|col| col.len() != n_rows) {
        return Err(value_error(format!(
            "the vectors of {name} must all have the same length"
        )));
    }
    Ok(DMatrix::from_iterator(
        n_rows,
        m.len(),
        m.into_iter().flatten(),
    ))
}

/// The problem, the options and the kind of invariants from the keywords.
fn parse(args: &Bound<'_, PyDict>) -> PyResult<(Problem, Options)> {
    let q: Vec<Vec<i32>> = get(args, "q")?.ok_or_else(|| value_error("q is required"))?;
    let q = to_matrix(q, "q", 0)?;
    let h11 = q.ncols();
    let matrix = |key: &str| -> PyResult<Option<DMatrix<i32>>> {
        get::<Vec<Vec<i32>>>(args, key)?
            .map(|m| to_matrix(m, key, h11))
            .transpose()
    };
    let intnums: HashMap<(usize, usize, usize), i32> =
        get(args, "intnums")?.ok_or_else(|| value_error("intnums is required"))?;
    let nefpart: Vec<Vec<usize>> = get(args, "nefpart")?.unwrap_or_default();
    let fan: Option<Vec<Vec<usize>>> = get(args, "fan")?;
    // Each stratum as its divisors and its class, a fraction per entry.
    type Strata = Vec<(Vec<usize>, Vec<(i64, i64)>)>;
    let strata: Option<Strata> = get(args, "strata")?;
    let ambient = match (fan, strata) {
        (Some(_), Some(_)) => {
            return Err(value_error(
                "give either the fan or the strata, which the fan gives, not both",
            ))
        }
        (Some(fan), None) => Ambient::Fan(fan),
        (None, Some(strata)) => Ambient::Strata(
            strata
                .into_iter()
                .map(|(divisors, class)| {
                    let divisors = divisors.try_into().map_err(|_| {
                        value_error("the divisors of a stratum must be three indices")
                    })?;
                    let class = class
                        .into_iter()
                        .map(|(n, d)| {
                            if d == 0 {
                                return Err(value_error(
                                    "the class of a stratum has a zero denominator",
                                ));
                            }
                            Ok(Rational::from(n) / Rational::from(d))
                        })
                        .collect::<PyResult<_>>()?;
                    Ok(Stratum { divisors, class })
                })
                .collect::<PyResult<_>>()?,
        ),
        (None, None) => Ambient::Plain,
    };
    let cone = match (
        matrix("generators")?,
        matrix("mori_hyperplanes")?,
        matrix("mori_equations")?,
    ) {
        (Some(_), Some(_), _) | (Some(_), _, Some(_)) => {
            return Err(value_error(
                "give either the generators or the inequalities of the Mori cone \
                 (mori_hyperplanes and mori_equations), not both",
            ))
        }
        (Some(g), None, None) => MoriCone::Generators(g),
        (None, Some(hyperplanes), equations) => MoriCone::Inequalities {
            hyperplanes,
            equations: equations.unwrap_or_else(|| DMatrix::zeros(h11, 0)),
        },
        (None, None, Some(_)) => return Err(value_error("mori_equations needs mori_hyperplanes")),
        (None, None, None) if matches!(ambient, Ambient::Fan(_)) => MoriCone::FromFan,
        (None, None, None) => {
            return Err(value_error(
                "the Mori cone is needed: give generators or mori_hyperplanes, or a fan",
            ))
        }
    };
    let max_deg: Option<u32> = get(args, "max_deg")?;
    let selection =
        match (
            matrix("target_points")?,
            get::<u32>(args, "min_points")?,
            matrix("curve_classes")?,
        ) {
            (Some(points), None, None) => Selection::TargetPoints { points, max_deg },
            (None, None, None) => Selection::MaxDeg(max_deg.ok_or_else(|| {
                value_error(
                    "which curve classes: give max_deg, target_points, min_points or curve_classes",
                )
            })?),
            (None, Some(n), None) if max_deg.is_none() => Selection::MinPoints(n),
            (None, None, Some(classes)) if max_deg.is_none() => Selection::CurveClasses(classes),
            _ => return Err(value_error(
                "give only one of max_deg, target_points, min_points and curve_classes (max_deg \
                 may accompany target_points)",
            )),
        };
    let method = match get::<String>(args, "method")? {
        None => Method::default(),
        Some(name) => Method::from_name(&name).ok_or_else(|| {
            value_error(format!(
                "unknown method \"{name}\", expected \"auto\", \"modular\" or \"reference\""
            ))
        })?,
    };
    // Each cone as its divisors and the columns of its Hilbert basis.
    type Bases = Vec<(Vec<usize>, Vec<Vec<i32>>)>;
    let hilbert_bases: Option<Bases> = get(args, "hilbert_bases")?;
    let hilbert_bases = hilbert_bases
        .map(|bases| {
            bases
                .into_iter()
                .map(|(divisors, basis)| {
                    Ok(CandidateCone {
                        divisors,
                        hilbert_basis: to_matrix(basis, "hilbert_bases", h11)?,
                    })
                })
                .collect::<PyResult<_>>()
        })
        .transpose()?;
    let problem = Problem {
        geometry: Geometry {
            q,
            intnums,
            nefpart,
            ambient,
        },
        cone,
        selection,
        grading_vector: get(args, "grading_vector")?,
    };
    let options = Options {
        method,
        n_threads: get(args, "n_threads")?,
        certify: get(args, "certify")?.unwrap_or(false),
        low_memory: get(args, "low_memory")?.unwrap_or(false),
        hilbert_bases,
        ..Options::default()
    };
    Ok((problem, options))
}

/// The Python exception for an error: a `ValueError` for invalid or unsupported input, and a
/// `ComputationError` for a valid request whose computation failed.
fn to_py_err(e: Error) -> PyErr {
    let computation = matches!(
        e,
        Error::MoriCone(ConeError::TooLarge | ConeError::Parallelepiped(_))
            | Error::Modular(
                ModularError::NotStable { .. }
                    | ModularError::NotInteger { .. }
                    | ModularError::CertifyPrimes { .. }
                    | ModularError::NotInvertible
            )
            | Error::Reference(reference::Error::SeriesInversion(
                SeriesInversionError::NonIntegerGVError
            ))
    );
    if computation {
        ComputationError::new_err(e.to_string())
    } else {
        value_error(e.to_string())
    }
}

/// The progress of the computations, kept for the waits to write to `sys.stderr`, by the
/// identifier of each computation's interrupt handle (`interrupt::current_id`), which the
/// threads of its pools have too.
struct Progress(Mutex<Option<HashMap<usize, Run>>>);

/// The progress of one computation.
struct Run {
    level: log::LevelFilter,
    start: Instant,
    lines: String,
}

impl log::Log for Progress {
    fn enabled(&self, metadata: &log::Metadata) -> bool {
        metadata.target().starts_with("cygv")
    }

    fn log(&self, record: &log::Record) {
        if !self.enabled(record.metadata()) {
            return;
        }
        let Some(id) = crate::interrupt::current_id() else {
            return;
        };
        if let Some(run) = self.0.lock().unwrap().as_mut().and_then(|r| r.get_mut(&id)) {
            if record.level() <= run.level {
                let secs = run.start.elapsed().as_secs_f64();
                let line = format!("[{secs:9.2}s] {}\n", record.args());
                run.lines.push_str(&line);
            }
        }
    }

    fn flush(&self) {}
}

static PROGRESS: Progress = Progress(Mutex::new(None));

/// Keeps the progress of the computation of the interrupt handle `id` at `verbosity` (see
/// `crate::log_to_stderr`), stamped with the seconds from now, until it is dropped. The level
/// of the process is only ever raised, so that a computation with less progress leaves the
/// others', and those of any other logger, alone.
struct ProgressRun(usize);

impl ProgressRun {
    fn start(id: usize, verbosity: u8) -> Self {
        // Fails only when a logger is already set, which keeps it.
        let _ = log::set_logger(&PROGRESS);
        let level = match verbosity {
            0 => log::LevelFilter::Off,
            1 => log::LevelFilter::Info,
            _ => log::LevelFilter::Debug,
        };
        if level > log::max_level() {
            log::set_max_level(level);
        }
        let run = Run {
            level,
            start: Instant::now(),
            lines: String::new(),
        };
        PROGRESS
            .0
            .lock()
            .unwrap()
            .get_or_insert_with(HashMap::new)
            .insert(id, run);
        Self(id)
    }
}

impl Drop for ProgressRun {
    fn drop(&mut self) {
        if let Some(runs) = PROGRESS.0.lock().unwrap().as_mut() {
            runs.remove(&self.0);
        }
    }
}

/// Writes the progress of the computation `run` kept so far to `sys.stderr`.
fn forward_progress(py: Python<'_>, run: &ProgressRun) {
    let lines = match PROGRESS
        .0
        .lock()
        .unwrap()
        .as_mut()
        .and_then(|r| r.get_mut(&run.0))
    {
        Some(r) => std::mem::take(&mut r.lines),
        None => return,
    };
    if !lines.is_empty() {
        // The progress is not worth an exception of its own, as for a closed stream.
        let _ = py.import("sys").and_then(|sys| {
            let stderr = sys.getattr("stderr")?;
            stderr.call_method1("write", (lines,))?;
            stderr.call_method0("flush")
        });
    }
}

/// The values of the invariants (or the numerators or denominators of GW invariants): in
/// an array when they all fit in 64 bits, as Python integers otherwise.
enum Values {
    Small(Vec<i64>),
    Big(Vec<Integer>),
}

impl Values {
    fn new(values: Vec<Integer>) -> Self {
        let small: Option<Vec<i64>> = values.iter().map(|v| i64::try_from(v).ok()).collect();
        small.map_or(Self::Big(values), Self::Small)
    }

    fn into_py(self, py: Python<'_>) -> PyResult<Bound<'_, PyAny>> {
        Ok(match self {
            Self::Small(v) => v.into_pyarray(py).into_any(),
            Self::Big(v) => PyList::new(py, v)?.into_any(),
        })
    }
}

/// The invariants in columns, built on the thread of the computation.
struct Columns {
    /// The entries of the curve classes, one class after the other.
    curves: Vec<i32>,
    /// The length of a curve class.
    width: usize,
    surfaces: Option<Vec<i64>>,
    degrees: Vec<i64>,
    values: Values,
    /// The denominators of GW invariants.
    denominators: Option<Values>,
    dimension: usize,
    grading_vector: Vec<i32>,
    report: (String, bool, bool, f64),
}

impl Columns {
    /// The columns of `inv`, whose values `split` gives as integers, or as numerators and
    /// denominators.
    fn new<T>(inv: Invariants<T>, split: impl Fn(T) -> (Integer, Option<Integer>)) -> Self {
        let n = inv.entries.len();
        let width = inv.grading_vector.len();
        let nfold = inv.dimension > 3;
        let mut curves = Vec::with_capacity(n * width);
        let mut surfaces = Vec::with_capacity(if nfold { n } else { 0 });
        let mut degrees = Vec::with_capacity(n);
        let mut values = Vec::with_capacity(n);
        let mut denominators = Vec::new();
        let mut fractions = false;
        for e in inv.entries {
            curves.extend_from_slice(&e.curve);
            if let Some(s) = e.surface {
                surfaces.push(s as i64);
            }
            degrees.push(i64::from(e.degree));
            let (v, d) = split(e.value);
            values.push(v);
            fractions = d.is_some();
            denominators.extend(d);
        }
        let report = &inv.report;
        Self {
            curves,
            width,
            surfaces: nfold.then_some(surfaces),
            degrees,
            values: Values::new(values),
            denominators: fractions.then(|| Values::new(denominators)),
            dimension: inv.dimension,
            grading_vector: inv.grading_vector,
            report: (
                report.method.as_str().to_owned(),
                report.certified,
                report.grading_vector_derived,
                report.elapsed.as_secs_f64(),
            ),
        }
    }

    /// The columns as a tuple of `curves`, an array of one row per curve class, `surfaces`
    /// (or `None` for a threefold), `degrees`, `values`, `denominators` (or `None` for GV
    /// invariants), and the dimension, the grading vector and the report.
    fn into_py(self, py: Python<'_>) -> PyResult<Bound<'_, PyTuple>> {
        let n = self.degrees.len();
        let curves = self.curves.into_pyarray(py).reshape([n, self.width])?;
        let surfaces = self.surfaces.map(|s| s.into_pyarray(py));
        let denominators = self.denominators.map(|d| d.into_py(py)).transpose()?;
        (
            curves,
            surfaces,
            self.degrees.into_pyarray(py),
            self.values.into_py(py)?,
            denominators,
            self.dimension,
            self.grading_vector,
            self.report,
        )
            .into_pyobject(py)
    }
}

/// The message of a panic.
fn panic_message(payload: &(dyn std::any::Any + Send)) -> String {
    if let Some(s) = payload.downcast_ref::<&str>() {
        (*s).to_owned()
    } else if let Some(s) = payload.downcast_ref::<String>() {
        s.clone()
    } else {
        "unknown cause".to_owned()
    }
}

/// Compute GV or GW invariants.
#[pyfunction]
#[pyo3(name = "_compute")]
pub fn compute_py<'py>(
    py: Python<'py>,
    kind: &str,
    args: &Bound<'py, PyDict>,
) -> PyResult<Bound<'py, PyTuple>> {
    let kind = match kind {
        "gv" => InvariantKind::GV,
        "gw" => InvariantKind::GW,
        _ => return Err(value_error("the kind must be \"gv\" or \"gw\"")),
    };
    let (problem, mut options) = parse(args)?;
    let verbosity = get::<u8>(args, "verbose")?.unwrap_or(0);
    if options.low_memory {
        crate::return_freed_memory();
    }
    let interrupt = Interrupt::new();
    options.interrupt = Some(interrupt.clone());
    let progress = ProgressRun::start(interrupt.id(), verbosity);
    let (tx, mut rx) = mpsc::channel();
    let worker = std::thread::Builder::new()
        .name("cygv".to_owned())
        .stack_size(STACK)
        .spawn(move || {
            let low_memory = options.low_memory;
            let columns = compute(&problem, kind, &options).map(|c| match c {
                Computed::Gv(gv) => Columns::new(gv, |v| (v, None)),
                Computed::Gw(gw) => Columns::new(gw, |v: Rational| {
                    let negative = v < 0u32;
                    let (n, d) = v.into_numerator_and_denominator();
                    (
                        Integer::from_sign_and_abs(!negative, n),
                        Some(Integer::from(d)),
                    )
                }),
            });
            drop((problem, options));
            if low_memory {
                // Several computations at once share the setting, which then only keeps a
                // later one's memory a little longer.
                crate::keep_freed_memory();
            }
            crate::trim_freed_memory();
            // The receiver only goes away with a panic of its own.
            let _ = tx.send(columns);
        })
        .map_err(|e| PyRuntimeError::new_err(format!("failed to start the computation: {e}")))?;
    let result = loop {
        // The receiver goes to the wait and back, since it cannot be shared.
        let (received, back) = py.detach(move || (rx.recv_timeout(POLL), rx));
        rx = back;
        forward_progress(py, &progress);
        match received {
            Ok(result) => break result,
            Err(RecvTimeoutError::Timeout) => {
                if let Err(e) = py.check_signals() {
                    // The computation stops at its next checkpoint; it is waited for, so that
                    // its threads and its memory are gone when the exception is raised, unless
                    // a second interrupt comes first, which leaves it to stop by itself.
                    interrupt.interrupt();
                    loop {
                        let (received, back) = py.detach(move || (rx.recv_timeout(POLL), rx));
                        rx = back;
                        forward_progress(py, &progress);
                        if !matches!(received, Err(RecvTimeoutError::Timeout)) {
                            let _ = worker.join();
                            break;
                        }
                        if py.check_signals().is_err() {
                            break;
                        }
                    }
                    return Err(e);
                }
            }
            Err(RecvTimeoutError::Disconnected) => {
                let msg = match worker.join() {
                    Err(payload) => panic_message(&*payload),
                    Ok(()) => "it ended without a result".to_owned(),
                };
                return Err(PyRuntimeError::new_err(format!(
                    "the computation failed: {msg}"
                )));
            }
        }
    };
    let _ = worker.join();
    result.map_err(to_py_err)?.into_py(py)
}

/// A Python module implemented in Rust.
#[pymodule(gil_used = false)]
pub fn cygv(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compute_py, m)?)?;
    m.add("ComputationError", m.py().get_type::<ComputationError>())?;
    Ok(())
}
