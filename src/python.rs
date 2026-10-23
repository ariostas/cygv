//! The Python bindings: one function, `_compute`, that takes the keywords of
//! `cygv.Problem` and of its `gv` and `gw` methods as a dictionary, which the Python wrapper
//! (`python/cygv/hkty.py`) has already turned into lists, integers and strings.

use crate::cone::ConeError;
use crate::modular::Error as ModularError;
use crate::reference::{self, series_inversion::error::SeriesInversionError};
use crate::{
    compute, Ambient, CandidateCone, Computed, Error, Geometry, InvariantKind, Invariants, Method,
    MoriCone, Options, Problem, Selection, Stratum,
};
use malachite::Rational;
use nalgebra::DMatrix;
use pyo3::conversion::FromPyObjectOwned;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::HashMap;
use std::sync::Once;

pyo3::create_exception!(
    cygv,
    ComputationError,
    PyRuntimeError,
    "A valid request whose computation failed, e.g. because the invariants did not stabilize \
     as integers or the cones were too large."
);

/// Installing a ctrl-c handler can only be done once per process, so guard it.
static CTRLC_HANDLER: Once = Once::new();

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

/// The entries of the invariants, with their values as strings, the dimension, the grading
/// vector and the report.
type Output = (
    Vec<((Vec<i32>, Option<usize>), String)>,
    usize,
    Vec<i32>,
    (String, bool, bool, f64),
);

fn output<T: ToString>(inv: Invariants<T>) -> Output {
    let report = &inv.report;
    (
        inv.entries
            .into_iter()
            .map(|e| ((e.curve, e.surface), e.value.to_string()))
            .collect(),
        inv.dimension,
        inv.grading_vector,
        (
            report.method.as_str().to_owned(),
            report.certified,
            report.grading_vector_derived,
            report.elapsed.as_secs_f64(),
        ),
    )
}

/// Compute GV or GW invariants.
#[pyfunction]
#[pyo3(name = "_compute")]
pub fn compute_py(kind: &str, args: &Bound<'_, PyDict>) -> PyResult<Output> {
    CTRLC_HANDLER.call_once(|| {
        ctrlc::set_handler(|| std::process::exit(1)).expect("failed to install ctrl-c handler");
    });
    let kind = match kind {
        "gv" => InvariantKind::GV,
        "gw" => InvariantKind::GW,
        _ => return Err(value_error("the kind must be \"gv\" or \"gw\"")),
    };
    let (problem, options) = parse(args)?;
    // The computation runs in a process of its own (see hkty.py).
    if let Some(verbosity) = get::<u8>(args, "verbose")? {
        crate::log_to_stderr(verbosity);
    }
    if options.low_memory {
        crate::return_freed_memory();
    }
    Ok(
        match compute(&problem, kind, &options).map_err(to_py_err)? {
            Computed::Gv(gv) => output(gv),
            Computed::Gw(gw) => output(gw),
        },
    )
}

/// A Python module implemented in Rust.
#[pymodule(gil_used = false)]
pub fn cygv(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compute_py, m)?)?;
    m.add("ComputationError", m.py().get_type::<ComputationError>())?;
    Ok(())
}
