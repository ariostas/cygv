use crate::hkty::compute_gvgw_strings;
use ctrlc;
use nalgebra::{DMatrix, DVector, RowDVector};
use pyo3::prelude::*;
use std::collections::HashMap;

/// Turn a list of columns into a matrix, requiring the columns to have equal
/// lengths.
fn to_matrix(m: Vec<Vec<i32>>, name: &str) -> PyResult<DMatrix<i32>> {
    let n_rows = m.first().map_or(0, Vec::len);
    if m.iter().any(|col| col.len() != n_rows) {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "the vectors of \"{name}\" must all have the same length"
        )));
    }
    Ok(DMatrix::from_iterator(
        n_rows,
        m.len(),
        m.into_iter().flatten(),
    ))
}

/// Compute GV or GW invariants
#[pyfunction]
#[pyo3(name = "_compute_gvgw")]
#[pyo3(signature = (generators, grading_vector, q, intnums, find_gv, is_threefold, max_deg=None, min_points=None, target_points=None, nefpart=None, n_threads=None, pool_size=1000, prec=None))]
#[allow(clippy::type_complexity, clippy::too_many_arguments)]
pub fn compute_gvgw(
    generators: Vec<Vec<i32>>,
    grading_vector: Vec<i32>,
    q: Vec<Vec<i32>>,
    intnums: HashMap<(usize, usize, usize), i32>,
    find_gv: bool,
    is_threefold: bool,
    max_deg: Option<u32>,
    min_points: Option<u32>,
    target_points: Option<Vec<Vec<i32>>>,
    nefpart: Option<Vec<Vec<i32>>>,
    n_threads: Option<u32>,
    pool_size: usize,
    prec: Option<u32>,
) -> PyResult<Vec<((Vec<i32>, usize), String)>> {
    ctrlc::set_handler(|| std::process::exit(1)).unwrap();
    let generators = to_matrix(generators, "generators")?;
    let grading_vector = RowDVector::from_vec(grading_vector);
    let q = to_matrix(q, "q")?;
    let target_points = target_points
        .map(|m| to_matrix(m, "target_points"))
        .transpose()?;
    let nefpart: Vec<_> = nefpart
        .unwrap_or_default()
        .into_iter()
        .map(DVector::from_vec)
        .collect();

    let res = compute_gvgw_strings(
        generators,
        grading_vector,
        q,
        nefpart,
        intnums,
        find_gv,
        is_threefold,
        max_deg,
        min_points,
        target_points,
        n_threads,
        pool_size,
        prec,
    );

    Ok(res
        .into_iter()
        .map(|((v, c), gvgw)| ((v.as_slice().to_vec(), c), gvgw))
        .collect())
}

/// A Python module implemented in Rust.
#[pymodule(gil_used = false)]
pub fn cygv(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compute_gvgw, m)?)?;
    Ok(())
}
