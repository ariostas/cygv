use crate::hkty::compute_gvgw_strings;
use nalgebra::{DMatrix, DVector, RowDVector};
use pyo3::prelude::*;
use std::collections::HashMap;
use std::sync::Once;

/// Installing a ctrl-c handler can only be done once per process, so guard it.
static CTRLC_HANDLER: Once = Once::new();

fn to_matrix(m: Vec<Vec<i32>>) -> DMatrix<i32> {
    let n_cols = m.len();
    let n_rows = if let Some(v) = m.first() { v.len() } else { 0 };
    let mut res = DMatrix::zeros(n_rows, n_cols);
    for (col_src, mut col_dst) in m.iter().zip(res.column_iter_mut()) {
        for (src, dst) in col_src.iter().zip(col_dst.iter_mut()) {
            *dst = *src;
        }
    }
    res
}

fn to_vector(v: Vec<i32>) -> DVector<i32> {
    let len = v.len();
    let mut res = DVector::zeros(len);
    for (src, dst) in v.iter().zip(res.iter_mut()) {
        *dst = *src;
    }
    res
}

fn to_rowvector(v: Vec<i32>) -> RowDVector<i32> {
    let len = v.len();
    let mut res = RowDVector::zeros(len);
    for (src, dst) in v.iter().zip(res.iter_mut()) {
        *dst = *src;
    }
    res
}

fn to_vec(v: DVector<i32>) -> Vec<i32> {
    let len = v.len();
    let mut res = vec![0; len];
    for (src, dst) in v.iter().zip(res.iter_mut()) {
        *dst = *src;
    }
    res
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
    CTRLC_HANDLER.call_once(|| {
        ctrlc::set_handler(|| std::process::exit(1)).expect("failed to install ctrl-c handler");
    });
    let generators = to_matrix(generators);
    let grading_vector = to_rowvector(grading_vector);
    let q = to_matrix(q);
    let target_points = target_points.map(to_matrix);
    let nefpart = nefpart.unwrap_or_default();
    let nefpart: Vec<_> = nefpart.into_iter().map(to_vector).collect();

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
        .map(|((v, c), gvgw)| ((to_vec(v), c), gvgw))
        .collect())
}

/// A Python module implemented in Rust.
#[pymodule(gil_used = false)]
pub fn cygv(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compute_gvgw, m)?)?;
    Ok(())
}
