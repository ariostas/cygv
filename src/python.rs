use crate::run_hkty;
use ctrlc;
use nalgebra::{DMatrix, DVector, RowDVector};
use pyo3::prelude::*;
use rug::{ops::PowAssign, Float, Rational};
use std::collections::HashMap;

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
#[pyo3(signature = (generators, grading_vector, q, intnums, find_gv, is_threefold, max_deg=None, min_points=None, nefpart=None, n_threads=None, pool_size=1000, prec=None))]
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
    nefpart: Option<Vec<Vec<i32>>>,
    n_threads: Option<u32>,
    pool_size: usize,
    prec: Option<u32>,
) -> PyResult<Vec<((Vec<i32>, usize), String)>> {
    ctrlc::set_handler(|| std::process::exit(1)).unwrap();
    let generators = to_matrix(generators);
    let grading_vector = to_rowvector(grading_vector);
    let q = to_matrix(q);
    let nefpart = nefpart.unwrap_or_default();
    let nefpart: Vec<_> = nefpart.into_iter().map(to_vector).collect();

    let final_res;
    if let Some(n_bits) = prec {
        let mut zero_cutoff = Float::with_val(n_bits, 10);
        zero_cutoff.pow_assign(-(n_bits as i32) / 3);
        let res = match (find_gv, is_threefold) {
            (true, true) => run_hkty::<Float, true, true>(
                generators,
                grading_vector,
                zero_cutoff,
                max_deg,
                min_points,
                q,
                nefpart,
                intnums,
                n_threads,
                pool_size,
            ),
            (false, true) => run_hkty::<Float, false, true>(
                generators,
                grading_vector,
                zero_cutoff,
                max_deg,
                min_points,
                q,
                nefpart,
                intnums,
                n_threads,
                pool_size,
            ),
            (true, false) => run_hkty::<Float, true, false>(
                generators,
                grading_vector,
                zero_cutoff,
                max_deg,
                min_points,
                q,
                nefpart,
                intnums,
                n_threads,
                pool_size,
            ),
            (false, false) => run_hkty::<Float, false, false>(
                generators,
                grading_vector,
                zero_cutoff,
                max_deg,
                min_points,
                q,
                nefpart,
                intnums,
                n_threads,
                pool_size,
            ),
        };
        final_res = res
            .into_iter()
            .map(|((v, c), gvgw)| {
                (
                    (to_vec(v), c),
                    if find_gv {
                        gvgw.to_integer().unwrap().to_string()
                    } else {
                        gvgw.to_string()
                    },
                )
            })
            .collect();
    } else {
        let zero_cutoff = Rational::new();
        let res = match (find_gv, is_threefold) {
            (true, true) => run_hkty::<Rational, true, true>(
                generators,
                grading_vector,
                zero_cutoff,
                max_deg,
                min_points,
                q,
                nefpart,
                intnums,
                n_threads,
                pool_size,
            ),
            (false, true) => run_hkty::<Rational, false, true>(
                generators,
                grading_vector,
                zero_cutoff,
                max_deg,
                min_points,
                q,
                nefpart,
                intnums,
                n_threads,
                pool_size,
            ),
            (true, false) => run_hkty::<Rational, true, false>(
                generators,
                grading_vector,
                zero_cutoff,
                max_deg,
                min_points,
                q,
                nefpart,
                intnums,
                n_threads,
                pool_size,
            ),
            (false, false) => run_hkty::<Rational, false, false>(
                generators,
                grading_vector,
                zero_cutoff,
                max_deg,
                min_points,
                q,
                nefpart,
                intnums,
                n_threads,
                pool_size,
            ),
        };
        final_res = res
            .into_iter()
            .map(|((v, c), gvgw)| {
                (
                    (to_vec(v), c),
                    if find_gv {
                        gvgw.into_numer_denom().0.to_string()
                    } else {
                        gvgw.to_string()
                    },
                )
            })
            .collect();
    };

    Ok(final_res)
}

/// A Python module implemented in Rust.
#[pymodule]
pub fn cygv(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compute_gvgw, m)?)?;
    Ok(())
}
