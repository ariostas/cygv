//! Entry points that pick a pipeline and run it.

use crate::reference::run_hkty;
use crate::{CYKind, InvariantKind};
use nalgebra::{DMatrix, DVector, RowDVector};
use rug::{Integer, Rational};
use std::collections::HashMap;

/// Compute GV or GW invariants, formatting them as strings.
///
/// This is meant for interfaces that only learn the coefficient type at runtime,
/// such as the Python bindings and the command line interface. GV invariants are
/// formatted as integers, and GW invariants as fractions.
#[allow(clippy::too_many_arguments)]
pub fn compute_gvgw_strings(
    generators: DMatrix<i32>,
    grading_vector: RowDVector<i32>,
    q: DMatrix<i32>,
    nefpart: Vec<DVector<i32>>,
    intnums: HashMap<(usize, usize, usize), i32>,
    invariant_kind: InvariantKind,
    cy_kind: CYKind,
    max_deg: Option<u32>,
    min_points: Option<u32>,
    target_points: Option<DMatrix<i32>>,
    n_threads: Option<u32>,
) -> Vec<((DVector<i32>, usize), String)> {
    run_hkty::<Rational>(
        generators,
        grading_vector,
        invariant_kind,
        cy_kind,
        max_deg,
        min_points,
        target_points,
        q,
        nefpart,
        intnums,
        n_threads,
    )
    .into_iter()
    .map(|(k, gvgw)| {
        (
            k,
            match invariant_kind {
                InvariantKind::GV => gvgw.into_numer_denom().0.to_string(),
                InvariantKind::GW => gvgw.to_string(),
            },
        )
    })
    .collect()
}

#[allow(clippy::too_many_arguments)]
pub fn compute_gv_rat_threefold(
    generators: DMatrix<i32>,
    grading_vector: RowDVector<i32>,
    max_deg: Option<u32>,
    min_points: Option<u32>,
    target_points: Option<DMatrix<i32>>,
    q: DMatrix<i32>,
    nefpart: Vec<DVector<i32>>,
    intnums: HashMap<(usize, usize, usize), i32>,
    n_threads: Option<u32>,
) -> Vec<(DVector<i32>, Integer)> {
    run_hkty::<Rational>(
        generators,
        grading_vector,
        InvariantKind::GV,
        CYKind::Threefold,
        max_deg,
        min_points,
        target_points,
        q,
        nefpart,
        intnums,
        n_threads,
    )
    .into_iter()
    .map(|((v, _), gv)| (v, gv.into_numer_denom().0))
    .collect()
}

#[allow(clippy::too_many_arguments)]
pub fn compute_gw_rat_threefold(
    generators: DMatrix<i32>,
    grading_vector: RowDVector<i32>,
    max_deg: Option<u32>,
    min_points: Option<u32>,
    target_points: Option<DMatrix<i32>>,
    q: DMatrix<i32>,
    nefpart: Vec<DVector<i32>>,
    intnums: HashMap<(usize, usize, usize), i32>,
    n_threads: Option<u32>,
) -> Vec<(DVector<i32>, Rational)> {
    run_hkty::<Rational>(
        generators,
        grading_vector,
        InvariantKind::GW,
        CYKind::Threefold,
        max_deg,
        min_points,
        target_points,
        q,
        nefpart,
        intnums,
        n_threads,
    )
    .into_iter()
    .map(|((v, _), gv)| (v, gv))
    .collect()
}

#[allow(clippy::too_many_arguments)]
pub fn compute_gv_rat_nfold(
    generators: DMatrix<i32>,
    grading_vector: RowDVector<i32>,
    max_deg: Option<u32>,
    min_points: Option<u32>,
    target_points: Option<DMatrix<i32>>,
    q: DMatrix<i32>,
    nefpart: Vec<DVector<i32>>,
    intnums: HashMap<(usize, usize, usize), i32>,
    n_threads: Option<u32>,
) -> Vec<((DVector<i32>, usize), Integer)> {
    run_hkty::<Rational>(
        generators,
        grading_vector,
        InvariantKind::GV,
        CYKind::Nfold,
        max_deg,
        min_points,
        target_points,
        q,
        nefpart,
        intnums,
        n_threads,
    )
    .into_iter()
    .map(|((v, c), gv)| ((v, c), gv.into_numer_denom().0))
    .collect()
}

#[allow(clippy::too_many_arguments)]
pub fn compute_gw_rat_nfold(
    generators: DMatrix<i32>,
    grading_vector: RowDVector<i32>,
    max_deg: Option<u32>,
    min_points: Option<u32>,
    target_points: Option<DMatrix<i32>>,
    q: DMatrix<i32>,
    nefpart: Vec<DVector<i32>>,
    intnums: HashMap<(usize, usize, usize), i32>,
    n_threads: Option<u32>,
) -> Vec<((DVector<i32>, usize), Rational)> {
    run_hkty::<Rational>(
        generators,
        grading_vector,
        InvariantKind::GW,
        CYKind::Nfold,
        max_deg,
        min_points,
        target_points,
        q,
        nefpart,
        intnums,
        n_threads,
    )
}
