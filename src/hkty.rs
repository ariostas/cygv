//! HKTY procedure.

use crate::polynomial::coefficient::PolynomialCoeff;
use crate::{
    fundamental_period, instanton, misc, series_inversion, PolynomialProperties, Semigroup,
};
use nalgebra::{DMatrix, DVector, RowDVector};
use rug::{ops::PowAssign, Float, Integer, Rational};
use std::collections::HashMap;

// TODO: Need to handle errors properly
#[allow(clippy::too_many_arguments)]
pub fn run_hkty<T, const FIND_GV: bool, const IS_THREEFOLD: bool>(
    generators: DMatrix<i32>,
    grading_vector: RowDVector<i32>,
    zero_cutoff: T,
    max_deg: Option<u32>,
    min_points: Option<u32>,
    q: DMatrix<i32>,
    nefpart: Vec<DVector<i32>>,
    intnums: HashMap<(usize, usize, usize), i32>,
) -> Vec<((DVector<i32>, usize), T)>
where
    T: PolynomialCoeff<T>,
{
    let sg = if let Some(d) = max_deg {
        Semigroup::with_max_degree(generators, grading_vector, d).unwrap()
    } else if let Some(n) = min_points {
        Semigroup::with_min_elements(generators, grading_vector, n as usize).unwrap()
    } else {
        Semigroup::from_data(generators, grading_vector).unwrap()
    };

    let poly_props = PolynomialProperties::new(&sg, &zero_cutoff);

    let (intnum_dict, intnum_idxpairs, n_indices) =
        misc::process_int_nums(intnums, IS_THREEFOLD).unwrap();

    let fp = fundamental_period::compute_omega(&poly_props, &sg, &q, &nefpart, &intnum_idxpairs)
        .unwrap();

    let inst_data = instanton::compute_instanton_data(
        fp,
        &poly_props,
        &intnum_idxpairs,
        n_indices,
        &intnum_dict,
        IS_THREEFOLD,
    )
    .unwrap();

    let gv = series_inversion::invert_series::<T, FIND_GV, IS_THREEFOLD>(inst_data, &poly_props)
        .unwrap();

    let mut gv_sorted: Vec<_> = gv.into_iter().collect();
    gv_sorted.sort_unstable_by_key(|c| c.0 .0);
    gv_sorted
        .into_iter()
        .map(|(k, gv)| {
            (
                (
                    DVector::from_column_slice(
                        poly_props.semigroup.elements.column(k.0).as_slice(),
                    ),
                    k.1,
                ),
                gv,
            )
        })
        .collect()
}

pub fn compute_gv_rat_threefold(
    generators: DMatrix<i32>,
    grading_vector: RowDVector<i32>,
    max_deg: Option<u32>,
    min_points: Option<u32>,
    q: DMatrix<i32>,
    nefpart: Vec<DVector<i32>>,
    intnums: HashMap<(usize, usize, usize), i32>,
) -> Vec<(DVector<i32>, Integer)> {
    let zero_cutoff = Rational::new();
    run_hkty::<Rational, true, true>(
        generators,
        grading_vector,
        zero_cutoff,
        max_deg,
        min_points,
        q,
        nefpart,
        intnums,
    )
    .into_iter()
    .map(|((v, _), gv)| (v, gv.into_numer_denom().0))
    .collect()
}

#[allow(clippy::too_many_arguments)]
pub fn compute_gv_float_threefold(
    generators: DMatrix<i32>,
    grading_vector: RowDVector<i32>,
    max_deg: Option<u32>,
    min_points: Option<u32>,
    q: DMatrix<i32>,
    nefpart: Vec<DVector<i32>>,
    intnums: HashMap<(usize, usize, usize), i32>,
    precision: u32,
) -> Vec<(DVector<i32>, Integer)> {
    let mut zero_cutoff = Float::with_val(precision, 10);
    zero_cutoff.pow_assign(-(precision as i32) / 3);
    run_hkty::<Float, true, true>(
        generators,
        grading_vector,
        zero_cutoff,
        max_deg,
        min_points,
        q,
        nefpart,
        intnums,
    )
    .into_iter()
    .map(|((v, _), gv)| (v, gv.to_integer().unwrap()))
    .collect()
}

pub fn compute_gw_rat_threefold(
    generators: DMatrix<i32>,
    grading_vector: RowDVector<i32>,
    max_deg: Option<u32>,
    min_points: Option<u32>,
    q: DMatrix<i32>,
    nefpart: Vec<DVector<i32>>,
    intnums: HashMap<(usize, usize, usize), i32>,
) -> Vec<(DVector<i32>, Rational)> {
    let zero_cutoff = Rational::new();
    run_hkty::<Rational, false, true>(
        generators,
        grading_vector,
        zero_cutoff,
        max_deg,
        min_points,
        q,
        nefpart,
        intnums,
    )
    .into_iter()
    .map(|((v, _), gv)| (v, gv))
    .collect()
}

#[allow(clippy::too_many_arguments)]
pub fn compute_gw_float_threefold(
    generators: DMatrix<i32>,
    grading_vector: RowDVector<i32>,
    max_deg: Option<u32>,
    min_points: Option<u32>,
    q: DMatrix<i32>,
    nefpart: Vec<DVector<i32>>,
    intnums: HashMap<(usize, usize, usize), i32>,
    precision: u32,
) -> Vec<(DVector<i32>, Float)> {
    let mut zero_cutoff = Float::with_val(precision, 10);
    zero_cutoff.pow_assign(-(precision as i32) / 3);
    run_hkty::<Float, false, true>(
        generators,
        grading_vector,
        zero_cutoff,
        max_deg,
        min_points,
        q,
        nefpart,
        intnums,
    )
    .into_iter()
    .map(|((v, _), gv)| (v, gv))
    .collect()
}

pub fn compute_gv_rat_nfold(
    generators: DMatrix<i32>,
    grading_vector: RowDVector<i32>,
    max_deg: Option<u32>,
    min_points: Option<u32>,
    q: DMatrix<i32>,
    nefpart: Vec<DVector<i32>>,
    intnums: HashMap<(usize, usize, usize), i32>,
) -> Vec<((DVector<i32>, usize), Integer)> {
    let zero_cutoff = Rational::new();
    run_hkty::<Rational, true, false>(
        generators,
        grading_vector,
        zero_cutoff,
        max_deg,
        min_points,
        q,
        nefpart,
        intnums,
    )
    .into_iter()
    .map(|((v, c), gv)| ((v, c), gv.into_numer_denom().0))
    .collect()
}

#[allow(clippy::too_many_arguments)]
pub fn compute_gv_float_nfold(
    generators: DMatrix<i32>,
    grading_vector: RowDVector<i32>,
    max_deg: Option<u32>,
    min_points: Option<u32>,
    q: DMatrix<i32>,
    nefpart: Vec<DVector<i32>>,
    intnums: HashMap<(usize, usize, usize), i32>,
    precision: u32,
) -> Vec<((DVector<i32>, usize), Integer)> {
    let mut zero_cutoff = Float::with_val(precision, 10);
    zero_cutoff.pow_assign(-(precision as i32) / 3);
    run_hkty::<Float, true, false>(
        generators,
        grading_vector,
        zero_cutoff,
        max_deg,
        min_points,
        q,
        nefpart,
        intnums,
    )
    .into_iter()
    .map(|((v, c), gv)| ((v, c), gv.to_integer().unwrap()))
    .collect()
}

pub fn compute_gw_rat_nfold(
    generators: DMatrix<i32>,
    grading_vector: RowDVector<i32>,
    max_deg: Option<u32>,
    min_points: Option<u32>,
    q: DMatrix<i32>,
    nefpart: Vec<DVector<i32>>,
    intnums: HashMap<(usize, usize, usize), i32>,
) -> Vec<((DVector<i32>, usize), Rational)> {
    let zero_cutoff = Rational::new();
    run_hkty::<Rational, false, false>(
        generators,
        grading_vector,
        zero_cutoff,
        max_deg,
        min_points,
        q,
        nefpart,
        intnums,
    )
}

#[allow(clippy::too_many_arguments)]
pub fn compute_gw_float_nfold(
    generators: DMatrix<i32>,
    grading_vector: RowDVector<i32>,
    max_deg: Option<u32>,
    min_points: Option<u32>,
    q: DMatrix<i32>,
    nefpart: Vec<DVector<i32>>,
    intnums: HashMap<(usize, usize, usize), i32>,
    precision: u32,
) -> Vec<((DVector<i32>, usize), Float)> {
    let mut zero_cutoff = Float::with_val(precision, 10);
    zero_cutoff.pow_assign(-(precision as i32) / 3);
    run_hkty::<Float, false, false>(
        generators,
        grading_vector,
        zero_cutoff,
        max_deg,
        min_points,
        q,
        nefpart,
        intnums,
    )
}
