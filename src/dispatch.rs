//! Entry points that pick a pipeline and run it.
//!
//! Every entry point takes a [`Method`]. The modular pipeline only supports truncation by a
//! maximum degree; with [`Method::Auto`] the other requests, and those whose curve classes do
//! not fit its keys, go to the reference pipeline instead. Both give the nonzero invariants
//! only.
//!
//! The reference pipeline still panics on invalid input, while the errors of the modular one
//! are returned.

use crate::modular::key::KeyError;
use crate::modular::{self, CurveInvariants, Geometry, Options};
use crate::reference::run_hkty;
use crate::{CYKind, InvariantKind, Method};
use core::fmt;
use nalgebra::{DMatrix, DVector, RowDVector};
use num_bigint::{BigInt, Sign};
use num_rational::BigRational;
use rug::integer::Order;
use rug::{Integer, Rational};
use std::collections::HashMap;

/// An invariant, labeled by its curve class and the index of its reference surface, which is
/// zero for a threefold.
pub type Labeled<T> = ((DVector<i32>, usize), T);

/// Errors of the entry points.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ComputeError {
    /// The modular method was requested for something it does not support, which is named.
    Unsupported(&'static str),
    /// The modular pipeline failed.
    Modular(modular::Error),
}

impl fmt::Display for ComputeError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::Unsupported(what) => write!(
                f,
                "the modular method does not support {what}; use the auto or reference method"
            ),
            Self::Modular(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for ComputeError {}

impl From<modular::Error> for ComputeError {
    fn from(e: modular::Error) -> Self {
        Self::Modular(e)
    }
}

/// What makes a request unsupported by the modular pipeline, if anything.
pub(crate) fn modular_unsupported(
    max_deg: Option<u32>,
    min_points: Option<u32>,
    has_target_points: bool,
) -> Option<&'static str> {
    if min_points.is_some() {
        Some("min_points")
    } else if has_target_points {
        Some("target_points")
    } else if max_deg.is_none() {
        Some("computations without max_deg")
    } else {
        None
    }
}

/// The invariants as one of the pipelines gives them.
enum Results {
    Reference(Vec<Labeled<Rational>>),
    ModularGv(Vec<CurveInvariants>),
    ModularGw(Vec<CurveInvariants<BigRational>>),
}

/// Runs the pipeline that `method` selects.
#[allow(clippy::too_many_arguments)]
fn compute(
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
    method: Method,
) -> Result<Results, ComputeError> {
    let unsupported = modular_unsupported(max_deg, min_points, target_points.is_some());
    let modular_max_deg = match (method, unsupported) {
        (Method::Reference, _) | (Method::Auto, Some(_)) => None,
        (Method::Modular, Some(what)) => return Err(ComputeError::Unsupported(what)),
        (_, None) => max_deg,
    };
    if let Some(max_deg) = modular_max_deg {
        let res = run_modular(
            &generators,
            &grading_vector,
            &q,
            &nefpart,
            &intnums,
            invariant_kind,
            cy_kind,
            max_deg,
        );
        match res {
            Err(ComputeError::Modular(modular::Error::Key(
                KeyError::TooWide { .. } | KeyError::BoundTooLarge,
            ))) if method == Method::Auto => {}
            res => return res,
        }
    }
    Ok(Results::Reference(run_hkty::<Rational>(
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
    )))
}

/// Runs the modular pipeline with its default options.
#[allow(clippy::too_many_arguments)]
fn run_modular(
    generators: &DMatrix<i32>,
    grading_vector: &RowDVector<i32>,
    q: &DMatrix<i32>,
    nefpart: &[DVector<i32>],
    intnums: &HashMap<(usize, usize, usize), i32>,
    invariant_kind: InvariantKind,
    cy_kind: CYKind,
    max_deg: u32,
) -> Result<Results, ComputeError> {
    let grading_vector = grading_vector.as_slice();
    let geometry = Geometry::from_toric(q, nefpart, intnums, cy_kind, grading_vector)?;
    let gens: Vec<&[i32]> = generators
        .as_slice()
        .chunks(generators.nrows().max(1))
        .collect();
    let options = Options::default();
    Ok(match invariant_kind {
        InvariantKind::GV => Results::ModularGv(modular::compute_gv(
            &gens,
            grading_vector,
            max_deg,
            &geometry,
            &options,
        )?),
        InvariantKind::GW => Results::ModularGw(modular::compute_gw(
            &gens,
            grading_vector,
            max_deg,
            &geometry,
            &options,
        )?),
    })
}

/// The nonzero invariants of the modular pipeline, labeled by curve class and component.
fn entries<T: PartialEq>(invariants: Vec<CurveInvariants<T>>, zero: &T) -> Vec<Labeled<T>> {
    let mut res = Vec::new();
    for inv in invariants {
        let curve = DVector::from_vec(inv.curve);
        for (k, v) in inv.values.into_iter().enumerate() {
            if v != *zero {
                res.push(((curve.clone(), k), v));
            }
        }
    }
    res
}

fn to_integer(x: &BigInt) -> Integer {
    let (sign, digits) = x.to_u32_digits();
    let x = Integer::from_digits(&digits, Order::Lsf);
    if sign == Sign::Minus {
        -x
    } else {
        x
    }
}

fn to_rational(x: BigRational) -> Rational {
    let (num, den) = x.into_raw();
    Rational::from((to_integer(&num), to_integer(&den)))
}

impl Results {
    /// The GV invariants, which are integers.
    fn into_gv(self) -> Vec<Labeled<Integer>> {
        match self {
            Self::Reference(r) => r
                .into_iter()
                .map(|(k, gv)| (k, gv.into_numer_denom().0))
                .collect(),
            Self::ModularGv(r) => entries(r, &BigInt::from(0))
                .into_iter()
                .map(|(k, gv)| (k, to_integer(&gv)))
                .collect(),
            Self::ModularGw(_) => unreachable!("GW invariants were computed instead of GV"),
        }
    }

    /// The GW invariants, which are rational.
    fn into_gw(self) -> Vec<Labeled<Rational>> {
        match self {
            Self::Reference(r) => r,
            Self::ModularGw(r) => entries(r, &BigRational::from_integer(0.into()))
                .into_iter()
                .map(|(k, gw)| (k, to_rational(gw)))
                .collect(),
            Self::ModularGv(_) => unreachable!("GV invariants were computed instead of GW"),
        }
    }
}

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
    method: Method,
) -> Result<Vec<Labeled<String>>, ComputeError> {
    let res = compute(
        generators,
        grading_vector,
        q,
        nefpart,
        intnums,
        invariant_kind,
        cy_kind,
        max_deg,
        min_points,
        target_points,
        n_threads,
        method,
    )?;
    Ok(match res {
        Results::Reference(r) => r
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
            .collect(),
        Results::ModularGv(r) => entries(r, &BigInt::from(0))
            .into_iter()
            .map(|(k, gv)| (k, gv.to_string()))
            .collect(),
        Results::ModularGw(r) => entries(r, &BigRational::from_integer(0.into()))
            .into_iter()
            .map(|(k, gw)| (k, gw.to_string()))
            .collect(),
    })
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
    method: Method,
) -> Result<Vec<(DVector<i32>, Integer)>, ComputeError> {
    let res = compute(
        generators,
        grading_vector,
        q,
        nefpart,
        intnums,
        InvariantKind::GV,
        CYKind::Threefold,
        max_deg,
        min_points,
        target_points,
        n_threads,
        method,
    )?;
    Ok(res
        .into_gv()
        .into_iter()
        .map(|((v, _), gv)| (v, gv))
        .collect())
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
    method: Method,
) -> Result<Vec<(DVector<i32>, Rational)>, ComputeError> {
    let res = compute(
        generators,
        grading_vector,
        q,
        nefpart,
        intnums,
        InvariantKind::GW,
        CYKind::Threefold,
        max_deg,
        min_points,
        target_points,
        n_threads,
        method,
    )?;
    Ok(res
        .into_gw()
        .into_iter()
        .map(|((v, _), gw)| (v, gw))
        .collect())
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
    method: Method,
) -> Result<Vec<Labeled<Integer>>, ComputeError> {
    let res = compute(
        generators,
        grading_vector,
        q,
        nefpart,
        intnums,
        InvariantKind::GV,
        CYKind::Nfold,
        max_deg,
        min_points,
        target_points,
        n_threads,
        method,
    )?;
    Ok(res.into_gv())
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
    method: Method,
) -> Result<Vec<Labeled<Rational>>, ComputeError> {
    let res = compute(
        generators,
        grading_vector,
        q,
        nefpart,
        intnums,
        InvariantKind::GW,
        CYKind::Nfold,
        max_deg,
        min_points,
        target_points,
        n_threads,
        method,
    )?;
    Ok(res.into_gw())
}
