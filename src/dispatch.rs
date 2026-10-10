//! Entry points that pick a pipeline and run it.
//!
//! Every entry point takes a [`Method`]. The modular pipeline supports truncation by a maximum
//! degree and by target points, but not by a minimum number of points, and no vex phases yet;
//! with [`Method::Auto`] the other requests, and those whose curve classes do not fit its keys,
//! go to the reference pipeline instead. Both give the nonzero invariants only.
//!
//! Target points may come with a maximum degree, which then truncates their backward
//! lightcones further; without one, the largest degree of the targets is the maximum.
//!
//! The optional [`ConeData`] lets the modular pipeline take its candidates from the Hilbert
//! bases of the cones $K_T$ instead of enumerating the semigroup of the generators, which is
//! much faster and uses every lattice point of the Mori cone. With it the generators may be
//! empty (no columns), but then only the modular pipeline can run. The reference pipeline
//! ignores the cone data.
//!
//! The reference pipeline still panics on invalid input, while the errors of the modular one
//! are returned.

use crate::cone::{check_cone_data, columns, ConeData};
use crate::modular::key::KeyError;
use crate::modular::{self, CandidateSource, CurveInvariants, Geometry, Options};
use crate::reference::run_hkty;
use crate::{check_strata, CYKind, InvariantKind, Method, Stratum};
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
    /// The strata do not fit the geometry, for the given reason.
    InvalidStrata(&'static str),
    /// The cone data does not fit the geometry or the generators, for the given reason.
    InvalidConeData(&'static str),
    /// The target points do not fit the grading vector, for the given reason.
    InvalidTargetPoints(&'static str),
    /// There are no generators, but the request needs them.
    MissingGenerators,
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
            Self::InvalidStrata(reason)
            | Self::InvalidConeData(reason)
            | Self::InvalidTargetPoints(reason) => write!(f, "{reason}"),
            Self::MissingGenerators => write!(
                f,
                "the generators of the Mori cone are needed, except by the modular method with \
                 cone data"
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
    has_strata: bool,
) -> Option<&'static str> {
    if has_strata {
        Some("strata")
    } else if min_points.is_some() {
        Some("min_points")
    } else if max_deg.is_none() && !has_target_points {
        Some("computations without max_deg or target_points")
    } else {
        None
    }
}

/// Checks that the target points, one per column, have one entry per Kähler parameter and
/// positive degree.
pub(crate) fn check_target_points(
    target_points: &DMatrix<i32>,
    grading_vector: &RowDVector<i32>,
) -> Result<(), &'static str> {
    if target_points.ncols() > 0 && target_points.nrows() != grading_vector.len() {
        return Err("the target points must have one entry per Kähler parameter");
    }
    let positive = target_points.column_iter().all(|p| {
        let deg: i64 = p
            .iter()
            .zip(grading_vector.iter())
            .map(|(&x, &w)| x as i64 * w as i64)
            .sum();
        deg > 0
    });
    if !positive {
        return Err("the target points must have positive degree under the grading vector");
    }
    Ok(())
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
    strata: Vec<Stratum>,
    cone_data: Option<ConeData>,
    invariant_kind: InvariantKind,
    cy_kind: CYKind,
    max_deg: Option<u32>,
    min_points: Option<u32>,
    target_points: Option<DMatrix<i32>>,
    n_threads: Option<u32>,
    method: Method,
) -> Result<Results, ComputeError> {
    check_strata(&strata, q.nrows(), q.ncols(), cy_kind, nefpart.is_empty())
        .map_err(ComputeError::InvalidStrata)?;
    if let Some(cone_data) = &cone_data {
        check_cone_data(cone_data, &q, &generators)
            .map_err(|(_, reason)| ComputeError::InvalidConeData(reason))?;
    }
    // No target points means no truncation by them.
    let target_points = target_points.filter(|t| t.ncols() > 0);
    if let Some(target_points) = &target_points {
        check_target_points(target_points, &grading_vector)
            .map_err(ComputeError::InvalidTargetPoints)?;
    }
    let has_generators = generators.ncols() > 0;
    if !has_generators && cone_data.is_none() {
        return Err(ComputeError::MissingGenerators);
    }
    let unsupported = modular_unsupported(
        max_deg,
        min_points,
        target_points.is_some(),
        !strata.is_empty(),
    );
    let use_modular = match (method, unsupported) {
        (Method::Reference, _) | (Method::Auto, Some(_)) => false,
        (Method::Modular, Some(what)) => return Err(ComputeError::Unsupported(what)),
        (_, None) => true,
    };
    if use_modular {
        let res = run_modular(
            &generators,
            &grading_vector,
            &q,
            &nefpart,
            &intnums,
            cone_data.as_ref(),
            invariant_kind,
            cy_kind,
            max_deg,
            target_points.as_ref(),
        );
        match res {
            Err(ComputeError::Modular(modular::Error::Key(
                KeyError::TooWide { .. } | KeyError::BoundTooLarge,
            ))) if method == Method::Auto && has_generators => {}
            res => return res,
        }
    }
    if !has_generators {
        return Err(ComputeError::MissingGenerators);
    }
    // The reference truncates by the target points alone, so a maximum degree is applied to
    // its results.
    let truncate = max_deg.filter(|_| target_points.is_some());
    let mut res = run_hkty::<Rational>(
        generators,
        grading_vector.clone(),
        invariant_kind,
        cy_kind,
        max_deg,
        min_points,
        target_points,
        q,
        nefpart,
        intnums,
        strata,
        n_threads,
    );
    if let Some(max_deg) = truncate {
        res.retain(|((c, _), _)| grading_vector.dot(&c.transpose()) <= max_deg as i32);
    }
    Ok(Results::Reference(res))
}

/// Runs the modular pipeline with its default options, with the candidates from the cone data
/// if there is any, and from the generators otherwise. One of `max_deg` and `target_points`
/// must be given.
#[allow(clippy::too_many_arguments)]
fn run_modular(
    generators: &DMatrix<i32>,
    grading_vector: &RowDVector<i32>,
    q: &DMatrix<i32>,
    nefpart: &[DVector<i32>],
    intnums: &HashMap<(usize, usize, usize), i32>,
    cone_data: Option<&ConeData>,
    invariant_kind: InvariantKind,
    cy_kind: CYKind,
    max_deg: Option<u32>,
    target_points: Option<&DMatrix<i32>>,
) -> Result<Results, ComputeError> {
    let grading_vector = grading_vector.as_slice();
    let targets: Vec<&[i32]> = target_points.map_or_else(Vec::new, |t| columns(t).collect());
    // The targets lower it to their largest degree.
    let max_deg = max_deg.unwrap_or(u32::MAX);
    let geometry = Geometry::from_toric(q, nefpart, intnums, cy_kind, grading_vector)?;
    let gens: Vec<&[i32]> = columns(generators).collect();
    let source = match cone_data {
        Some(cone_data) => CandidateSource::Cones(cone_data),
        None => CandidateSource::Generators(&gens),
    };
    let options = Options::default();
    Ok(match invariant_kind {
        InvariantKind::GV => Results::ModularGv(modular::compute_gv(
            source,
            grading_vector,
            max_deg,
            &targets,
            &geometry,
            &options,
        )?),
        InvariantKind::GW => Results::ModularGw(modular::compute_gw(
            source,
            grading_vector,
            max_deg,
            &targets,
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
    strata: Vec<Stratum>,
    cone_data: Option<ConeData>,
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
        strata,
        cone_data,
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
    strata: Vec<Stratum>,
    cone_data: Option<ConeData>,
    n_threads: Option<u32>,
    method: Method,
) -> Result<Vec<(DVector<i32>, Integer)>, ComputeError> {
    let res = compute(
        generators,
        grading_vector,
        q,
        nefpart,
        intnums,
        strata,
        cone_data,
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
    strata: Vec<Stratum>,
    cone_data: Option<ConeData>,
    n_threads: Option<u32>,
    method: Method,
) -> Result<Vec<(DVector<i32>, Rational)>, ComputeError> {
    let res = compute(
        generators,
        grading_vector,
        q,
        nefpart,
        intnums,
        strata,
        cone_data,
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
    strata: Vec<Stratum>,
    cone_data: Option<ConeData>,
    n_threads: Option<u32>,
    method: Method,
) -> Result<Vec<Labeled<Integer>>, ComputeError> {
    let res = compute(
        generators,
        grading_vector,
        q,
        nefpart,
        intnums,
        strata,
        cone_data,
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
    strata: Vec<Stratum>,
    cone_data: Option<ConeData>,
    n_threads: Option<u32>,
    method: Method,
) -> Result<Vec<Labeled<Rational>>, ComputeError> {
    let res = compute(
        generators,
        grading_vector,
        q,
        nefpart,
        intnums,
        strata,
        cone_data,
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
