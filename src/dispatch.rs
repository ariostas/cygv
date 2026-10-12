//! Entry points that pick a pipeline and run it.
//!
//! Every entry point takes a [`Method`]. The modular pipeline supports truncation by a maximum
//! degree and by target points, but not by a minimum number of points; with [`Method::Auto`]
//! those requests, and those whose curve classes do not fit its keys, go to the reference
//! pipeline instead. Both support the strata of vex phases. Both give the nonzero invariants only.
//!
//! Target points may come with a maximum degree, which then truncates their backward
//! lightcones further; without one, the largest degree of the targets is the maximum.
//!
//! Both pipelines use every lattice point of the Mori cone, even when the generators do not
//! generate them all (decision 6.2 of `PLAN.md`). The Mori cone is described by the
//! inequalities of the optional [`ConeData`], or else by the facets of the cone the generators
//! span ([`crate::cone::facets`]). The reference pipeline takes the Hilbert basis of its
//! lattice points as its generators ([`crate::cone::hilbert_basis`]), and the modular one
//! takes its candidates from the Hilbert bases of the cones $K_T$: those of the cone data, or
//! else computed ([`crate::cone::cone_data`]). With cone data the generators may be empty (no
//! columns).
//!
//! The optional fan of the ambient toric variety (its maximal cones) gives the strata, which
//! then cannot be given too, and defaults for the Mori cone, spanned by the classes of its
//! walls, used when the generators are empty and there is no cone data, and for the grading
//! vector, used when it is empty ([`crate::phase::Fan`]). The Mori cone is only computed for
//! these defaults, since it takes a double description in dimension $h^{1,1}$. The
//! intersection numbers and a grading vector that is given are checked against the fan.
//!
//! With `certify`, the results of the modular pipeline are proven exact (see
//! [`modular::Options::certify`]); the reference pipeline is exact anyway. Vex phases cannot
//! be certified yet.
//!
//! The input is checked up front (the nef partition, the dimension of the CY, the intersection
//! numbers, the strata, the cone data, the target points, the fan, and the generators of a
//! computation without truncation, which get zero if they miss it), so that invalid input is
//! an error with either pipeline, although [`run_hkty`] itself still panics on it.

use crate::cone::{check_cone_data, columns, ConeData, ConeError};
use crate::modular::key::KeyError;
use crate::modular::{self, CandidateSource, CurveInvariants, Geometry, Options};
use crate::phase::{Fan, FanError};
use crate::reference::{self, run_hkty};
use crate::{check_intnums, check_nefpart, check_strata, CYKind, InvariantKind, Method, Stratum};
use core::fmt;
use malachite::{Integer, Rational};
use nalgebra::{DMatrix, DVector, RowDVector};
use std::collections::HashMap;

/// An invariant, labeled by its curve class and the index of its reference surface, which is
/// zero for a threefold.
pub type Labeled<T> = ((DVector<i32>, usize), T);

/// Errors of the entry points.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ComputeError {
    /// The modular method was requested for something it does not support, which is named.
    Unsupported(&'static str),
    /// The nef partition does not split the rows of the charges, for the given reason.
    InvalidNefPartition(&'static str),
    /// The CY would have the given dimension, which is below three.
    InvalidDimension(i64),
    /// The intersection numbers do not fit the geometry, for the given reason.
    InvalidIntersectionNumbers(&'static str),
    /// The generators cannot be the curve classes of a computation without truncation, for
    /// the given reason.
    InvalidGenerators(&'static str),
    /// The strata do not fit the geometry, for the given reason.
    InvalidStrata(&'static str),
    /// The cone data does not fit the geometry or the generators, for the given reason.
    InvalidConeData(&'static str),
    /// The target points do not fit the grading vector, for the given reason.
    InvalidTargetPoints(&'static str),
    /// There are no generators, but the request needs them: there is no cone data or fan
    /// either, or nothing truncates the computation, which then uses the generators as the
    /// curve classes.
    MissingGenerators,
    /// There is no grading vector (it has no entries), and no fan to derive one from.
    MissingGradingVector,
    /// The fan, or the intersection numbers that come with it, were rejected.
    InvalidFan(FanError),
    /// The Hilbert bases of the Mori cone could not be computed.
    MoriCone(ConeError),
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
            Self::InvalidNefPartition(reason)
            | Self::InvalidIntersectionNumbers(reason)
            | Self::InvalidGenerators(reason)
            | Self::InvalidStrata(reason)
            | Self::InvalidConeData(reason)
            | Self::InvalidTargetPoints(reason) => write!(f, "{reason}"),
            Self::InvalidDimension(dim) => write!(
                f,
                "the dimension of the CY must be at least three, but it is {dim}"
            ),
            Self::MissingGenerators => write!(
                f,
                "the generators are needed, without cone data or a fan, or without max_deg, \
                 min_points and target_points"
            ),
            Self::MissingGradingVector => write!(f, "the grading vector is needed without a fan"),
            Self::InvalidFan(e) => write!(f, "{e}"),
            Self::MoriCone(e) => write!(f, "{e}"),
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
    Reference(Vec<Labeled<reference::Rational>>),
    ModularGv(Vec<CurveInvariants>),
    ModularGw(Vec<CurveInvariants<Rational>>),
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
    fan: Option<Vec<Vec<usize>>>,
    cone_data: Option<ConeData>,
    invariant_kind: InvariantKind,
    cy_kind: CYKind,
    max_deg: Option<u32>,
    min_points: Option<u32>,
    target_points: Option<DMatrix<i32>>,
    n_threads: Option<u32>,
    method: Method,
    certify: bool,
) -> Result<Results, ComputeError> {
    check_nefpart(&nefpart, q.nrows()).map_err(ComputeError::InvalidNefPartition)?;
    let dim = q.nrows() as i64 - q.ncols() as i64 - nefpart.len().max(1) as i64;
    if dim < 3 {
        return Err(ComputeError::InvalidDimension(dim));
    }
    check_intnums(&intnums, q.ncols(), cy_kind)
        .map_err(ComputeError::InvalidIntersectionNumbers)?;
    // From the fan: the strata, and, when needed, the Mori cone, which is expensive at large
    // h11: for a default grading vector, or as the cone without cone data or generators.
    let (strata, grading_vector, mori_cone) = match fan {
        None => (strata, grading_vector, None),
        Some(_) if !cy_kind.is_threefold() || !nefpart.is_empty() => {
            return Err(ComputeError::InvalidFan(FanError::NotHypersurface))
        }
        Some(_) if !strata.is_empty() => {
            return Err(ComputeError::InvalidStrata(
                "the strata are derived from the fan, so they cannot be given with it",
            ))
        }
        Some(cones) => {
            let fan = Fan::new(&cones, &q).map_err(ComputeError::InvalidFan)?;
            fan.check_intnums(&intnums)
                .map_err(ComputeError::InvalidFan)?;
            let strata = fan.strata().map_err(ComputeError::InvalidFan)?;
            let needed =
                grading_vector.is_empty() || (cone_data.is_none() && generators.ncols() == 0);
            let mori_cone = if needed {
                Some(fan.mori_cone().map_err(ComputeError::InvalidFan)?)
            } else {
                None
            };
            let grading_vector = match &mori_cone {
                Some(m) if grading_vector.is_empty() => m.grading_vector.clone(),
                _ => {
                    fan.check_grading(&grading_vector)
                        .map_err(ComputeError::InvalidFan)?;
                    grading_vector
                }
            };
            (strata, grading_vector, mori_cone)
        }
    };
    if grading_vector.is_empty() {
        return Err(ComputeError::MissingGradingVector);
    }
    check_strata(&strata, q.nrows(), q.ncols(), cy_kind, nefpart.is_empty())
        .map_err(ComputeError::InvalidStrata)?;
    if let Some(cone_data) = &cone_data {
        check_cone_data(cone_data, &q, &generators, &strata)
            .map_err(|(_, reason)| ComputeError::InvalidConeData(reason))?;
    }
    // No target points means no truncation by them.
    let target_points = target_points.filter(|t| t.ncols() > 0);
    if let Some(target_points) = &target_points {
        check_target_points(target_points, &grading_vector)
            .map_err(ComputeError::InvalidTargetPoints)?;
    }
    let h11 = grading_vector.len();
    if generators.ncols() > 0 && generators.nrows() != h11 {
        return Err(modular::Error::DimensionMismatch.into());
    }
    // The Mori cone by inequalities: the cone data's, or the fan's without generators, or the
    // facets of the generators.
    let (hyperplanes, equations) = match (&cone_data, mori_cone) {
        (Some(cone_data), _) => (cone_data.hyperplanes.clone(), cone_data.equations.clone()),
        (None, Some(m)) if generators.ncols() == 0 => {
            (m.inequalities.hyperplanes, m.inequalities.equations)
        }
        (None, _) if generators.ncols() == 0 => return Err(ComputeError::MissingGenerators),
        (None, _) => {
            let f = crate::cone::facets(&generators)
                .map_err(|_| ComputeError::MoriCone(ConeError::TooLarge))?;
            (f.hyperplanes, f.equations)
        }
    };
    let unsupported = modular_unsupported(max_deg, min_points, target_points.is_some());
    let use_modular = match (method, unsupported) {
        (Method::Reference, _) | (Method::Auto, Some(_)) => false,
        (Method::Modular, Some(what)) => return Err(ComputeError::Unsupported(what)),
        (_, None) => true,
    };
    if use_modular {
        let computed;
        let cone_data = match &cone_data {
            Some(cone_data) => cone_data,
            None => {
                // The candidates up to the maximum degree, which the targets lower to theirs.
                let targets: Vec<&[i32]> = target_points
                    .as_ref()
                    .map_or_else(Vec::new, |t| columns(t).collect());
                let top = match max_degree(&targets, grading_vector.as_slice()) {
                    Some(d) => max_deg.map_or(d, |m| m.min(d)),
                    // Without targets the modular pipeline has a maximum degree.
                    None => max_deg.unwrap_or(u32::MAX),
                };
                computed = crate::cone::cone_data(
                    &hyperplanes,
                    &equations,
                    &q,
                    &strata,
                    grading_vector.as_slice(),
                    top,
                )
                .map_err(ComputeError::MoriCone)?;
                &computed
            }
        };
        let res = run_modular(
            &grading_vector,
            &q,
            &nefpart,
            &intnums,
            &strata,
            cone_data,
            invariant_kind,
            cy_kind,
            max_deg,
            target_points.as_ref(),
            n_threads,
            certify,
        );
        match res {
            Err(ComputeError::Modular(modular::Error::Key(
                KeyError::TooWide { .. } | KeyError::BoundTooLarge,
            ))) if method == Method::Auto => {}
            res => return res,
        }
    }
    let generators = if max_deg.is_none() && min_points.is_none() && target_points.is_none() {
        // Without any truncation the generators are the curve classes to use, as given.
        if generators.ncols() == 0 {
            return Err(ComputeError::MissingGenerators);
        }
        untruncated_classes(generators, &grading_vector)?
    } else {
        // Every lattice point of the Mori cone, from the Hilbert basis.
        let basis =
            crate::cone::hilbert_basis(&hyperplanes, &equations).map_err(ComputeError::MoriCone)?;
        if basis
            .column_iter()
            .any(|g| grading_vector.dot(&g.transpose()) <= 0)
        {
            return Err(ComputeError::MoriCone(ConeError::NonpositiveDegree));
        }
        basis
    };
    // The reference truncates by the target points alone, so a maximum degree is applied to
    // its results.
    let truncate = max_deg.filter(|_| target_points.is_some());
    let mut res = run_hkty::<reference::Rational>(
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

/// The curve classes of a computation without truncation: the generators, once each, after
/// zero, which is added if it is missing. The others must have positive degree.
fn untruncated_classes(
    generators: DMatrix<i32>,
    grading_vector: &RowDVector<i32>,
) -> Result<DMatrix<i32>, ComputeError> {
    let h11 = generators.nrows();
    let mut seen = std::collections::HashSet::new();
    let mut classes = vec![0; h11];
    for g in generators.as_slice().chunks(h11) {
        if g.iter().all(|&x| x == 0) || !seen.insert(g) {
            continue;
        }
        let deg: i64 = g
            .iter()
            .zip(grading_vector.iter())
            .map(|(&x, &w)| x as i64 * w as i64)
            .sum();
        if deg <= 0 {
            return Err(ComputeError::InvalidGenerators(
                "the nonzero generators must have positive degree under the grading vector",
            ));
        }
        classes.extend_from_slice(g);
    }
    Ok(DMatrix::from_vec(h11, classes.len() / h11, classes))
}

/// The largest degree of the targets, if there are any.
fn max_degree(targets: &[&[i32]], grading_vector: &[i32]) -> Option<u32> {
    targets
        .iter()
        .map(|p| {
            let d: i64 = p
                .iter()
                .zip(grading_vector)
                .map(|(&x, &w)| x as i64 * w as i64)
                .sum();
            d.clamp(0, u32::MAX as i64) as u32
        })
        .max()
}

/// Runs the modular pipeline with its default options, `n_threads` threads and certified if
/// `certify`, with the candidates from the cone data. One of `max_deg` and `target_points`
/// must be given.
#[allow(clippy::too_many_arguments)]
fn run_modular(
    grading_vector: &RowDVector<i32>,
    q: &DMatrix<i32>,
    nefpart: &[DVector<i32>],
    intnums: &HashMap<(usize, usize, usize), i32>,
    strata: &[Stratum],
    cone_data: &ConeData,
    invariant_kind: InvariantKind,
    cy_kind: CYKind,
    max_deg: Option<u32>,
    target_points: Option<&DMatrix<i32>>,
    n_threads: Option<u32>,
    certify: bool,
) -> Result<Results, ComputeError> {
    let grading_vector = grading_vector.as_slice();
    let targets: Vec<&[i32]> = target_points.map_or_else(Vec::new, |t| columns(t).collect());
    // The targets lower it to their largest degree.
    let max_deg = max_deg.unwrap_or(u32::MAX);
    let geometry = Geometry::from_toric(q, nefpart, intnums, cy_kind, grading_vector)?
        .with_strata(strata, grading_vector)?;
    let source = CandidateSource::Cones(cone_data);
    let options = Options {
        // Zero threads means one, as for the reference.
        threads: n_threads.map(|n| n.max(1) as usize),
        certify,
        ..Options::default()
    };
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

impl Results {
    /// The GV invariants, which are integers.
    fn into_gv(self) -> Vec<Labeled<Integer>> {
        match self {
            Self::Reference(r) => r.into_iter().map(|(k, gv)| (k, gv.numerator())).collect(),
            Self::ModularGv(r) => entries(r, &Integer::default()),
            Self::ModularGw(_) => unreachable!("GW invariants were computed instead of GV"),
        }
    }

    /// The GW invariants, which are rational.
    fn into_gw(self) -> Vec<Labeled<Rational>> {
        match self {
            Self::Reference(r) => r.into_iter().map(|(k, gw)| (k, gw.into_inner())).collect(),
            Self::ModularGw(r) => entries(r, &Rational::default()),
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
    fan: Option<Vec<Vec<usize>>>,
    cone_data: Option<ConeData>,
    invariant_kind: InvariantKind,
    cy_kind: CYKind,
    max_deg: Option<u32>,
    min_points: Option<u32>,
    target_points: Option<DMatrix<i32>>,
    n_threads: Option<u32>,
    method: Method,
    certify: bool,
) -> Result<Vec<Labeled<String>>, ComputeError> {
    let res = compute(
        generators,
        grading_vector,
        q,
        nefpart,
        intnums,
        strata,
        fan,
        cone_data,
        invariant_kind,
        cy_kind,
        max_deg,
        min_points,
        target_points,
        n_threads,
        method,
        certify,
    )?;
    Ok(match res {
        Results::Reference(r) => r
            .into_iter()
            .map(|(k, gvgw)| {
                (
                    k,
                    match invariant_kind {
                        InvariantKind::GV => gvgw.numerator().to_string(),
                        InvariantKind::GW => gvgw.to_string(),
                    },
                )
            })
            .collect(),
        Results::ModularGv(r) => entries(r, &Integer::default())
            .into_iter()
            .map(|(k, gv)| (k, gv.to_string()))
            .collect(),
        Results::ModularGw(r) => entries(r, &Rational::default())
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
    fan: Option<Vec<Vec<usize>>>,
    cone_data: Option<ConeData>,
    n_threads: Option<u32>,
    method: Method,
    certify: bool,
) -> Result<Vec<(DVector<i32>, Integer)>, ComputeError> {
    let res = compute(
        generators,
        grading_vector,
        q,
        nefpart,
        intnums,
        strata,
        fan,
        cone_data,
        InvariantKind::GV,
        CYKind::Threefold,
        max_deg,
        min_points,
        target_points,
        n_threads,
        method,
        certify,
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
    fan: Option<Vec<Vec<usize>>>,
    cone_data: Option<ConeData>,
    n_threads: Option<u32>,
    method: Method,
    certify: bool,
) -> Result<Vec<(DVector<i32>, Rational)>, ComputeError> {
    let res = compute(
        generators,
        grading_vector,
        q,
        nefpart,
        intnums,
        strata,
        fan,
        cone_data,
        InvariantKind::GW,
        CYKind::Threefold,
        max_deg,
        min_points,
        target_points,
        n_threads,
        method,
        certify,
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
    fan: Option<Vec<Vec<usize>>>,
    cone_data: Option<ConeData>,
    n_threads: Option<u32>,
    method: Method,
    certify: bool,
) -> Result<Vec<Labeled<Integer>>, ComputeError> {
    let res = compute(
        generators,
        grading_vector,
        q,
        nefpart,
        intnums,
        strata,
        fan,
        cone_data,
        InvariantKind::GV,
        CYKind::Nfold,
        max_deg,
        min_points,
        target_points,
        n_threads,
        method,
        certify,
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
    fan: Option<Vec<Vec<usize>>>,
    cone_data: Option<ConeData>,
    n_threads: Option<u32>,
    method: Method,
    certify: bool,
) -> Result<Vec<Labeled<Rational>>, ComputeError> {
    let res = compute(
        generators,
        grading_vector,
        q,
        nefpart,
        intnums,
        strata,
        fan,
        cone_data,
        InvariantKind::GW,
        CYKind::Nfold,
        max_deg,
        min_points,
        target_points,
        n_threads,
        method,
        certify,
    )?;
    Ok(res.into_gw())
}
