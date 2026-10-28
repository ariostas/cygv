//! The entry points: [`compute_gv`], [`compute_gw`] and [`compute`], which check a
//! [`Problem`], resolve its defaults, pick a pipeline and run it.
//!
//! Both pipelines use every lattice point of the Mori cone, even when the generators do not
//! generate them all: the cone is described by its inequalities,
//! given or else the facets of the generators ([`crate::cone::facets`]) or of the fan's wall
//! classes. The reference pipeline takes the Hilbert basis of its lattice points as its
//! generators ([`crate::cone::hilbert_basis`]), and the modular one takes its candidates from
//! the Hilbert bases of the cones $K_T$: those of [`Options::hilbert_bases`], or else computed
//! ([`crate::cone::cone_data`]).
//!
//! [`Method::Auto`] runs the modular pipeline for every selection but
//! [`Selection::MinPoints`], which only the reference one supports. Curve classes
//! that do not fit the modular keys are an error, since the reference pipeline, which takes
//! them, is too slow where they occur (high $h^{1,1}$) to run unasked.

use crate::cone::{columns, ConeData, ConeError};
use crate::modular::key::KeyError;
use crate::modular::{self, CandidateSource, CurveInvariants};
use crate::problem::{
    check_hilbert_bases, check_target_degrees, cone_field, input_err, Checked, MoriCone, Options,
    Selection,
};
use crate::reference::{self, run_hkty};
use crate::{CYKind, InvariantKind, Method, Problem, Stratum};
use core::fmt;
use malachite::{Integer, Rational};
use nalgebra::{DMatrix, DVector, RowDVector};
use std::time::{Duration, Instant};

/// The errors of the entry points.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Error {
    /// Invalid input: the field at fault, named as the Python keyword and the YAML key, and
    /// what is wrong with it.
    Input { field: &'static str, reason: String },
    /// The input is valid, but the method cannot compute it, for the reason given.
    Unsupported(String),
    /// The cone computations failed.
    MoriCone(ConeError),
    /// The modular pipeline failed.
    Modular(modular::Error),
    /// The reference pipeline failed.
    Reference(reference::Error),
    /// The computation was interrupted through [`Options::interrupt`].
    Interrupted,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::Input { field, reason } => write!(f, "invalid {field}: {reason}"),
            Self::Unsupported(reason) => write!(f, "{reason}"),
            Self::MoriCone(e) => write!(f, "{e}"),
            Self::Modular(e) => write!(f, "{e}"),
            Self::Reference(e) => write!(f, "{e}"),
            Self::Interrupted => write!(f, "the computation was interrupted"),
        }
    }
}

impl std::error::Error for Error {}

impl From<modular::Error> for Error {
    fn from(e: modular::Error) -> Self {
        match e {
            modular::Error::Interrupted => Self::Interrupted,
            e => Self::Modular(e),
        }
    }
}

impl From<reference::Error> for Error {
    fn from(e: reference::Error) -> Self {
        match e {
            reference::Error::Interrupted => Self::Interrupted,
            e => Self::Reference(e),
        }
    }
}

impl From<ConeError> for Error {
    fn from(e: ConeError) -> Self {
        match e {
            ConeError::Interrupted => Self::Interrupted,
            e => Self::MoriCone(e),
        }
    }
}

/// The invariants of a [`Problem`], with what is needed to read them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Invariants<T> {
    /// The dimension of the Calabi-Yau.
    pub dimension: usize,
    /// The grading vector the degrees are taken with: the one given, or the default.
    pub grading_vector: Vec<i32>,
    /// The nonzero invariants, sorted by degree, then by curve class, then by surface.
    pub entries: Vec<Invariant<T>>,
    /// How they were computed.
    pub report: Report,
}

/// The invariant of one curve class (and, for an n-fold, one reference surface).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Invariant<T> {
    pub curve: Vec<i32>,
    /// The reference surface of an n-fold, `None` for a threefold.
    pub surface: Option<usize>,
    pub degree: u32,
    pub value: T,
}

/// How the invariants were computed.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct Report {
    /// The pipeline that ran, never [`Method::Auto`].
    pub method: Method,
    /// Whether the results are proven exact: by the reference pipeline, which is exact, or by
    /// the modular one with [`Options::certify`].
    pub certified: bool,
    /// Whether the grading vector is the default rather than given.
    pub grading_vector_derived: bool,
    /// The wall-clock time of the computation.
    pub elapsed: Duration,
}

/// GV or GW invariants, for callers that only know which at run time.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Computed {
    Gv(Invariants<Integer>),
    Gw(Invariants<Rational>),
}

/// Computes the GV invariants of the problem.
pub fn compute_gv(problem: &Problem, options: &Options) -> Result<Invariants<Integer>, Error> {
    match compute(problem, InvariantKind::GV, options)? {
        Computed::Gv(gv) => Ok(gv),
        Computed::Gw(_) => unreachable!("GW invariants were computed instead of GV"),
    }
}

/// Computes the GW invariants of the problem.
pub fn compute_gw(problem: &Problem, options: &Options) -> Result<Invariants<Rational>, Error> {
    match compute(problem, InvariantKind::GW, options)? {
        Computed::Gw(gw) => Ok(gw),
        Computed::Gv(_) => unreachable!("GV invariants were computed instead of GW"),
    }
}

/// Computes the invariants of the given kind.
pub fn compute(
    problem: &Problem,
    kind: InvariantKind,
    options: &Options,
) -> Result<Computed, Error> {
    crate::interrupt::with(options.interrupt.as_ref(), || {
        compute_in(problem, kind, options)
    })
}

/// [`compute`], with the interrupt of the options current.
fn compute_in(
    problem: &Problem,
    kind: InvariantKind,
    options: &Options,
) -> Result<Computed, Error> {
    let start = Instant::now();
    let Checked {
        cy_kind,
        strata,
        fan,
    } = problem.check(options)?;
    let geometry = &problem.geometry;
    let q = &geometry.q;
    let fan_err = |field| {
        move |e: crate::phase::FanError| match e {
            crate::phase::FanError::Interrupted => Error::Interrupted,
            e => input_err(field, e.to_string()),
        }
    };
    // The Mori cone by inequalities.
    let (hyperplanes, equations) = match &problem.cone {
        MoriCone::Generators(g) => {
            let f = crate::cone::facets(g).map_err(ConeError::from)?;
            (f.hyperplanes, f.equations)
        }
        MoriCone::Inequalities {
            hyperplanes,
            equations,
        } => {
            let equations = if equations.ncols() == 0 {
                DMatrix::zeros(q.ncols(), 0)
            } else {
                equations.clone()
            };
            (hyperplanes.clone(), equations)
        }
        MoriCone::FromFan => {
            let fan = fan.as_ref().expect("the fan was checked");
            let m = fan.mori_cone().map_err(fan_err("fan"))?;
            (m.inequalities.hyperplanes, m.inequalities.equations)
        }
    };
    log::info!(
        "Mori cone: {} inequalities and {} equations{}",
        hyperplanes.ncols(),
        equations.ncols(),
        match &problem.cone {
            MoriCone::Generators(_) => ", from the generators",
            MoriCone::Inequalities { .. } => "",
            MoriCone::FromFan => ", from the fan",
        }
    );
    // A Mori cone given with the fan, used as given, must contain the fan's.
    if let (Some(fan), MoriCone::Generators(_) | MoriCone::Inequalities { .. }) =
        (&fan, &problem.cone)
    {
        if let Some(c) = fan
            .wall_outside(&hyperplanes, &equations)
            .map_err(fan_err("fan"))?
        {
            return Err(input_err(
                cone_field(&problem.cone),
                format!(
                    "the Mori cone does not contain the curve class {c:?} of a wall of the \
                     fan; leave it out to derive it from the fan"
                ),
            ));
        }
    }
    let grading_vector = match &problem.grading_vector {
        // Checked with the problem.
        Some(g) => RowDVector::from_row_slice(g),
        None => default_grading(&hyperplanes).ok_or_else(|| {
            input_err(
                cone_field(&problem.cone),
                "the default grading vector, the sum of the inequalities of the Mori cone, \
                 does not fit in 32 bits or vanishes; give a grading vector",
            )
        })?,
    };
    // Curve classes are computed as the targets of their lightcones, which hold every class
    // their invariants depend on, and only theirs are kept.
    let mut listed = None;
    let selection = match &problem.selection {
        Selection::CurveClasses(classes) => {
            let points = curve_classes(classes, &grading_vector)?;
            listed = Some(
                columns(&points)
                    .map(<[i32]>::to_vec)
                    .collect::<std::collections::HashSet<_>>(),
            );
            std::borrow::Cow::Owned(Selection::TargetPoints {
                points,
                max_deg: None,
            })
        }
        selection => std::borrow::Cow::Borrowed(selection),
    };
    if let Selection::TargetPoints { points, .. } = &*selection {
        check_target_degrees(points, grading_vector.as_slice())?;
        let in_cone = |p: &[i32]| {
            columns(&hyperplanes).all(|h| dot(h, p) >= 0)
                && columns(&equations).all(|e| dot(e, p) == 0)
        };
        if let Some(p) = columns(points).find(|p| !in_cone(p)) {
            let field = problem.selection.name();
            return Err(input_err(
                field,
                format!(
                    "the {} {p:?} lies outside the Mori cone",
                    if listed.is_some() {
                        "curve class"
                    } else {
                        "target point"
                    }
                ),
            ));
        }
    }
    let use_modular = options.method != Method::Reference && selection.modular_supported();
    log::info!(
        "grading vector {:?}{}, {} method",
        grading_vector.as_slice(),
        if problem.grading_vector.is_none() {
            " (derived)"
        } else {
            ""
        },
        if use_modular { "modular" } else { "reference" }
    );
    let run = || -> Result<Entries, Error> {
        Ok(if use_modular {
            let generators = match &problem.cone {
                MoriCone::Generators(g) => g.clone(),
                _ => DMatrix::zeros(q.ncols(), 0),
            };
            let computed;
            let cone_data = match &options.hilbert_bases {
                Some(cones) => {
                    computed = check_hilbert_bases(
                        cones,
                        &hyperplanes,
                        &equations,
                        q,
                        Some(&generators),
                        &strata,
                    )?;
                    &computed
                }
                None => {
                    computed = crate::cone::cone_data(
                        &hyperplanes,
                        &equations,
                        q,
                        &strata,
                        grading_vector.as_slice(),
                        top_degree(&selection, &grading_vector),
                    )?;
                    &computed
                }
            };
            log::info!(
                "cone data: {} cones, {} Hilbert basis elements{}",
                cone_data.cones.len(),
                cone_data
                    .cones
                    .iter()
                    .map(|c| c.hilbert_basis.ncols())
                    .sum::<usize>(),
                if options.hilbert_bases.is_some() {
                    ", as given"
                } else {
                    ""
                }
            );
            run_modular(
                problem,
                &selection,
                cy_kind,
                &strata,
                &grading_vector,
                cone_data,
                kind,
                options,
            )
            .map_err(|e| match e {
                Error::Modular(modular::Error::Key(
                    e @ (KeyError::TooWide { .. } | KeyError::BoundTooLarge),
                )) => Error::Unsupported(format!(
                    "the curve classes are too large for the modular method ({e}); the \
                         reference method can compute them, but slowly"
                )),
                e => e,
            })?
        } else {
            // The reference does not use the Hilbert bases, but they are checked all the same.
            if let Some(cones) = &options.hilbert_bases {
                check_hilbert_bases(cones, &hyperplanes, &equations, q, None, &strata)?;
            }
            run_reference(
                problem,
                &selection,
                cy_kind,
                strata,
                &grading_vector,
                &hyperplanes,
                &equations,
                kind,
                options,
            )?
        })
    };
    let mut res = run().map_err(|e| cone_input(problem, e))?;
    if let Some(listed) = &listed {
        match &mut res {
            Entries::Gv(e) => e.retain(|i| listed.contains(&i.curve)),
            Entries::Gw(e) => e.retain(|i| listed.contains(&i.curve)),
        }
    }
    log::info!(
        "done: {} nonzero invariants",
        match &res {
            Entries::Gv(e) => e.len(),
            Entries::Gw(e) => e.len(),
        }
    );
    let report = Report {
        method: if use_modular {
            Method::Modular
        } else {
            Method::Reference
        },
        certified: !use_modular || options.certify,
        grading_vector_derived: problem.grading_vector.is_none(),
        elapsed: start.elapsed(),
    };
    let dimension = geometry.dimension() as usize;
    let grading_vector = grading_vector.as_slice().to_vec();
    Ok(match &mut res {
        Entries::Gv(e) => Computed::Gv(Invariants {
            dimension,
            grading_vector,
            entries: sorted(std::mem::take(e)),
            report,
        }),
        Entries::Gw(e) => Computed::Gw(Invariants {
            dimension,
            grading_vector,
            entries: sorted(std::mem::take(e)),
            report,
        }),
    })
}

/// An error of the cone computations that comes from the input, as an input error of its
/// field: the Mori cone's, or the grading vector's if it is given and not positive.
fn cone_input(problem: &Problem, e: Error) -> Error {
    match e {
        Error::MoriCone(e @ ConeError::NonpositiveDegree) if problem.grading_vector.is_some() => {
            input_err("grading_vector", e.to_string())
        }
        Error::MoriCone(e @ (ConeError::NotPointed | ConeError::NonpositiveDegree)) => {
            input_err(cone_field(&problem.cone), e.to_string())
        }
        e => e,
    }
}

/// The invariants as a pipeline gives them, before they are sorted.
enum Entries {
    Gv(Vec<Invariant<Integer>>),
    Gw(Vec<Invariant<Rational>>),
}

/// The entries sorted by degree, then by curve class, then by surface.
fn sorted<T>(mut entries: Vec<Invariant<T>>) -> Vec<Invariant<T>> {
    entries.sort_unstable_by(|a, b| {
        (a.degree, &a.curve, a.surface).cmp(&(b.degree, &b.curve, b.surface))
    });
    entries
}

/// The degree of a curve class under the grading vector.
fn degree(c: &[i32], grading_vector: &RowDVector<i32>) -> i64 {
    c.iter()
        .zip(grading_vector.iter())
        .map(|(&x, &w)| x as i64 * w as i64)
        .sum()
}

/// The default grading vector: the primitive sum of the hyperplanes, if it fits in `i32` and
/// does not vanish.
fn default_grading(hyperplanes: &DMatrix<i32>) -> Option<RowDVector<i32>> {
    let mut sum = vec![0i64; hyperplanes.nrows()];
    for h in columns(hyperplanes) {
        for (s, &x) in sum.iter_mut().zip(h) {
            *s += x as i64;
        }
    }
    let g = sum.iter().fold(0u64, |g, &x| gcd(g, x.unsigned_abs()));
    if g == 0 {
        return None;
    }
    let entries: Option<Vec<i32>> = sum
        .iter()
        .map(|&x| i32::try_from(x / g as i64).ok())
        .collect();
    Some(RowDVector::from_vec(entries?))
}

fn gcd(a: u64, b: u64) -> u64 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

/// The degree up to which the modular pipeline computes: the maximum degree, lowered to the
/// largest degree of the targets.
fn top_degree(selection: &Selection, grading_vector: &RowDVector<i32>) -> u32 {
    match selection {
        Selection::MaxDeg(d) => *d,
        Selection::TargetPoints { points, max_deg } => {
            let top = columns(points)
                .map(|p| degree(p, grading_vector).clamp(0, u32::MAX as i64) as u32)
                .max()
                .unwrap_or(0);
            max_deg.map_or(top, |m| m.min(top))
        }
        _ => unreachable!("the modular pipeline only runs for max_deg and target_points"),
    }
}

/// The nef partition as the pipelines take it.
fn nefpart_vectors(nefpart: &[Vec<usize>]) -> Vec<DVector<i32>> {
    nefpart
        .iter()
        .map(|part| DVector::from_iterator(part.len(), part.iter().map(|&r| r as i32)))
        .collect()
}

/// Runs the modular pipeline, with the candidates from the cone data.
#[allow(clippy::too_many_arguments)]
fn run_modular(
    problem: &Problem,
    selection: &Selection,
    cy_kind: CYKind,
    strata: &[Stratum],
    grading_vector: &RowDVector<i32>,
    cone_data: &ConeData,
    kind: InvariantKind,
    options: &Options,
) -> Result<Entries, Error> {
    let g = grading_vector.as_slice();
    let (max_deg, targets): (u32, Vec<&[i32]>) = match selection {
        Selection::MaxDeg(d) => (*d, Vec::new()),
        // The targets lower it to their largest degree.
        Selection::TargetPoints { points, max_deg } => {
            (max_deg.unwrap_or(u32::MAX), columns(points).collect())
        }
        _ => unreachable!("the modular pipeline only runs for max_deg and target_points"),
    };
    let geometry = &problem.geometry;
    let geometry = modular::Geometry::from_toric(
        &geometry.q,
        &nefpart_vectors(&geometry.nefpart),
        &geometry.intnums,
        cy_kind,
        g,
    )?
    .with_strata(strata, g)?;
    let source = CandidateSource::Cones(cone_data);
    let modular_options = modular::Options {
        threads: options.n_threads,
        certify: options.certify,
        low_memory: options.low_memory,
        tuning: options.tuning.clone(),
    };
    let threefold = cy_kind.is_threefold();
    Ok(match kind {
        InvariantKind::GV => Entries::Gv(entries(
            modular::compute_gv(source, g, max_deg, &targets, &geometry, &modular_options)?,
            threefold,
        )),
        InvariantKind::GW => Entries::Gw(entries(
            modular::compute_gw(source, g, max_deg, &targets, &geometry, &modular_options)?,
            threefold,
        )),
    })
}

/// The nonzero invariants of the modular pipeline, one per curve class and component.
fn entries<T: PartialEq + Default>(
    invariants: Vec<CurveInvariants<T>>,
    threefold: bool,
) -> Vec<Invariant<T>> {
    let zero = T::default();
    let n = invariants
        .iter()
        .map(|inv| inv.values.iter().filter(|&v| *v != zero).count())
        .sum();
    let mut res = Vec::with_capacity(n);
    // Converted from the end, the input shrinking as it goes, so that the two are never both
    // whole.
    let mut invariants = invariants;
    let step = (invariants.len() / 8).max(1 << 16);
    while let Some(inv) = invariants.pop() {
        if invariants.capacity() - invariants.len() >= step {
            invariants.shrink_to_fit();
        }
        let mut values: Vec<(usize, T)> = inv
            .values
            .into_iter()
            .enumerate()
            .filter(|(_, v)| *v != zero)
            .collect();
        // Backwards, like the curves; the first component takes the curve, the others a copy.
        let first = (!values.is_empty()).then(|| values.remove(0));
        for (k, value) in values.into_iter().rev() {
            res.push(Invariant {
                curve: inv.curve.clone(),
                surface: (!threefold).then_some(k),
                degree: inv.degree,
                value,
            });
        }
        if let Some((k, value)) = first {
            res.push(Invariant {
                curve: inv.curve,
                surface: (!threefold).then_some(k),
                degree: inv.degree,
                value,
            });
        }
    }
    res.reverse();
    res
}

/// Runs the reference pipeline, on the Hilbert basis of the Mori cone.
#[allow(clippy::too_many_arguments)]
fn run_reference(
    problem: &Problem,
    selection: &Selection,
    cy_kind: CYKind,
    strata: Vec<Stratum>,
    grading_vector: &RowDVector<i32>,
    hyperplanes: &DMatrix<i32>,
    equations: &DMatrix<i32>,
    kind: InvariantKind,
    options: &Options,
) -> Result<Entries, Error> {
    // Every lattice point of the Mori cone, from the Hilbert basis.
    let generators = crate::cone::hilbert_basis(hyperplanes, equations)?;
    log::info!(
        "Hilbert basis of the Mori cone: {} elements",
        generators.ncols()
    );
    if columns(&generators).any(|g| degree(g, grading_vector) <= 0) {
        return Err(ConeError::NonpositiveDegree.into());
    }
    let (max_deg, min_points, target_points) = match selection {
        Selection::MaxDeg(d) => (Some(*d), None, None),
        Selection::TargetPoints { points, max_deg } => (*max_deg, None, Some(points.clone())),
        Selection::MinPoints(n) => (None, Some(*n), None),
        Selection::CurveClasses(_) => unreachable!("curve classes run as target points"),
    };
    // The reference truncates by the target points alone, so a maximum degree is applied to
    // its results.
    let truncate = max_deg.filter(|_| target_points.is_some());
    let geometry = &problem.geometry;
    let mut res = run_hkty::<reference::Rational>(
        generators,
        grading_vector.clone(),
        kind,
        cy_kind,
        max_deg,
        min_points,
        target_points,
        geometry.q.clone(),
        nefpart_vectors(&geometry.nefpart),
        geometry.intnums.clone(),
        strata,
        options.n_threads.map(|n| n as u32),
    )?;
    if let Some(max_deg) = truncate {
        res.retain(|((c, _), _)| degree(c.as_slice(), grading_vector) <= max_deg as i64);
    }
    let threefold = cy_kind.is_threefold();
    let label = |c: DVector<i32>, k: usize| {
        let degree = degree(c.as_slice(), grading_vector) as u32;
        (c.as_slice().to_vec(), (!threefold).then_some(k), degree)
    };
    Ok(match kind {
        InvariantKind::GV => Entries::Gv(
            res.into_iter()
                .map(|((c, k), v)| {
                    let (curve, surface, degree) = label(c, k);
                    Invariant {
                        curve,
                        surface,
                        degree,
                        value: v.numerator(),
                    }
                })
                .collect(),
        ),
        InvariantKind::GW => Entries::Gw(
            res.into_iter()
                .map(|((c, k), v)| {
                    let (curve, surface, degree) = label(c, k);
                    Invariant {
                        curve,
                        surface,
                        degree,
                        value: v.into_inner(),
                    }
                })
                .collect(),
        ),
    })
}

fn dot(a: &[i32], b: &[i32]) -> i64 {
    a.iter().zip(b).map(|(&x, &y)| x as i64 * y as i64).sum()
}

/// The curve classes of [`Selection::CurveClasses`]: the nonzero ones, once each, which must
/// have positive degree.
fn curve_classes(
    classes: &DMatrix<i32>,
    grading_vector: &RowDVector<i32>,
) -> Result<DMatrix<i32>, Error> {
    let h11 = classes.nrows();
    let mut seen = std::collections::HashSet::new();
    let mut res = Vec::new();
    for c in columns(classes) {
        if c.iter().all(|&x| x == 0) || !seen.insert(c) {
            continue;
        }
        if degree(c, grading_vector) <= 0 {
            return Err(input_err(
                "curve_classes",
                "the nonzero curve classes must have positive degree under the grading vector",
            ));
        }
        res.extend_from_slice(c);
    }
    if res.is_empty() {
        return Err(input_err(
            "curve_classes",
            "there must be at least one nonzero curve class",
        ));
    }
    Ok(DMatrix::from_vec(h11, res.len() / h11, res))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::dmatrix;

    #[test]
    fn test_default_grading() {
        let g = default_grading(&dmatrix![2, 0; 0, 2; 4, 2]).unwrap();
        assert_eq!(g.as_slice(), &[1, 1, 3]);
        assert_eq!(default_grading(&dmatrix![1, -1; 0, 0]), None);
    }
}
