//! The modular pipeline, a faster reorganization of the HKTY procedure (work in progress).
//!
//! It computes the invariants modulo a few word-sized primes and recovers them by the Chinese
//! remainder theorem, stores curve classes as packed integer keys, and never enumerates the
//! whole semigroup. The design follows the C program in PR #94 and is described in `PLAN.md`.
//!
//! The primitives the stages are built from are the arithmetic ([`field`]), the curve-class
//! keys ([`key`]) and the tables keyed by them ([`table`]). A pass, i.e. the computation
//! modulo a few primes, runs the stages in turn: the curve classes that carry period data
//! ([`candidates`]), the fundamental period at each of them ([`period`]), the series that the
//! invariants are extracted from ([`instanton`]), its coefficients in the basis of curves
//! ([`extraction`]) and the invariants from those ([`multicover`]). [`compute_gv`] repeats
//! passes with new primes until the lift of the residues ([`crt`]) is stable, and
//! [`compute_gw`] derives the GW invariants from the result.
//!
//! The entry points of [`dispatch`](crate::dispatch) reach it through
//! [`Method`](crate::Method). The candidates come either from the Hilbert bases of the cones
//! of [`crate::cone`] or from the enumeration of the semigroup of the generators. The stages
//! run on the threads of a rayon pool that each run builds for itself (see
//! [`Options::threads`]); the extraction, which dominates, applies the curves of each degree
//! in parallel (see [`extraction::extract`]).

pub mod candidates;
pub mod crt;
pub mod extraction;
pub mod field;
pub mod instanton;
pub mod key;
pub mod lightcone;
pub mod multicover;
pub mod period;
pub mod table;

#[cfg(test)]
mod reference_tests;

use crate::{check_strata, CYKind, Stratum, NON_INTEGER_CAUSES};
pub use candidates::CandidateSource;
use core::fmt;
use crt::Residues;
use field::{LaneField, Lanes};
use key::{KeyError, KeyLayout, PackedLayout};
use lightcone::Lightcone;
use malachite::base::num::logic::traits::SignificantBits;
use malachite::{Integer, Rational};
use nalgebra::{DMatrix, DVector};
use std::collections::{HashMap, HashSet};
use table::SparsePoly;

/// Errors of the modular pipeline.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    /// The curve classes do not fit the keys, or the grading vector is not positive.
    Key(KeyError),
    /// The input matrices and vectors do not have consistent dimensions.
    DimensionMismatch,
    /// A part of the nef partition names a divisor that does not exist.
    InvalidNefPartition,
    /// An intersection number has an index out of range, or is given more than once.
    InvalidIntersectionNumbers,
    /// A curve has negative intersection with an anticanonical class but fewer than two
    /// negative intersections with the toric divisors.
    NegativeAnticanonicalDegree,
    /// A value that has to be inverted is not invertible modulo one of the primes.
    NotInvertible,
    /// The options are out of range.
    InvalidOptions,
    /// The strata do not fit the geometry: it is not a threefold hypersurface, or their
    /// divisors are not three distinct ones, or their classes are out of range.
    InvalidStrata,
    /// The lift of the residues did not stabilize within the allowed number of passes.
    NotStable,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::Key(e) => write!(f, "{e}"),
            Self::DimensionMismatch => write!(f, "the dimensions of the input data do not match"),
            Self::InvalidNefPartition => write!(f, "the nef partition is inconsistent"),
            Self::InvalidIntersectionNumbers => write!(
                f,
                "the intersection numbers must have indices in range and be given only once"
            ),
            Self::NegativeAnticanonicalDegree => write!(
                f,
                "a curve has negative anticanonical degree but fewer than two negative \
                 intersections with the toric divisors"
            ),
            Self::NotInvertible => write!(f, "a value is not invertible modulo one of the primes"),
            Self::InvalidOptions => write!(f, "the options are out of range"),
            Self::InvalidStrata => write!(
                f,
                "the strata must be on three distinct divisors of a threefold hypersurface"
            ),
            Self::NotStable => write!(
                f,
                "the GV invariants did not stabilize as more primes were used, so they are \
                 probably not integers; {NON_INTEGER_CAUSES}"
            ),
        }
    }
}

impl std::error::Error for Error {}

impl From<KeyError> for Error {
    fn from(e: KeyError) -> Self {
        Self::Key(e)
    }
}

/// The integer data of the Calabi-Yau that the stages need, in the basis of curve classes.
///
/// - The GLSM charges $Q_r$, one row of $h^{1,1}$ integers per toric divisor.
/// - The anticanonical charges $Q_{0,j}$, one row per part of the nef partition, which is the
///   sum of the rows of its divisors; a hypersurface has a single row, the sum of all of them.
/// - $M$ symmetric bilinear forms $K^{(k)}_{ab}$ on the curve classes, which the second
///   derivatives of the period are contracted with, each giving one component of the
///   instanton series. A threefold has $M = 1$ with $K_{ab} = \sum_t w_t \kappa_{tab}$, $w$
///   the grading vector; an n-fold has one form $K^{(k)}_{ab} = \kappa_{kab}$ per reference
///   surface $k$.
///
/// - In a vex phase of a threefold hypersurface, the strata (see [`Stratum`]): each one's three
///   divisors and its class contracted with the grading vector, $w \cdot d_S$, which is what
///   the curves of negative anticanonical degree that meet exactly those divisors negatively
///   add to $S_2$ (see [`period`]).
///
/// The kind of the Calabi-Yau also decides how the invariants are read off the instanton
/// series (see [`multicover`]).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Geometry {
    kind: CYKind,
    h11: usize,
    q: Vec<i32>,
    q0: Vec<i32>,
    forms: Vec<i64>,
    strata: Vec<VexStratum>,
}

/// A stratum of a vex phase, as the period step uses it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VexStratum {
    /// The three divisors, in increasing order.
    pub divisors: [usize; 3],
    /// $w \cdot d_S = \sum_t w_t\, n_t / d_t$, as the terms $(w_t n_t, d_t)$ with $w_t n_t \ne 0$.
    pub weighted_class: Vec<(i64, i64)>,
}

impl Geometry {
    /// Assembles the data from flat row-major arrays: `q` and `q0` with rows of `h11` charges,
    /// and `forms` with `h11 * h11` entries per form, of which a threefold has exactly one.
    pub fn new(
        kind: CYKind,
        h11: usize,
        q: Vec<i32>,
        q0: Vec<i32>,
        forms: Vec<i64>,
    ) -> Result<Self, Error> {
        let n2 = h11 * h11;
        if h11 == 0
            || !q.len().is_multiple_of(h11)
            || q0.is_empty()
            || !q0.len().is_multiple_of(h11)
            || forms.is_empty()
            || !forms.len().is_multiple_of(n2)
            || (kind.is_threefold() && forms.len() != n2)
        {
            return Err(Error::DimensionMismatch);
        }
        let symmetric = forms
            .chunks_exact(n2)
            .all(|k| (0..h11).all(|a| (0..a).all(|b| k[a * h11 + b] == k[b * h11 + a])));
        if !symmetric {
            return Err(Error::InvalidIntersectionNumbers);
        }
        Ok(Self {
            kind,
            h11,
            q,
            q0,
            forms,
            strata: Vec::new(),
        })
    }

    /// Adds the strata of a vex phase, contracted with the grading vector. Only threefold
    /// hypersurfaces can have them, with distinct divisors that exist and positive
    /// denominators, and no two on the same divisors.
    pub fn with_strata(
        mut self,
        strata: &[Stratum],
        grading_vector: &[i32],
    ) -> Result<Self, Error> {
        if strata.is_empty() {
            return Ok(self);
        }
        if grading_vector.len() != self.h11 {
            return Err(Error::DimensionMismatch);
        }
        check_strata(
            strata,
            self.n_divisors(),
            self.h11,
            self.kind,
            self.n_anticanonical() == 1,
        )
        .map_err(|_| Error::InvalidStrata)?;
        for stratum in strata {
            let mut divisors = stratum.divisors;
            divisors.sort_unstable();
            let weighted_class = stratum
                .class
                .iter()
                .zip(grading_vector)
                .filter(|&(&(num, _), &w)| num != 0 && w != 0)
                .map(|(&(num, den), &w)| Some((num.checked_mul(w as i64)?, den)))
                .collect::<Option<_>>()
                .ok_or(Error::InvalidStrata)?;
            self.strata.push(VexStratum {
                divisors,
                weighted_class,
            });
        }
        Ok(self)
    }

    /// Derives the data from the inputs of the reference pipeline: the charge matrix `q` with
    /// one row per toric divisor, the nef partition (empty for a hypersurface), and the
    /// intersection numbers, with the grading vector for the threefold contraction.
    ///
    /// Each set of indices may appear only once, in any order for a threefold and in any
    /// order of the last two for an n-fold, whose first index is the reference surface.
    pub fn from_toric(
        q: &DMatrix<i32>,
        nefpart: &[DVector<i32>],
        intnums: &HashMap<(usize, usize, usize), i32>,
        cy_kind: CYKind,
        grading_vector: &[i32],
    ) -> Result<Self, Error> {
        let (ndiv, h11) = q.shape();
        if grading_vector.len() != h11 {
            return Err(Error::DimensionMismatch);
        }
        let rows: Vec<i32> = (0..ndiv)
            .flat_map(|r| (0..h11).map(move |a| q[(r, a)]))
            .collect();

        let mut q0 = Vec::new();
        if nefpart.is_empty() {
            q0.extend((0..h11).map(|a| q.column(a).sum()));
        } else {
            for part in nefpart {
                let mut row = vec![0; h11];
                for &r in part.iter() {
                    let r = usize::try_from(r)
                        .ok()
                        .filter(|&r| r < ndiv)
                        .ok_or(Error::InvalidNefPartition)?;
                    for (x, &c) in row.iter_mut().zip(&rows[r * h11..(r + 1) * h11]) {
                        *x += c;
                    }
                }
                q0.extend(row);
            }
        }

        let n2 = h11 * h11;
        let mut seen = HashSet::new();
        let forms = match cy_kind {
            CYKind::Threefold => {
                let mut form = vec![0i64; n2];
                for (&(i, j, k), &v) in intnums {
                    if i.max(j).max(k) >= h11 {
                        return Err(Error::InvalidIntersectionNumbers);
                    }
                    let mut idx = [i, j, k];
                    idx.sort_unstable();
                    if v == 0 {
                        continue;
                    }
                    if !seen.insert(idx) {
                        return Err(Error::InvalidIntersectionNumbers);
                    }
                    // Every distinct ordering (t, a, b) of the indices adds w_t kappa to K_ab.
                    let mut orderings = HashSet::new();
                    for [t, a, b] in [
                        [0, 1, 2],
                        [0, 2, 1],
                        [1, 0, 2],
                        [1, 2, 0],
                        [2, 0, 1],
                        [2, 1, 0],
                    ] {
                        if orderings.insert([idx[t], idx[a], idx[b]]) {
                            form[idx[a] * h11 + idx[b]] += grading_vector[idx[t]] as i64 * v as i64;
                        }
                    }
                }
                form
            }
            CYKind::Nfold => {
                let m = intnums.keys().map(|&(k, _, _)| k + 1).max().unwrap_or(1);
                let mut forms = vec![0i64; m * n2];
                for (&(k, a, b), &v) in intnums {
                    if a.max(b) >= h11 {
                        return Err(Error::InvalidIntersectionNumbers);
                    }
                    if v == 0 {
                        continue;
                    }
                    if !seen.insert([k, a.min(b), a.max(b)]) {
                        return Err(Error::InvalidIntersectionNumbers);
                    }
                    forms[k * n2 + a * h11 + b] = v as i64;
                    forms[k * n2 + b * h11 + a] = v as i64;
                }
                forms
            }
        };
        Self::new(cy_kind, h11, rows, q0, forms)
    }

    /// Whether the Calabi-Yau is a threefold or of higher dimension.
    #[inline]
    pub fn kind(&self) -> CYKind {
        self.kind
    }

    /// The dimension of the space of curve classes.
    #[inline]
    pub fn h11(&self) -> usize {
        self.h11
    }

    /// The number of toric divisors.
    #[inline]
    pub fn n_divisors(&self) -> usize {
        self.q.len() / self.h11
    }

    /// The charges of divisor `r`.
    #[inline]
    pub fn q(&self, r: usize) -> &[i32] {
        &self.q[r * self.h11..(r + 1) * self.h11]
    }

    /// The number of anticanonical rows, i.e. of parts of the nef partition.
    #[inline]
    pub fn n_anticanonical(&self) -> usize {
        self.q0.len() / self.h11
    }

    /// The anticanonical charges of part `j`.
    #[inline]
    pub fn q0(&self, j: usize) -> &[i32] {
        &self.q0[j * self.h11..(j + 1) * self.h11]
    }

    /// The number $M$ of bilinear forms, i.e. of components of the instanton series.
    #[inline]
    pub fn n_forms(&self) -> usize {
        self.forms.len() / (self.h11 * self.h11)
    }

    /// Form `k`, as a row-major $h^{1,1} \times h^{1,1}$ matrix.
    #[inline]
    pub fn form(&self, k: usize) -> &[i64] {
        let n2 = self.h11 * self.h11;
        &self.forms[k * n2..(k + 1) * n2]
    }

    /// The number of negative intersections of the curve class `c` with the toric divisors.
    pub fn n_negative(&self, c: &[i32]) -> usize {
        (0..self.n_divisors())
            .filter(|&r| dot(c, self.q(r)) < 0)
            .count()
    }

    /// The strata of a vex phase, if any.
    #[inline]
    pub fn strata(&self) -> &[VexStratum] {
        &self.strata
    }

    /// The index of the stratum on the divisors `neg`, given in increasing order, if any.
    pub fn stratum(&self, neg: &[usize]) -> Option<usize> {
        self.strata.iter().position(|s| s.divisors[..] == *neg)
    }

    /// Whether the fundamental period can be nonzero at the curve class `c`: when it has at most
    /// two negative intersections with the toric divisors, or when it has negative
    /// anticanonical degree and its negative intersections are with the divisors of a stratum.
    pub fn carries_period(&self, c: &[i32]) -> bool {
        let mut neg = [0; 3];
        let mut n = 0;
        for r in 0..self.n_divisors() {
            if dot(c, self.q(r)) < 0 {
                if n == 3 {
                    return false;
                }
                neg[n] = r;
                n += 1;
            }
        }
        n <= 2
            || (self.stratum(&neg).is_some()
                && (0..self.n_anticanonical()).any(|j| dot(c, self.q0(j)) < 0))
    }
}

/// The options of [`compute_gv`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Options {
    /// The number of primes that only check the lift, at least one: the lift from the other
    /// primes is accepted when it agrees with these at every curve.
    pub check_primes: usize,
    /// The number of primes per pass, from 2 to 4, or `None` to choose it from the size of the
    /// invariants: up to degree 12 two, and above that as many as the largest invariant needs,
    /// predicted from a run at 3/4 of the degree, plus the check primes, as the C program's
    /// driver does.
    pub lanes: Option<usize>,
    /// The minimum number of primes, even if fewer give a stable lift.
    pub min_primes: usize,
    /// The number of passes after which an unstable lift is an error.
    pub max_passes: usize,
    /// The number of worker threads, at least one, or `None` for one per core. The run
    /// builds a thread pool of its own rather than use rayon's global one.
    pub threads: Option<usize>,
    /// How the extraction parallelizes.
    pub extraction: extraction::Tuning,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            check_primes: 1,
            lanes: None,
            min_primes: 0,
            max_passes: 8,
            threads: None,
            extraction: extraction::Tuning::default(),
        }
    }
}

/// The invariants of one curve class, integers for GV and rationals for GW.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CurveInvariants<T = Integer> {
    pub curve: Vec<i32>,
    pub degree: u32,
    /// One value per component of the instanton series: a single one for a threefold, and
    /// one per reference surface for an n-fold.
    pub values: Vec<T>,
}

/// Computes the GV invariants of the curve classes of degree at most `max_deg` that the
/// candidates of `source` give, in increasing order of degree and then key, leaving out those
/// where all of them vanish.
///
/// With `targets`, only the curve classes in their backward lightcones are computed, i.e.
/// those that some target minus the class leaves in the Mori cone (see [`lightcone`]), and
/// `max_deg` is lowered to the largest degree of the targets.
pub fn compute_gv(
    source: CandidateSource<'_>,
    grading_vector: &[i32],
    max_deg: u32,
    targets: &[&[i32]],
    geometry: &Geometry,
    options: &Options,
) -> Result<Vec<CurveInvariants>, Error> {
    run(
        source,
        grading_vector,
        max_deg,
        targets,
        geometry,
        options,
        |gv, _, _| gv,
    )
}

/// Computes the GW invariants of the curve classes of degree at most `max_deg` that the
/// candidates of `source` give, in increasing order of degree and then curve class, leaving
/// out those where all of them vanish, and restricted to the backward lightcones of `targets`
/// as in [`compute_gv`]. They come from the GV invariants of [`compute_gv`], by
/// [`multicover::gw_from_gv`].
pub fn compute_gw(
    source: CandidateSource<'_>,
    grading_vector: &[i32],
    max_deg: u32,
    targets: &[&[i32]],
    geometry: &Geometry,
    options: &Options,
) -> Result<Vec<CurveInvariants<Rational>>, Error> {
    run(
        source,
        grading_vector,
        max_deg,
        targets,
        geometry,
        options,
        |gv, max_deg, in_lightcone| {
            let mut gw = multicover::gw_from_gv(geometry.kind(), &gv, max_deg);
            // Every divisor of a curve in the lightcones is in them, so the GW invariants are
            // right there, but not at the multiples outside them.
            gw.retain(|inv| in_lightcone(&inv.curve));
            gw
        },
    )
}

/// Computes the GV invariants as [`compute_gv`] describes, and hands them to `finish` with the
/// maximum degree, lowered to the largest degree of the targets if there are any, and a test
/// of whether a curve class of at most that degree is in their backward lightcones.
fn run<T>(
    source: CandidateSource<'_>,
    grading_vector: &[i32],
    max_deg: u32,
    targets: &[&[i32]],
    geometry: &Geometry,
    options: &Options,
    finish: impl FnOnce(Vec<CurveInvariants>, u32, &dyn Fn(&[i32]) -> bool) -> T,
) -> Result<T, Error> {
    if options.check_primes == 0
        || options.lanes.is_some_and(|l| !(2..=4).contains(&l))
        || options.threads == Some(0)
    {
        return Err(Error::InvalidOptions);
    }
    let threads = options.threads.unwrap_or_else(|| {
        std::thread::available_parallelism().map_or(1, std::num::NonZeroUsize::get)
    });
    // The stages parallelize with rayon and run in whichever pool they are installed in, so
    // a private one honours `threads` and keeps out of a library consumer's own use of rayon.
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("failed to build the worker thread pool");
    let max_deg = if targets.is_empty() {
        max_deg
    } else {
        max_deg.min(lightcone::max_degree(targets, grading_vector)?)
    };
    let bounds = source.coordinate_bounds(grading_vector, max_deg, targets)?;
    let layout = PackedLayout::new(&bounds)?;
    let candidates::Candidates {
        points: cands,
        lightcone,
    } = pool.install(|| source.candidates(&layout, grading_vector, max_deg, targets, geometry))?;
    let in_lightcone = |c: &[i32]| layout.fits(c) && lightcone.contains(layout.pack(c), c);
    // Without zero, which is the only class of degree zero, nothing is below the targets.
    if cands.first().is_none_or(|&(_, deg)| deg != 0) {
        return Ok(finish(Vec::new(), max_deg, &in_lightcone));
    }
    let run = Run {
        layout: &layout,
        geometry,
        cands: &cands,
        max_deg,
        lightcone: &lightcone,
        tuning: &options.extraction,
    };
    let gv = pool.install(|| {
        let lanes = match options.lanes {
            Some(lanes) => lanes,
            None => choose_lanes(&run, options)?,
        };
        let residues = match lanes {
            2 => passes::<2, _>(&run, options),
            3 => passes::<3, _>(&run, options),
            _ => passes::<4, _>(&run, options),
        }?;
        let mut curve = vec![0; layout.dim()];
        Ok::<_, Error>(
            residues
                .lift()
                .into_iter()
                .map(|(key, degree, values)| {
                    layout.unpack(key, &mut curve);
                    CurveInvariants {
                        curve: curve.clone(),
                        degree,
                        values,
                    }
                })
                .collect(),
        )
    })?;
    Ok(finish(gv, max_deg, &in_lightcone))
}

/// The number of primes per pass for a run whose options leave it open, chosen the way the C
/// program's driver does: enough for the largest invariant, plus the check primes.
///
/// Up to degree 12 that is taken to be two. Above, a probe computes the invariants up to 3/4 of
/// the maximum degree, with two primes per pass, and the largest number of bits per degree
/// over its upper half is extrapolated linearly to the maximum degree, with a margin that grows
/// with the scatter of the fit; $k$ primes of 62 bits lift values below about $2^{62k - 2}$.
fn choose_lanes<L: KeyLayout>(run: &Run<'_, L>, options: &Options) -> Result<usize, Error> {
    let max_deg = run.max_deg;
    if max_deg <= 12 {
        return Ok(2);
    }
    let probe_deg = 3 * max_deg / 4;
    if probe_deg < 8 {
        return Ok(3);
    }
    let n = run.cands.partition_point(|&(_, deg)| deg <= probe_deg);
    let probe = Run {
        cands: &run.cands[..n],
        max_deg: probe_deg,
        ..*run
    };
    let probe_options = Options {
        min_primes: 0,
        ..options.clone()
    };
    let residues = passes::<2, _>(&probe, &probe_options)?;
    // The largest number of bits of an invariant of each degree, i.e. floor(log2 |GV|).
    let mut bits: Vec<Option<f64>> = vec![None; probe_deg as usize + 1];
    for (_, deg, values) in residues.lift() {
        for v in values {
            let n = v.significant_bits();
            if n > 0 {
                let b = (n - 1) as f64;
                let max = &mut bits[deg as usize];
                *max = Some(max.map_or(b, |m: f64| m.max(b)));
            }
        }
    }
    let points: Vec<(f64, f64)> = (probe_deg / 2..=probe_deg)
        .filter_map(|d| bits[d as usize].map(|b| (d as f64, b)))
        .collect();
    let Some(&(last_deg, last_bits)) = points.last() else {
        return Ok(2);
    };
    // The least-squares line through the points, and the scatter around it.
    let n = points.len() as f64;
    let (sx, sy) = points
        .iter()
        .fold((0.0, 0.0), |(sx, sy), &(x, y)| (sx + x, sy + y));
    let (sxx, sxy) = points
        .iter()
        .fold((0.0, 0.0), |(sxx, sxy), &(x, y)| (sxx + x * x, sxy + x * y));
    let slope = if points.len() >= 2 {
        ((n * sxy - sx * sy) / (n * sxx - sx * sx)).max(0.0)
    } else {
        1.0
    };
    let intercept = (sy - slope * sx) / n;
    let rmse = (points
        .iter()
        .map(|&(x, y)| (y - intercept - slope * x).powi(2))
        .sum::<f64>()
        / n)
        .sqrt();
    let predicted = last_bits + slope * (max_deg as f64 - last_deg) + 4.0 + 2.0 * rmse;
    let mut k = 1;
    while 62.0 * k as f64 - 2.0 < predicted {
        k += 1;
    }
    Ok((k + options.check_primes).clamp(2, 4))
}

/// What the passes of one run share.
struct Run<'a, L: KeyLayout> {
    layout: &'a L,
    geometry: &'a Geometry,
    cands: &'a [(L::Key, u32)],
    max_deg: u32,
    lightcone: &'a Lightcone<L::Key>,
    tuning: &'a extraction::Tuning,
}

/// Runs passes over `NL` primes each, taken in order from [`field::primes`], until the lift
/// is stable.
fn passes<const NL: usize, L: KeyLayout>(
    run: &Run<'_, L>,
    options: &Options,
) -> Result<Residues<L::Key>, Error> {
    let mut residues = Residues::new(run.geometry.n_forms());
    let mut primes = field::primes();
    for _ in 0..options.max_passes {
        let ring = LaneField::<NL>::new(core::array::from_fn(|_| primes.next().unwrap()));
        residues.add_pass(&ring, &pass(&ring, run)?);
        if residues.primes().len() >= options.min_primes && residues.is_stable(options.check_primes)
        {
            return Ok(residues);
        }
    }
    Err(Error::NotStable)
}

/// The GV invariants modulo the primes of `ring`.
fn pass<const NL: usize, L: KeyLayout>(
    ring: &LaneField<NL>,
    run: &Run<'_, L>,
) -> Result<SparsePoly<L::Key, Lanes<NL>>, Error> {
    let Run {
        layout,
        geometry,
        cands,
        max_deg,
        lightcone,
        tuning,
    } = *run;
    let pd = period::period(ring, layout, geometry, cands, max_deg)?;
    let inst = instanton::instanton(ring, geometry, &pd, max_deg)?;
    let inv = &pd.tables.inv;
    let coeffs = extraction::extract(
        ring,
        layout,
        &inst.alpha,
        &inst.inst,
        inv,
        max_deg,
        lightcone,
        tuning,
    )?;
    let mut gv = multicover::multicover(ring, layout, geometry.kind(), &coeffs, inv, max_deg)?;
    // The products of curves in the lightcones reach curves outside them, where the values
    // miss the curves that were dropped.
    if !lightcone.is_all() {
        let mut curve = vec![0; layout.dim()];
        gv.retain(|key| {
            layout.unpack(key, &mut curve);
            lightcone.contains(key, &curve)
        });
    }
    Ok(gv)
}

/// The intersection of a curve class with a divisor, computed in `i64`.
#[inline]
pub(crate) fn dot(c: &[i32], d: &[i32]) -> i64 {
    c.iter().zip(d).map(|(&x, &y)| x as i64 * y as i64).sum()
}

/// A small deterministic pseudorandom generator for the unit tests (SplitMix64).
#[cfg(test)]
pub(crate) struct TestRng(u64);

#[cfg(test)]
impl TestRng {
    pub(crate) fn new(seed: u64) -> Self {
        Self(seed)
    }

    pub(crate) fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// A uniformly distributed integer in `lo..=hi`.
    pub(crate) fn range(&mut self, lo: i64, hi: i64) -> i64 {
        let span = (hi - lo) as u64 + 1;
        lo + (self.next_u64() % span) as i64
    }
}

/// The residue modulo a prime `p` of a rational whose denominator `p` does not divide.
#[cfg(test)]
pub(crate) fn residue(x: &Rational, p: u64) -> u64 {
    use malachite::base::num::arithmetic::traits::{ModInverse, ModMul};
    use malachite::Natural;
    let reduce = |n: &Natural| u64::try_from(&(n % Natural::from(p))).unwrap();
    let num = reduce(x.numerator_ref());
    let num = if *x < 0u32 { (p - num) % p } else { num };
    let den = reduce(x.denominator_ref());
    num.mod_mul(den.mod_inverse(p).unwrap(), p)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn toy_q() -> DMatrix<i32> {
        DMatrix::from_column_slice(6, 2, &[1, 1, 1, 0, 1, 2, 0, 0, -1, 1, 1, -1])
    }

    #[test]
    fn test_geometry_threefold() {
        let intnums = HashMap::from([
            ((0, 0, 0), 2),
            ((1, 0, 0), 1),
            ((0, 1, 1), -1),
            ((1, 1, 1), 5),
        ]);
        let g = Geometry::from_toric(&toy_q(), &[], &intnums, CYKind::Threefold, &[3, -1]).unwrap();
        assert_eq!(
            (g.h11(), g.n_divisors(), g.n_anticanonical(), g.n_forms()),
            (2, 6, 1, 1)
        );
        assert_eq!(g.q(2), &[1, -1]);
        assert_eq!(g.q0(0), &[6, 0]);
        // K_ab = sum_t w_t kappa_tab, e.g. K_01 = 3 kappa_001 - kappa_101 = 3 + 1.
        assert_eq!(g.form(0), &[5, 4, 4, -8]);
        assert_eq!(g.n_negative(&[1, 2]), 1);
        assert_eq!(g.n_negative(&[0, -1]), 2);

        let nefpart = [
            DVector::from_column_slice(&[0, 2]),
            DVector::from_column_slice(&[5]),
        ];
        let g = Geometry::from_toric(&toy_q(), &nefpart, &intnums, CYKind::Threefold, &[3, -1])
            .unwrap();
        assert_eq!((g.q0(0), g.q0(1)), (&[2, -1][..], &[2, -1][..]));
    }

    #[test]
    fn test_geometry_nfold() {
        let intnums = HashMap::from([((0, 0, 1), 2), ((2, 1, 1), 3), ((2, 0, 0), 0)]);
        let g = Geometry::from_toric(&toy_q(), &[], &intnums, CYKind::Nfold, &[3, -1]).unwrap();
        assert_eq!(g.n_forms(), 3);
        assert_eq!(g.form(0), &[0, 2, 2, 0]);
        assert_eq!(g.form(1), &[0; 4]);
        assert_eq!(g.form(2), &[0, 0, 0, 3]);
    }

    #[test]
    fn test_geometry_errors() {
        let q = toy_q();
        let ok = HashMap::from([((0, 0, 0), 1)]);
        let three = CYKind::Threefold;
        let err = |nefpart: &[DVector<i32>], intnums: &HashMap<_, _>, kind, w: &[i32]| {
            Geometry::from_toric(&q, nefpart, intnums, kind, w).unwrap_err()
        };
        assert_eq!(err(&[], &ok, three, &[1]), Error::DimensionMismatch);
        assert_eq!(
            err(&[DVector::from_column_slice(&[6])], &ok, three, &[1, 1]),
            Error::InvalidNefPartition
        );
        let repeated = HashMap::from([((0, 0, 1), 1), ((0, 1, 0), 1)]);
        assert_eq!(
            err(&[], &repeated, three, &[1, 1]),
            Error::InvalidIntersectionNumbers
        );
        let repeated = HashMap::from([((0, 0, 1), 1), ((0, 1, 0), 1)]);
        assert_eq!(
            err(&[], &repeated, CYKind::Nfold, &[1, 1]),
            Error::InvalidIntersectionNumbers
        );
        let out_of_range = HashMap::from([((0, 0, 2), 1)]);
        assert_eq!(
            err(&[], &out_of_range, three, &[1, 1]),
            Error::InvalidIntersectionNumbers
        );
        // The same order with a zero value is ignored.
        assert!(Geometry::from_toric(
            &q,
            &[],
            &HashMap::from([((0, 0, 1), 1), ((0, 1, 0), 0)]),
            three,
            &[1, 1]
        )
        .is_ok());

        assert_eq!(
            Geometry::new(CYKind::Nfold, 2, vec![1, 0], vec![1, 0], vec![0, 1, 2, 0]),
            Err(Error::InvalidIntersectionNumbers)
        );
        assert_eq!(
            Geometry::new(CYKind::Nfold, 2, vec![1, 0, 1], vec![1, 0], vec![0; 4]),
            Err(Error::DimensionMismatch)
        );
        assert_eq!(
            Geometry::new(CYKind::Nfold, 2, vec![1, 0], vec![], vec![0; 4]),
            Err(Error::DimensionMismatch)
        );
        assert_eq!(
            Geometry::new(CYKind::Nfold, 2, vec![1, 0], vec![1, 0], vec![0; 3]),
            Err(Error::DimensionMismatch)
        );
        assert_eq!(
            Geometry::new(CYKind::Nfold, 0, vec![], vec![], vec![]),
            Err(Error::DimensionMismatch)
        );
    }
}
