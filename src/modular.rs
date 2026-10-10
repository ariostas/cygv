//! The modular pipeline, a faster reorganization of the HKTY procedure (work in progress).
//!
//! It computes the invariants modulo a few word-sized primes and recovers them by the Chinese
//! remainder theorem, stores curve classes as packed integer keys, and never enumerates the
//! whole semigroup. The design follows the C program in PR #94 and is described in `PLAN.md`.
//!
//! So far this module holds the primitives the stages are built from, the arithmetic
//! ([`field`]), the curve-class keys ([`key`]) and the tables keyed by them ([`table`]), and
//! the first three stages: the curve classes that carry period data ([`candidates`]), the
//! fundamental period at each of them ([`period`]), and the series that the invariants are
//! extracted from ([`instanton`]).

pub mod candidates;
pub mod field;
pub mod instanton;
pub mod key;
pub mod period;
pub mod table;

#[cfg(test)]
mod reference_tests;

use crate::CYKind;
use core::fmt;
use key::KeyError;
use nalgebra::{DMatrix, DVector};
use std::collections::{HashMap, HashSet};

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
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Geometry {
    h11: usize,
    q: Vec<i32>,
    q0: Vec<i32>,
    forms: Vec<i64>,
}

impl Geometry {
    /// Assembles the data from flat row-major arrays: `q` and `q0` with rows of `h11` charges,
    /// and `forms` with `h11 * h11` entries per form.
    pub fn new(h11: usize, q: Vec<i32>, q0: Vec<i32>, forms: Vec<i64>) -> Result<Self, Error> {
        let n2 = h11 * h11;
        if h11 == 0
            || !q.len().is_multiple_of(h11)
            || q0.is_empty()
            || !q0.len().is_multiple_of(h11)
            || forms.is_empty()
            || !forms.len().is_multiple_of(n2)
        {
            return Err(Error::DimensionMismatch);
        }
        let symmetric = forms
            .chunks_exact(n2)
            .all(|k| (0..h11).all(|a| (0..a).all(|b| k[a * h11 + b] == k[b * h11 + a])));
        if !symmetric {
            return Err(Error::InvalidIntersectionNumbers);
        }
        Ok(Self { h11, q, q0, forms })
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
        Self::new(h11, rows, q0, forms)
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
            Geometry::new(2, vec![1, 0], vec![1, 0], vec![0, 1, 2, 0]),
            Err(Error::InvalidIntersectionNumbers)
        );
        assert_eq!(
            Geometry::new(2, vec![1, 0, 1], vec![1, 0], vec![0; 4]),
            Err(Error::DimensionMismatch)
        );
        assert_eq!(
            Geometry::new(2, vec![1, 0], vec![], vec![0; 4]),
            Err(Error::DimensionMismatch)
        );
        assert_eq!(
            Geometry::new(2, vec![1, 0], vec![1, 0], vec![0; 3]),
            Err(Error::DimensionMismatch)
        );
        assert_eq!(
            Geometry::new(0, vec![], vec![], vec![]),
            Err(Error::DimensionMismatch)
        );
    }
}
