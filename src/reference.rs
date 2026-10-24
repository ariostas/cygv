//! The reference pipeline, a faithful implementation of the HKTY procedure.

pub mod factorial;
pub mod fundamental_period;
pub mod instanton;
pub mod misc;
pub mod polynomial;
pub mod rational;
pub mod semigroup;
pub mod series_inversion;

use crate::{CYKind, InvariantKind, Stratum};
use core::fmt;
use fundamental_period::error::FundamentalPeriodError;
use misc::error::MiscError;
use nalgebra::{DMatrix, DVector, RowDVector};
use polynomial::error::PolynomialError;
use polynomial::{coefficient::PolynomialCoeff, properties::PolynomialProperties};
use semigroup::error::SemigroupError;
use semigroup::Semigroup;
use series_inversion::error::SeriesInversionError;
use std::collections::HashMap;
use std::thread;

#[doc(inline)]
pub use rational::Rational;

/// The errors of [`run_hkty`], by the stage that failed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    Semigroup(SemigroupError),
    IntersectionNumbers(MiscError),
    FundamentalPeriod(FundamentalPeriodError),
    Instanton(PolynomialError),
    SeriesInversion(SeriesInversionError),
    /// The computation was interrupted ([`crate::Interrupt`]).
    Interrupted,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::Semigroup(e) => write!(f, "{e}"),
            Self::IntersectionNumbers(e) => write!(f, "{e}"),
            Self::FundamentalPeriod(e) => write!(f, "{e}"),
            Self::Instanton(e) => write!(f, "{e}"),
            Self::SeriesInversion(e) => write!(f, "{e}"),
            Self::Interrupted => write!(f, "the computation was interrupted"),
        }
    }
}

impl std::error::Error for Error {}

/// The invariants of each curve class and reference surface (0 for a threefold).
pub type Invariants<T> = Vec<((DVector<i32>, usize), T)>;

/// Runs the HKTY procedure on the semigroup of `generators`, truncated by the first of
/// `target_points`, `max_deg` and `min_points` that is given (or not at all, the generators
/// then being every curve class), and returns the invariants of every curve class and, for an
/// n-fold, reference surface, in order of the semigroup's elements.
///
/// The public entry points check the input before calling it; an inconsistent geometry is an
/// error, from the stage that finds it.
#[allow(clippy::too_many_arguments)]
pub fn run_hkty<T>(
    generators: DMatrix<i32>,
    grading_vector: RowDVector<i32>,
    invariant_kind: InvariantKind,
    cy_kind: CYKind,
    max_deg: Option<u32>,
    min_points: Option<u32>,
    target_points: Option<DMatrix<i32>>,
    q: DMatrix<i32>,
    nefpart: Vec<DVector<i32>>,
    intnums: HashMap<(usize, usize, usize), i32>,
    strata: Vec<Stratum>,
    n_threads: Option<u32>,
) -> Result<Invariants<T>, Error>
where
    T: PolynomialCoeff<T>,
{
    let sg = if let Some(target_points) = target_points {
        Semigroup::with_target_points(generators, grading_vector, target_points)
    } else if let Some(d) = max_deg {
        Semigroup::with_max_degree(generators, grading_vector, d)
    } else if let Some(n) = min_points {
        Semigroup::with_min_elements(generators, grading_vector, n as usize)
    } else {
        Semigroup::from_data(generators, grading_vector)
    }
    .map_err(Error::Semigroup)?;
    log::info!(
        "semigroup: {} curve classes up to degree {}",
        sg.elements.ncols(),
        sg.max_degree
    );

    let poly_props = PolynomialProperties::new(&sg);

    let (intnum_dict, intnum_idxpairs, n_indices) =
        misc::process_int_nums(intnums, cy_kind).map_err(Error::IntersectionNumbers)?;
    // A threefold's instanton corrections are indexed by every coordinate of the curves,
    // whichever the intersection numbers reach.
    let n_indices = match cy_kind {
        CYKind::Threefold => n_indices.max(sg.elements.nrows()),
        CYKind::Nfold => n_indices,
    };

    let n_threads = match n_threads {
        None => thread::available_parallelism()
            .unwrap_or(std::num::NonZeroUsize::new(1).unwrap())
            .get(),
        Some(0) => 1,
        Some(n) => n as usize,
    };

    // The stages parallelize with rayon and so run in whichever pool they are
    // installed in. A private pool is built here rather than reaching for the
    // global one, both to honour `n_threads` and to keep a library consumer's
    // own use of rayon out of the picture.
    let pool = crate::interrupt::pool(n_threads)
        .build()
        .expect("failed to build the worker thread pool");

    let gv = pool.install(|| {
        let fp = fundamental_period::compute_omega(
            &poly_props,
            &sg,
            &q,
            &nefpart,
            &intnum_idxpairs,
            cy_kind,
            &strata,
        )
        .map_err(Error::FundamentalPeriod)?;
        log::info!("fundamental period done");
        crate::interrupt::check(Error::Interrupted)?;

        let inst_data = instanton::compute_instanton_data(
            fp,
            &poly_props,
            &intnum_idxpairs,
            n_indices,
            &intnum_dict,
            cy_kind,
        )
        .map_err(Error::Instanton)?;
        log::info!("instanton corrections done");
        crate::interrupt::check(Error::Interrupted)?;

        series_inversion::invert_series(inst_data, &poly_props, invariant_kind, cy_kind).map_err(
            |e| match e {
                SeriesInversionError::Interrupted => Error::Interrupted,
                e => Error::SeriesInversion(e),
            },
        )
    })?;

    let mut gv_sorted: Vec<_> = gv.into_iter().collect();
    gv_sorted.sort_unstable_by_key(|c| c.0 .0);
    Ok(gv_sorted
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
        .collect())
}
