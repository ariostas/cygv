//! The reference pipeline, a faithful implementation of the HKTY procedure.

pub mod factorial;
pub mod fundamental_period;
pub mod instanton;
pub mod misc;
pub mod polynomial;
pub mod semigroup;
pub mod series_inversion;

use crate::{CYKind, InvariantKind, Stratum};
use nalgebra::{DMatrix, DVector, RowDVector};
use polynomial::{coefficient::PolynomialCoeff, properties::PolynomialProperties};
use semigroup::Semigroup;
use std::collections::HashMap;
use std::thread;

// TODO: Need to handle errors properly
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
) -> Vec<((DVector<i32>, usize), T)>
where
    T: PolynomialCoeff<T>,
{
    let sg = if let Some(target_points) = target_points {
        Semigroup::with_target_points(generators, grading_vector, target_points).unwrap()
    } else if let Some(d) = max_deg {
        Semigroup::with_max_degree(generators, grading_vector, d).unwrap()
    } else if let Some(n) = min_points {
        Semigroup::with_min_elements(generators, grading_vector, n as usize).unwrap()
    } else {
        Semigroup::from_data(generators, grading_vector).unwrap()
    };

    let poly_props = PolynomialProperties::new(&sg);

    let (intnum_dict, intnum_idxpairs, n_indices) =
        misc::process_int_nums(intnums, cy_kind).unwrap();

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
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
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
        .unwrap();

        let inst_data = instanton::compute_instanton_data(
            fp,
            &poly_props,
            &intnum_idxpairs,
            n_indices,
            &intnum_dict,
            cy_kind,
        )
        .unwrap();

        series_inversion::invert_series(inst_data, &poly_props, invariant_kind, cy_kind).unwrap()
    });

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
