//! Shared fixtures and scenario definitions for the `cygv` benchmarks.
//!
//! This module is compiled into every bench target, and not every target uses
//! every item, hence the blanket `dead_code` allowance.

#![allow(dead_code)]

use cygv::hkty::{
    compute_gv_float_nfold, compute_gv_float_threefold, compute_gv_rat_nfold,
    compute_gv_rat_threefold, compute_gw_float_nfold, compute_gw_float_threefold,
    compute_gw_rat_nfold, compute_gw_rat_threefold,
};
use nalgebra::{dmatrix, dvector, DMatrix, DVector, RowDVector};
use std::collections::HashMap;

/// Size of the [`cygv::NumberPool`]s used by the benchmarks.
pub const POOL_SIZE: usize = 1000;

/// Number of threads the benchmarks run with.
///
/// Defaults to `None`, i.e. one thread per available core. Set
/// `CYGV_BENCH_THREADS` to pin it, which makes results less noisy and easier to
/// compare across machines.
pub fn n_threads() -> Option<u32> {
    match std::env::var("CYGV_BENCH_THREADS") {
        Ok(s) => Some(
            s.parse()
                .expect("CYGV_BENCH_THREADS must be a non-negative integer"),
        ),
        Err(_) => None,
    }
}

/// The input data describing a Calabi-Yau geometry.
#[derive(Clone)]
pub struct Model {
    pub name: &'static str,
    pub generators: DMatrix<i32>,
    pub grading_vector: RowDVector<i32>,
    pub q: DMatrix<i32>,
    pub nefpart: Vec<DVector<i32>>,
    pub intnums: HashMap<(usize, usize, usize), i32>,
    pub is_threefold: bool,
}

/// Looks up one of the benchmark models by name.
pub fn model(name: &str) -> Model {
    match name {
        "threefold" => threefold(),
        "fourfold" => fourfold(),
        _ => panic!("unknown model {name:?}"),
    }
}

/// A hypersurface CY threefold with $h^{1,1} = 2$.
pub fn threefold() -> Model {
    let generators = dmatrix![
         0, 1;
        -1, 2;
    ];
    let grading_vector = RowDVector::from_row_slice(&[3, -1]);

    let q = dmatrix![
        1,  0;
        1,  0;
        1, -1;
        0,  1;
        1,  1;
        2, -1;
    ];
    let nefpart = Vec::new();

    let intnums = HashMap::from([
        ((0, 0, 0), 2),
        ((0, 0, 1), 1),
        ((0, 1, 1), -1),
        ((1, 1, 1), -5),
    ]);

    Model {
        name: "threefold",
        generators,
        grading_vector,
        q,
        nefpart,
        intnums,
        is_threefold: true,
    }
}

/// A complete intersection CY fourfold with $h^{1,1} = 6$.
pub fn fourfold() -> Model {
    let generators = dmatrix![
        0,1,0,0,0,0;
        0,0,0,1,-2,0;
        0,-4,1,0,0,0;
        -2,0,1,0,1,0;
        3,0,0,1,0,-2;
        0,0,0,-2,2,1;
    ];
    let grading_vector = RowDVector::from_row_slice(&[73, -30, 18, -17, -11, -21]);

    let q = dmatrix![
        1,0,0,0,0,0;
        1,0,0,0,0,0;
        1,0,0,0,0,0;
        1,0,0,0,0,0;
        0,1,0,0,0,0;
        0,0,1,0,0,0;
        8,-4,2,-2,-2,-3;
        12,-6,3,-3,-2,-4;
        4,-3,0,-2,-1,-2;
        0,0,0,1,0,0;
        0,0,0,0,1,0;
        0,0,0,0,0,1;
    ];
    let nefpart = vec![dvector![0, 1, 2, 3], dvector![4, 5, 6, 7, 8, 9, 10, 11]];

    let intnums = HashMap::from([
        ((0, 0, 2), 1),
        ((0, 1, 1), -8),
        ((0, 1, 3), 4),
        ((0, 1, 5), 8),
        ((0, 2, 2), -4),
        ((0, 3, 3), -8),
        ((0, 4, 4), -16),
        ((0, 4, 5), 8),
        ((0, 5, 5), -16),
        ((1, 0, 1), -8),
        ((1, 0, 3), 4),
        ((1, 0, 5), 8),
        ((1, 1, 3), 16),
        ((1, 1, 5), -32),
        ((1, 3, 3), -16),
        ((1, 5, 5), 64),
        ((2, 0, 0), 1),
        ((2, 0, 2), -4),
        ((2, 2, 2), 16),
        ((3, 0, 1), 4),
        ((3, 0, 3), -8),
        ((3, 1, 1), 16),
        ((3, 1, 3), -16),
        ((4, 0, 4), -16),
        ((4, 0, 5), 8),
        ((4, 4, 4), -128),
        ((4, 4, 5), 32),
        ((5, 0, 1), 8),
        ((5, 0, 4), 8),
        ((5, 0, 5), -16),
        ((5, 1, 1), -32),
        ((5, 1, 5), 64),
        ((5, 4, 4), 32),
        ((5, 5, 5), -128),
        ((6, 0, 0), -8),
        ((6, 0, 3), 16),
        ((6, 0, 5), -32),
        ((6, 1, 1), -128),
        ((6, 1, 3), 64),
        ((6, 1, 5), 128),
        ((6, 3, 3), -64),
        ((6, 5, 5), -256),
        ((7, 0, 0), 4),
        ((7, 0, 1), 16),
        ((7, 0, 3), -16),
        ((7, 1, 1), 64),
        ((7, 1, 3), -64),
        ((7, 3, 3), 64),
        ((8, 0, 0), 8),
        ((8, 0, 1), -32),
        ((8, 0, 5), 64),
        ((8, 1, 1), 128),
        ((8, 1, 5), -256),
        ((8, 5, 5), 512),
        ((9, 0, 0), -4),
        ((9, 0, 2), 16),
        ((9, 2, 2), -64),
        ((10, 0, 0), -8),
        ((10, 0, 1), -16),
        ((10, 1, 1), -64),
        ((10, 1, 3), 64),
        ((10, 3, 3), -128),
        ((11, 0, 0), -16),
        ((11, 0, 4), -128),
        ((11, 0, 5), 32),
        ((11, 4, 4), -768),
        ((11, 4, 5), 128),
        ((12, 0, 0), 8),
        ((12, 0, 4), 32),
        ((12, 4, 4), 128),
        ((13, 0, 0), -16),
        ((13, 0, 1), 64),
        ((13, 0, 5), -128),
        ((13, 1, 1), -256),
        ((13, 1, 5), 512),
        ((13, 5, 5), -1024),
    ]);

    Model {
        name: "fourfold",
        generators,
        grading_vector,
        q,
        nefpart,
        intnums,
        is_threefold: false,
    }
}

/// The invariants to compute and the coefficient type to compute them with.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Variant {
    GvRat,
    GwRat,
    GvFloat,
    GwFloat,
}

impl Variant {
    pub const ALL: [Variant; 4] = [
        Variant::GvRat,
        Variant::GwRat,
        Variant::GvFloat,
        Variant::GwFloat,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Variant::GvRat => "gv-rational",
            Variant::GwRat => "gw-rational",
            Variant::GvFloat => "gv-float",
            Variant::GwFloat => "gw-float",
        }
    }
}

/// A single benchmark case: one model, up to one degree, with one variant.
pub struct Scenario {
    pub model: &'static str,
    pub max_deg: u32,
    pub variant: Variant,
    /// Precision of the floating-point variants, in bits. Ignored by the exact
    /// ones.
    pub precision: u32,
    /// Number of samples criterion should collect for this case.
    pub sample_size: usize,
    /// Rough time of a single run, in seconds, measured on a 6-core machine.
    /// Only used to size criterion's measurement window; being off by a factor
    /// of a few costs time but does not affect what is reported.
    pub expected_secs: f64,
    /// Whether this case is only run when `CYGV_BENCH_HEAVY` is set.
    pub heavy: bool,
}

impl Scenario {
    /// Name of the criterion group this case belongs to.
    pub fn group(&self) -> String {
        format!("{}/deg{}", self.model, self.max_deg)
    }

    /// Fully qualified name, unique across all scenarios.
    pub fn id(&self) -> String {
        format!("{}/{}", self.group(), self.variant.name())
    }

    /// How long criterion should keep sampling this case.
    ///
    /// Criterion rounds the number of iterations per sample up, so the window is
    /// kept just under what one run per sample needs; asking for slightly more
    /// than that would double the time spent on the case. Criterion therefore
    /// suggests raising the target time on every slow case, which is expected.
    pub fn measurement_time(&self) -> std::time::Duration {
        let secs = 0.9 * self.expected_secs * self.sample_size as f64;
        std::time::Duration::from_secs_f64(secs.max(1.0))
    }

    /// Runs the case, returning the number of invariants that were computed.
    ///
    /// The model is passed in so that constructing it is not part of what is
    /// being measured; the entry points consume their arguments, so it still
    /// has to be cloned on every call.
    pub fn run(&self, m: &Model) -> usize {
        let gens = m.generators.clone();
        let grading = m.grading_vector.clone();
        let deg = Some(self.max_deg);
        let q = m.q.clone();
        let nefpart = m.nefpart.clone();
        let intnums = m.intnums.clone();
        let threads = n_threads();
        let prec = self.precision;

        match (self.variant, m.is_threefold) {
            (Variant::GvRat, true) => compute_gv_rat_threefold(
                gens, grading, deg, None, None, q, nefpart, intnums, threads, POOL_SIZE,
            )
            .len(),
            (Variant::GwRat, true) => compute_gw_rat_threefold(
                gens, grading, deg, None, None, q, nefpart, intnums, threads, POOL_SIZE,
            )
            .len(),
            (Variant::GvFloat, true) => compute_gv_float_threefold(
                gens, grading, deg, None, None, q, nefpart, intnums, threads, POOL_SIZE, prec,
            )
            .len(),
            (Variant::GwFloat, true) => compute_gw_float_threefold(
                gens, grading, deg, None, None, q, nefpart, intnums, threads, POOL_SIZE, prec,
            )
            .len(),
            (Variant::GvRat, false) => compute_gv_rat_nfold(
                gens, grading, deg, None, None, q, nefpart, intnums, threads, POOL_SIZE,
            )
            .len(),
            (Variant::GwRat, false) => compute_gw_rat_nfold(
                gens, grading, deg, None, None, q, nefpart, intnums, threads, POOL_SIZE,
            )
            .len(),
            (Variant::GvFloat, false) => compute_gv_float_nfold(
                gens, grading, deg, None, None, q, nefpart, intnums, threads, POOL_SIZE, prec,
            )
            .len(),
            (Variant::GwFloat, false) => compute_gw_float_nfold(
                gens, grading, deg, None, None, q, nefpart, intnums, threads, POOL_SIZE, prec,
            )
            .len(),
        }
    }
}

/// All benchmark cases, in the order they should be run.
///
/// The low-degree cases run fast enough for criterion to sample them properly.
/// The heavy ones are where the coefficients grow large enough for the bignum
/// arithmetic, the number pools and the sliding window of the series inversion
/// to dominate, which is what most changes to this crate are about; they are
/// skipped unless `CYGV_BENCH_HEAVY` is set, since sampling them takes minutes.
pub fn scenarios() -> Vec<Scenario> {
    // Per-run times are for one variant on a 6-core machine, in the order of
    // `Variant::ALL`. The `gw-float` n-fold case is the outlier: without an
    // exact zero to compare against it keeps every near-zero invariant, so it
    // computes two orders of magnitude more of them than its rational
    // counterpart.
    //
    //  model,       deg, prec, samples, per-run secs,             heavy
    let cases = [
        ("threefold", 10, 200, 100, [0.01, 0.01, 0.01, 0.01], false),
        ("threefold", 30, 500, 10, [0.76, 0.66, 0.75, 0.59], false),
        ("fourfold", 10, 200, 20, [0.09, 0.11, 0.08, 0.33], false),
        ("threefold", 50, 500, 10, [9.6, 7.8, 8.1, 6.2], true),
        ("fourfold", 15, 500, 10, [1.13, 1.50, 1.16, 15.07], true),
    ];

    let heavy_enabled = std::env::var_os("CYGV_BENCH_HEAVY").is_some();

    cases
        .into_iter()
        .filter(|(_, _, _, _, _, heavy)| heavy_enabled || !heavy)
        .flat_map(
            |(model, max_deg, precision, sample_size, expected_secs, heavy)| {
                Variant::ALL
                    .into_iter()
                    .zip(expected_secs)
                    .map(move |(variant, expected_secs)| Scenario {
                        model,
                        max_deg,
                        variant,
                        precision,
                        sample_size,
                        expected_secs,
                        heavy,
                    })
            },
        )
        .collect()
}

/// Returns the scenarios whose id contains one of the given filters. An empty
/// filter list matches everything.
pub fn filtered_scenarios(filters: &[String]) -> Vec<Scenario> {
    scenarios()
        .into_iter()
        .filter(|s| filters.is_empty() || filters.iter().any(|f| s.id().contains(f.as_str())))
        .collect()
}
