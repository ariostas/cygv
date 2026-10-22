//! Shared fixtures and scenario definitions for the `cygv` benchmarks.
//!
//! This module is compiled into every bench target, and not every target uses
//! every item, hence the blanket `dead_code` allowance.

#![allow(dead_code)]

use cygv::cone::ConeData;
use cygv::io::Input;
use cygv::{
    compute_gv, compute_gw, CYKind, Geometry, Method, MoriCone, Options, Problem, Selection,
};
use nalgebra::{dmatrix, DMatrix, RowDVector};
use std::collections::HashMap;

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
    pub nefpart: Vec<Vec<usize>>,
    pub intnums: HashMap<(usize, usize, usize), i32>,
    pub cy_kind: CYKind,
    /// The cone data, for the realistic models, with which the modular pipeline takes its
    /// candidates from the Hilbert bases of the cones $K_T$.
    pub cone_data: Option<ConeData>,
}

/// Looks up one of the benchmark models by name.
pub fn model(name: &str) -> Model {
    match name {
        "threefold" => threefold(),
        "fourfold" => fourfold(),
        "h11_8" => fixture(
            "h11_8",
            include_str!("../data/h11_8.yaml"),
            include_str!("../data/cones/h11_8.yaml"),
        ),
        "h11_9" => h11_9(),
        "h11_9_plike" => h11_9_plike(),
        "h11_10" => fixture(
            "h11_10",
            include_str!("../data/h11_10.yaml"),
            include_str!("../data/cones/h11_10.yaml"),
        ),
        "h11_11" => fixture(
            "h11_11",
            include_str!("../data/h11_11.yaml"),
            include_str!("../data/cones/h11_11.yaml"),
        ),
        _ => panic!("unknown model {name:?}"),
    }
}

/// Loads a model from one of the files in `benches/data`, which are ordinary
/// inputs to the command line interface, so each of them can also be run with
/// `cygv --file`. See the comments at the top of each file for where it comes
/// from.
///
/// These are the realistic models: CY threefold hypersurfaces with
/// $h^{1,1}$ between 8 and 11, whose generators include thousands of lattice
/// points of the Mori cone, the way CYTools passes them. Their cone data, from
/// `benches/data/cones`, is appended to them; it was computed with Normaliz
/// and is only complete up to the degree the scenarios reach (see the comments
/// at the top of those files). It is left out of the published crate.
fn fixture(name: &'static str, data: &str, cones: &str) -> Model {
    let load = |data: &str| {
        let mut inputs =
            Input::load_all(data).unwrap_or_else(|e| panic!("invalid fixture {name:?}: {e}"));
        assert_eq!(inputs.len(), 1, "fixture {name:?} must hold one document");
        inputs.remove(0)
    };
    let problem = load(data).problem;
    // The cone data replaces the generators, which come last.
    let head = &data[..data
        .find("\ngenerators:")
        .expect("the generators come last")
        + 1];
    let with_cones = load(&format!("{head}{cones}"));
    let MoriCone::Inequalities {
        hyperplanes,
        equations,
    } = with_cones.problem.cone
    else {
        panic!("the cone data of {name:?} has no inequalities")
    };
    let MoriCone::Generators(generators) = problem.cone else {
        panic!("fixture {name:?} has no generators")
    };
    let geometry = problem.geometry;
    Model {
        name,
        generators,
        grading_vector: RowDVector::from_vec(problem.grading_vector.unwrap()),
        cy_kind: CYKind::from_is_threefold(geometry.dimension() == 3),
        q: geometry.q,
        nefpart: geometry.nefpart,
        intnums: geometry.intnums,
        cone_data: Some(ConeData {
            hyperplanes,
            equations,
            cones: with_cones.options.hilbert_bases.unwrap(),
        }),
    }
}

/// A CY threefold hypersurface with $h^{1,1} = 9$ whose semigroup grows
/// quickly with the degree.
pub fn h11_9() -> Model {
    fixture(
        "h11_9",
        include_str!("../data/h11_9.yaml"),
        include_str!("../data/cones/h11_9.yaml"),
    )
}

/// The [`h11_9`] model under a grading vector with large entries of both signs.
///
/// The curve classes are the same, but their degrees are spread over hundreds
/// of values, so the series inversion walks through many more degree levels,
/// each holding only a few classes.
pub fn h11_9_plike() -> Model {
    Model {
        name: "h11_9_plike",
        grading_vector: RowDVector::from_row_slice(&[
            148, -121, -108, -29, -147, -114, -89, -141, -136,
        ]),
        // The cone data leaves out elements by their degree under the other grading.
        cone_data: None,
        ..h11_9()
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
        cy_kind: CYKind::Threefold,
        cone_data: None,
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
    let nefpart = vec![vec![0, 1, 2, 3], vec![4, 5, 6, 7, 8, 9, 10, 11]];

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
        cy_kind: CYKind::Nfold,
        cone_data: None,
    }
}

/// The invariants to compute.
///
/// The names keep the `-rational` suffix from when floating-point variants
/// existed, so that criterion baselines saved before then still match.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Variant {
    GvRat,
    GwRat,
}

impl Variant {
    pub const ALL: [Variant; 2] = [Variant::GvRat, Variant::GwRat];

    pub fn name(self) -> &'static str {
        match self {
            Variant::GvRat => "gv-rational",
            Variant::GwRat => "gw-rational",
        }
    }
}

/// The pipeline a scenario runs, and where the modular one takes its candidates from.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Pipeline {
    /// The reference pipeline.
    Reference,
    /// The modular pipeline, which computes the cone data from the generators.
    Modular,
    /// The modular pipeline, with the supplied cone data.
    ModularCones,
    /// The modular pipeline in its low-memory mode, which computes the cone data from the
    /// generators.
    ModularLowMemory,
}

impl Pipeline {
    /// The suffix of the scenario ids, empty for the reference so that the ids
    /// of baselines saved before the modular pipeline had scenarios still match.
    pub fn suffix(self) -> &'static str {
        match self {
            Pipeline::Reference => "",
            Pipeline::Modular => "/modular",
            Pipeline::ModularCones => "/modular-cones",
            Pipeline::ModularLowMemory => "/modular-lowmem",
        }
    }
}

/// How expensive a scenario is, which decides when it is run.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Tier {
    /// Always run.
    Quick,
    /// Run when `CYGV_BENCH_HEAVY` is set. Sampling each of these takes minutes.
    Heavy,
    /// Run when `CYGV_BENCH_HUGE` is set, and only by the memory benchmark: a
    /// single run takes several minutes and gigabytes of memory, which is too
    /// much for criterion to sample repeatedly.
    Huge,
}

/// A single benchmark case: one model, up to one degree, with one variant, in
/// one pipeline.
pub struct Scenario {
    pub model: &'static str,
    pub max_deg: u32,
    pub variant: Variant,
    pub pipeline: Pipeline,
    /// Number of samples criterion should collect for this case.
    pub sample_size: usize,
    /// Rough time of a single run, in seconds. Only used to size criterion's
    /// measurement window; being off by a factor of a few costs time but does
    /// not affect what is reported.
    pub expected_secs: f64,
    /// When this case is run.
    pub tier: Tier,
}

impl Scenario {
    /// Name of the criterion group this case belongs to.
    pub fn group(&self) -> String {
        format!("{}/deg{}", self.model, self.max_deg)
    }

    /// Fully qualified name, unique across all scenarios.
    pub fn id(&self) -> String {
        format!(
            "{}/{}{}",
            self.group(),
            self.variant.name(),
            self.pipeline.suffix()
        )
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
        let (method, cone, hilbert_bases) = match self.pipeline {
            Pipeline::Reference => (Method::Reference, None, None),
            Pipeline::Modular | Pipeline::ModularLowMemory => (Method::Modular, None, None),
            Pipeline::ModularCones => {
                let cd = m.cone_data.clone().expect("the model has no cone data");
                let cone = MoriCone::Inequalities {
                    hyperplanes: cd.hyperplanes,
                    equations: cd.equations,
                };
                (Method::Modular, Some(cone), Some(cd.cones))
            }
        };
        let problem = Problem {
            geometry: Geometry {
                q: m.q.clone(),
                intnums: m.intnums.clone(),
                nefpart: m.nefpart.clone(),
                ..Geometry::default()
            },
            cone: cone.unwrap_or_else(|| MoriCone::Generators(m.generators.clone())),
            selection: Selection::MaxDeg(self.max_deg),
            grading_vector: Some(m.grading_vector.iter().copied().collect()),
        };
        let options = Options {
            method,
            n_threads: n_threads().map(|n| n.max(1) as usize),
            low_memory: self.pipeline == Pipeline::ModularLowMemory,
            hilbert_bases,
            ..Options::default()
        };
        match self.variant {
            Variant::GvRat => compute_gv(&problem, &options).unwrap().entries.len(),
            Variant::GwRat => compute_gw(&problem, &options).unwrap().entries.len(),
        }
    }
}

/// All benchmark cases, in the order they should be run.
///
/// The quick cases run fast enough for criterion to sample them properly. The
/// heavy ones are where the coefficients grow large enough for the bignum
/// arithmetic and the sliding window of the series inversion to dominate,
/// which is what most changes to this crate are about. The huge ones are there
/// to track peak memory at the scale real computations reach.
pub fn scenarios() -> Vec<Scenario> {
    use Tier::{Heavy, Huge, Quick};

    // Per-run times are for one variant, in the order of `Variant::ALL`, and
    // `None` skips that variant. For the reference, the `threefold` and
    // `fourfold` times were measured on a 6-core machine and the rest on a
    // 16-core one; for the modular pipeline, all of them on 4 cores with
    // hyperthreading.
    //
    // The realistic models only run `gv-rational`, which is what they are used
    // for in practice; `gw-rational` would double their cost for little extra
    // coverage, since the toy models already exercise it.
    //
    // The modular pipeline is so much faster that its scenarios go to higher
    // degrees, mostly with the cone data, which is how it is meant to be used.
    // Its scenarios without the cone data also measure computing it from the
    // generators.
    let both = |t: [f64; 2]| t.map(Some);
    let gv_rat = |t: f64| [Some(t), None];
    let (r, m, mc, ml) = (
        Pipeline::Reference,
        Pipeline::Modular,
        Pipeline::ModularCones,
        Pipeline::ModularLowMemory,
    );
    //   model,          deg, pipeline, samples, per-run secs,         tier
    #[rustfmt::skip]
    let cases = [
        ("threefold",    10, r,     100, both([0.01, 0.01]),  Quick),
        ("threefold",    30, r,      10, both([0.76, 0.66]),  Quick),
        ("fourfold",     10, r,      20, both([0.09, 0.11]),  Quick),
        ("h11_9",        10, r,      10, gv_rat(1.1),         Quick),
        ("h11_9_plike", 200, r,      10, gv_rat(3.1),         Quick),
        ("h11_10",        8, r,      10, gv_rat(2.2),         Quick),
        ("h11_11",       10, r,      10, gv_rat(1.9),         Quick),
        ("fourfold",     15, m,      20, both([0.1, 0.1]),    Quick),
        ("h11_9",        12, m,      10, gv_rat(0.5),         Quick),
        ("h11_9_plike", 400, m,      10, gv_rat(0.56),        Quick),
        ("h11_8",        16, mc,     10, gv_rat(0.55),        Quick),
        ("h11_9",        16, mc,     10, gv_rat(0.5),         Quick),
        ("h11_10",       12, mc,     10, gv_rat(0.45),        Quick),
        ("h11_11",       16, mc,     10, gv_rat(0.7),         Quick),
        ("threefold",    50, r,      10, both([9.6, 7.8]),    Heavy),
        ("fourfold",     15, r,      10, both([1.13, 1.50]),  Heavy),
        ("h11_8",        14, r,      10, gv_rat(11.0),        Heavy),
        ("h11_9",        14, r,      10, gv_rat(23.0),        Heavy),
        ("h11_9_plike", 240, r,      10, gv_rat(14.0),        Heavy),
        ("h11_10",       10, r,      10, gv_rat(8.5),         Heavy),
        ("h11_11",       12, r,      10, gv_rat(12.0),        Heavy),
        ("h11_9",        14, m,      10, gv_rat(1.8),         Heavy),
        ("h11_9_plike", 500, m,      10, gv_rat(3.7),         Heavy),
        ("h11_8",        20, mc,     10, gv_rat(9.5),         Heavy),
        ("h11_9",        20, mc,     10, gv_rat(5.0),         Heavy),
        ("h11_10",       16, mc,     10, gv_rat(10.0),        Heavy),
        ("h11_11",       18, mc,     10, gv_rat(1.8),         Heavy),
        ("h11_8",        22, m,      10, gv_rat(4.2),         Heavy),
        ("h11_8",        22, ml,     10, gv_rat(11.6),        Heavy),
        ("h11_8",        20, r,       1, gv_rat(1600.0),      Huge),
        ("h11_9",        18, r,       1, gv_rat(800.0),       Huge),
        ("h11_9_plike", 300, r,       1, gv_rat(123.0),       Huge),
        ("h11_11",       16, r,       1, gv_rat(240.0),       Huge),
        ("h11_10",       20, m,       1, gv_rat(25.0),        Huge),
        ("h11_10",       20, ml,      1, gv_rat(65.0),        Huge),
    ];

    let heavy_enabled = std::env::var_os("CYGV_BENCH_HEAVY").is_some();
    let huge_enabled = std::env::var_os("CYGV_BENCH_HUGE").is_some();

    cases
        .into_iter()
        .filter(|&(_, _, _, _, _, tier)| match tier {
            Quick => true,
            Heavy => heavy_enabled,
            Huge => huge_enabled,
        })
        .flat_map(
            |(model, max_deg, pipeline, sample_size, expected_secs, tier)| {
                Variant::ALL.into_iter().zip(expected_secs).filter_map(
                    move |(variant, expected_secs)| {
                        Some(Scenario {
                            model,
                            max_deg,
                            variant,
                            pipeline,
                            sample_size,
                            expected_secs: expected_secs?,
                            tier,
                        })
                    },
                )
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
