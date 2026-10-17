//! A stress test of the parallel code of the modular pipeline: repeated runs at several thread
//! counts, with the extraction along trees of two and three rows, and curve by curve, its
//! layers applied curve by curve and level-synchronously, must all give the results of a run
//! on one thread.
//!
//! It runs each case 40 times, which takes half a minute in a release build and much longer
//! in a debug one, so it is ignored by default; run it before releases and after changes to
//! the parallel code with `cargo test --release --test stress -- --ignored`.

#![cfg(feature = "cli")]

use cygv::cone::ConeData;
use cygv::io::Input;
use cygv::modular::extraction::Tuning;
use cygv::modular::{compute_gv, CandidateSource, CurveInvariants, Geometry, Options};
use cygv::MoriCone;
use nalgebra::{DMatrix, DVector};

const H11_9: &str = include_str!("../benches/data/h11_9.yaml");
const H11_9_CONES: &str = include_str!("data/h11_9_cones.yaml");
const FOURFOLD: &str = include_str!("../examples/fourfold.yaml");

/// What the modular pipeline needs of an input of the command line interface, with the cone
/// data appended to it if given.
struct Case {
    geometry: Geometry,
    grading_vector: Vec<i32>,
    generators: DMatrix<i32>,
    cone_data: Option<ConeData>,
}

fn case(yaml: &str, cones: Option<&str>) -> Case {
    let input = Input::load_all(yaml).unwrap().remove(0);
    let problem = input.problem;
    let w = problem.grading_vector.unwrap();
    let g = &problem.geometry;
    let nefpart: Vec<DVector<i32>> = g
        .nefpart
        .iter()
        .map(|p| DVector::from_iterator(p.len(), p.iter().map(|&r| r as i32)))
        .collect();
    let kind = cygv::CYKind::from_is_threefold(g.dimension() == 3);
    let geometry = Geometry::from_toric(&g.q, &nefpart, &g.intnums, kind, &w).unwrap();
    let MoriCone::Generators(generators) = problem.cone else {
        panic!("no generators")
    };
    // The cone data replaces the generators of the input.
    let cone_data = cones.map(|cones| {
        let head = &yaml[..yaml.find("\ngenerators:").unwrap() + 1];
        let input = Input::load_all(&format!("{head}{cones}"))
            .unwrap()
            .remove(0);
        let MoriCone::Inequalities {
            hyperplanes,
            equations,
        } = input.problem.cone
        else {
            panic!("no inequalities")
        };
        ConeData {
            hyperplanes,
            equations,
            cones: input.options.hilbert_bases.unwrap(),
        }
    });
    Case {
        geometry,
        grading_vector: w,
        generators,
        cone_data,
    }
}

fn run(case: &Case, max_deg: u32, options: &Options) -> Vec<CurveInvariants> {
    let w = &case.grading_vector;
    let gens: Vec<&[i32]> = case.generators.as_slice().chunks(w.len()).collect();
    let source = match &case.cone_data {
        Some(cone_data) => CandidateSource::Cones(cone_data),
        None => CandidateSource::Generators(&gens),
    };
    compute_gv(source, w, max_deg, &[], &case.geometry, options).unwrap()
}

#[test]
#[ignore = "slow; run before releases and after changes to the parallel code"]
fn test_stress() {
    let cases = [
        ("h11_9, cones", case(H11_9, Some(H11_9_CONES)), 10),
        ("h11_9, generators", case(H11_9, None), 12),
        ("fourfold", case(FOURFOLD, None), 15),
    ];
    let tree = |depth| Tuning {
        tree_depth: Some(depth),
        ..Tuning::default()
    };
    let by_level = Tuning {
        sync_max_curves: usize::MAX,
        sync_min_terms: 0,
        ..tree(0)
    };
    for (name, case, max_deg) in &cases {
        let one = Options {
            threads: Some(1),
            ..Options::default()
        };
        let expected = run(case, *max_deg, &one);
        assert!(!expected.is_empty());
        for threads in [2, 3, 4, 8] {
            for extraction in [tree(2), tree(3), tree(0), by_level.clone()] {
                let options = Options {
                    threads: Some(threads),
                    tuning: cygv::modular::Tuning {
                        extraction,
                        ..Default::default()
                    },
                    ..Options::default()
                };
                for _ in 0..5 {
                    let got = run(case, *max_deg, &options);
                    assert!(
                        got == expected,
                        "{name}: {options:?} differs from one thread"
                    );
                }
            }
        }
    }
}
