//! A stress test of the parallel code of the modular pipeline: repeated runs at several thread
//! counts, with the extraction's layers applied curve by curve and level-synchronously, must
//! all give the results of a run on one thread.
//!
//! It runs each case 40 times, which takes half a minute in a release build and much longer
//! in a debug one, so it is ignored by default; run it before releases and after changes to
//! the parallel code with `cargo test --release --test stress -- --ignored`.

#![cfg(feature = "cli")]

use cygv::io::Input;
use cygv::modular::extraction::Tuning;
use cygv::modular::{compute_gv, CandidateSource, CurveInvariants, Geometry, Options};

const H11_9: &str = include_str!("../benches/data/h11_9.yaml");
const H11_9_CONES: &str = include_str!("data/h11_9_cones.yaml");
const FOURFOLD: &str = include_str!("../examples/fourfold.yaml");

fn run(input: &Input, max_deg: u32, cones: bool, options: &Options) -> Vec<CurveInvariants> {
    let w: Vec<i32> = input.grading_vector.iter().copied().collect();
    let geometry =
        Geometry::from_toric(&input.q, &input.nefpart, &input.intnums, input.cy_kind, &w).unwrap();
    let h11 = w.len();
    let gens: Vec<&[i32]> = input.generators.as_slice().chunks(h11).collect();
    let source = match (&input.cone_data, cones) {
        (Some(cone_data), true) => CandidateSource::Cones(cone_data),
        _ => CandidateSource::Generators(&gens),
    };
    compute_gv(source, &w, max_deg, &[], &geometry, options).unwrap()
}

#[test]
#[ignore = "slow; run before releases and after changes to the parallel code"]
fn test_stress() {
    let h11_9 = Input::load_all(&format!("{H11_9}\n{H11_9_CONES}"))
        .unwrap()
        .remove(0);
    let fourfold = Input::load_all(FOURFOLD).unwrap().remove(0);
    let cases = [
        ("h11_9, cones", &h11_9, 10, true),
        ("h11_9, generators", &h11_9, 12, false),
        ("fourfold", &fourfold, 15, false),
    ];
    let by_level = Tuning {
        sync_max_curves: usize::MAX,
        sync_min_terms: 0,
    };
    for (name, input, max_deg, cones) in cases {
        let one = Options {
            threads: Some(1),
            ..Options::default()
        };
        let expected = run(input, max_deg, cones, &one);
        assert!(!expected.is_empty());
        for threads in [2, 3, 4, 8] {
            for extraction in [Tuning::default(), by_level.clone()] {
                let options = Options {
                    threads: Some(threads),
                    extraction,
                    ..Options::default()
                };
                for _ in 0..5 {
                    let got = run(input, max_deg, cones, &options);
                    assert!(
                        got == expected,
                        "{name}: {options:?} differs from one thread"
                    );
                }
            }
        }
    }
}
