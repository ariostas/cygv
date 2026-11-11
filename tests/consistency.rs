//! Inconsistent input must be an error rather than wrong invariants, with both methods and at
//! low degree, where the invariants of the modular pipeline are often still integers (the
//! cases of the review of PR #94). The geometry is the resolved P(1,1,2,2,2)[8] of
//! `examples/cones.yaml`, its vex flop of `examples/vex.yaml` and the fourfold of
//! `examples/fourfold.yaml`. What no check can see is tested as such at the end.

#![cfg(feature = "cli")]

use cygv::io::Input;
use cygv::{modular, reference, Computed, Consistency, Error, Method, Verify, MAX_DEGREE};

/// The resolved P(1,1,2,2,2)[8] with the given intersection numbers, Mori cone and max_deg.
fn p11222(intnums: &str, cone: &str, max_deg: u32) -> String {
    format!(
        "q: [[0, 0, 1, 1, 1, 1], [1, 1, 0, 0, 0, -2]]\nintnums: {intnums}\n{cone}\n\
         grading_vector: [1, 1]\nmax_deg: {max_deg}\n"
    )
}

const INTNUMS: &str = "[[0, 0, 0, 8], [0, 0, 1, 4]]";
const GENERATORS: &str = "generators: [[1, 0], [0, 1]]";

/// The result of the document with the method and the consistency check given.
fn run(doc: &str, method: Method, verify: Verify) -> Result<Computed, Error> {
    // Parsed without the checks, which the computation makes too.
    let mut input = Input::parse_all(doc, None)
        .unwrap_or_else(|e| panic!("{e}"))
        .remove(0);
    input.options.method = method;
    input.options.verify = verify;
    input.compute()
}

/// Whether an error is the consistency check of either method, or a GV invariant that is not
/// an integer, which inconsistent input also gives.
fn inconsistent(e: &Error) -> bool {
    matches!(
        e,
        Error::Modular(modular::Error::Inconsistent { .. } | modular::Error::NotInteger { .. })
            | Error::Reference(reference::Error::SeriesInversion(
                reference::SeriesInversionError::Inconsistent
                    | reference::SeriesInversionError::NonIntegerGVError
            ))
    )
}

/// The field of an input error.
fn field(res: Result<Computed, Error>) -> &'static str {
    match res {
        Err(Error::Input { field, .. }) => field,
        Err(e) => panic!("not an input error: {e}"),
        Ok(_) => panic!("not an error"),
    }
}

/// The GV invariants as (curve, value) strings.
fn gv(res: Result<Computed, Error>) -> Vec<(Vec<i32>, String)> {
    match res.unwrap() {
        Computed::Gv(gv) => gv
            .entries
            .into_iter()
            .map(|e| (e.curve, e.value.to_string()))
            .collect(),
        Computed::Gw(_) => panic!("GW invariants instead of GV"),
    }
}

#[test]
fn test_wrong_intnums() {
    for intnums in [
        // A spurious intersection number.
        "[[0, 0, 0, 8], [0, 0, 1, 4], [0, 1, 1, 1]]",
        "[[0, 0, 0, 8], [0, 0, 1, 4], [1, 1, 1, 2]]",
        // A wrong one.
        "[[0, 0, 0, 8], [0, 0, 1, 2]]",
        "[[0, 0, 0, 4], [0, 0, 1, 4]]",
    ] {
        for max_deg in [1, 2] {
            let doc = p11222(intnums, GENERATORS, max_deg);
            for method in [Method::Modular, Method::Reference] {
                let e = run(&doc, method, Verify::Quick).unwrap_err();
                assert!(inconsistent(&e), "{intnums} at {max_deg}, {method:?}: {e}");
            }
            // Without the check the modular pipeline gives integers, which are wrong.
            assert!(run(&doc, Method::Modular, Verify::Off).is_ok(), "{intnums}");
        }
    }
}

#[test]
fn test_cone_too_small() {
    // A cone that misses the class (1, 0), which meets every divisor nonnegatively, is refused
    // before computing.
    let doc = p11222(INTNUMS, "generators: [[1, 1], [0, 1]]", 2);
    for method in [Method::Modular, Method::Reference] {
        assert_eq!(field(run(&doc, method, Verify::Off)), "generators");
    }
    // One that only misses (0, 1) passes that, and is caught at degree 2, where it matters.
    let doc = p11222(INTNUMS, "generators: [[1, 0], [1, 1]]", 2);
    for method in [Method::Modular, Method::Reference] {
        let e = run(&doc, method, Verify::Quick).unwrap_err();
        assert!(inconsistent(&e), "{method:?}: {e}");
    }
}

#[test]
fn test_incomplete_hilbert_bases() {
    // The cones of divisor 5 with (0, 1) replaced by (1, 1), as in examples/cones.yaml
    // otherwise.
    let mut bases = String::from("mori_hyperplanes: [[1, 0], [0, 1]]\nhilbert_bases:\n");
    for t in (0..6)
        .map(|r| vec![r])
        .chain((0..6).flat_map(|r| (r + 1..6).map(move |s| vec![r, s])))
    {
        let basis = if t.contains(&5) {
            "[[1, 0], [1, 1]]"
        } else {
            "[[1, 0], [2, 1]]"
        };
        bases.push_str(&format!(
            "  - {{divisors: {t:?}, hilbert_basis: {basis}}}\n"
        ));
    }
    let doc = p11222(INTNUMS, &bases, 2);
    assert_eq!(
        field(run(&doc, Method::Modular, Verify::Off)),
        "hilbert_bases"
    );
    // The reference does not use them.
    let truth = gv(run(
        &p11222(INTNUMS, GENERATORS, 2),
        Method::Reference,
        Verify::Off,
    ));
    assert_eq!(gv(run(&doc, Method::Reference, Verify::Off)), truth);
}

#[test]
fn test_vex_without_strata() {
    let vex = "q: [[1, 0, -1, 1, 0, 1], [0, 1, 1, 1, 1, 0]]\n\
               intnums: [[0, 1, 1, 4], [1, 1, 1, 2]]\ngrading_vector: [2, 1]\nmax_deg: 2\n";
    for cone in [
        "generators: [[0, 1], [1, -1]]",
        "mori_hyperplanes: [[1, 1], [1, 0]]",
    ] {
        let doc = format!("{vex}{cone}\n");
        for method in [Method::Modular, Method::Reference] {
            assert_eq!(field(run(&doc, method, Verify::Off)), "fan", "{cone}");
        }
        let doc = format!("{doc}strata: [{{divisors: [1, 2, 4], class: [1, -1]}}]\n");
        for method in [Method::Modular, Method::Reference] {
            assert!(run(&doc, method, Verify::Quick).is_ok(), "{cone}");
        }
    }
}

const FOURFOLD: &str = include_str!("../examples/fourfold.yaml");

#[test]
fn test_dimension() {
    // A complete intersection without its nef partition is a hypersurface of one dimension
    // more, which only the dimension, when given, catches.
    let nefpart = "nefpart:\n  - [0, 1, 2, 3]\n  - [4, 5, 6, 7, 8, 9, 10, 11]\n";
    assert!(FOURFOLD.contains(nefpart));
    let without = FOURFOLD
        .replace(nefpart, "")
        .replace("min_points: 100", "max_deg: 1");
    match run(&without, Method::Modular, Verify::Quick).unwrap() {
        Computed::Gv(gv) => assert_eq!(gv.dimension, 5),
        Computed::Gw(_) => panic!("GW invariants instead of GV"),
    }
    assert_eq!(
        field(run(
            &format!("{without}dimension: 4\n"),
            Method::Modular,
            Verify::Quick
        )),
        "dimension"
    );
    let with = FOURFOLD.replace("min_points: 100", "max_deg: 1");
    assert!(run(
        &format!("{with}dimension: 4\n"),
        Method::Modular,
        Verify::Quick
    )
    .is_ok());
}

#[test]
fn test_nfold_intnums() {
    // An n-fold has no consistency check: a changed intersection number only shows when an
    // invariant is not an integer, here at degree 4, while degree 1 gives a wrong integer.
    assert!(FOURFOLD.contains("  - [0, 1, 1, -8]\n"));
    let wrong = FOURFOLD.replace("  - [0, 1, 1, -8]\n", "  - [0, 1, 1, -7]\n");
    let low = run(
        &wrong.replace("min_points: 100", "max_deg: 1"),
        Method::Modular,
        Verify::Full,
    )
    .unwrap();
    let Computed::Gv(low) = low else {
        panic!("GW invariants instead of GV")
    };
    assert_eq!(low.report.consistency, Consistency::Unchecked);
    let e = run(
        &wrong.replace("min_points: 100", "max_deg: 4"),
        Method::Modular,
        Verify::Quick,
    )
    .unwrap_err();
    assert!(inconsistent(&e), "{e}");
}

#[test]
fn test_scaling() {
    // Intersection numbers scaled by a common factor keep the gradient condition, so neither
    // method can see them: every invariant is scaled instead. Only a fan, which the
    // intersection numbers are checked against, catches them.
    let truth = gv(run(
        &p11222(INTNUMS, GENERATORS, 3),
        Method::Modular,
        Verify::Full,
    ));
    let doubled = p11222("[[0, 0, 0, 16], [0, 0, 1, 8]]", GENERATORS, 3);
    for method in [Method::Modular, Method::Reference] {
        let got = gv(run(&doubled, method, Verify::Full));
        let scaled: Vec<_> = truth
            .iter()
            .map(|(c, v)| (c.clone(), (2 * v.parse::<i64>().unwrap()).to_string()))
            .collect();
        assert_eq!(got, scaled, "{method:?}");
    }
}

#[test]
fn test_degree_limit() {
    let doc = p11222(INTNUMS, GENERATORS, MAX_DEGREE + 1);
    assert_eq!(field(run(&doc, Method::Auto, Verify::Quick)), "max_deg");
    let doc = p11222(INTNUMS, GENERATORS, 2)
        .replace("max_deg: 2", "target_points: [[1000000000, 1000000000]]");
    for method in [Method::Modular, Method::Reference] {
        assert_eq!(field(run(&doc, method, Verify::Quick)), "target_points");
    }
}

#[test]
fn test_report() {
    let doc = p11222(INTNUMS, GENERATORS, 3);
    for (method, verify, expected) in [
        (Method::Modular, Verify::Off, Consistency::Unchecked),
        (Method::Modular, Verify::Quick, Consistency::Quick),
        (Method::Modular, Verify::Full, Consistency::Full),
        (Method::Reference, Verify::Off, Consistency::Full),
    ] {
        let Computed::Gv(gv) = run(&doc, method, verify).unwrap() else {
            panic!("GW invariants instead of GV")
        };
        assert_eq!(gv.report.consistency, expected, "{method:?}, {verify:?}");
    }
}
