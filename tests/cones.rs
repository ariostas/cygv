//! Tests of the inequalities of the Mori cone and of the Hilbert bases of the cones $K_T$, from
//! which the modular pipeline takes its candidates.
//!
//! The data in `tests/data` was computed with Normaliz for `benches/data/h11_9.yaml` (see the
//! headers of the files), and is appended to its input.

#![cfg(feature = "cli")]

use cygv::cone::ConeData;
use cygv::io::Input;
use cygv::{
    compute, CandidateCone, Computed, Error, Geometry, InvariantKind, Method, MoriCone, Options,
    Problem, Selection,
};
use nalgebra::{dmatrix, DMatrix};
use std::collections::HashMap;

const H11_8: &str = include_str!("../benches/data/h11_8.yaml");
const H11_9: &str = include_str!("../benches/data/h11_9.yaml");
const H11_10: &str = include_str!("../benches/data/h11_10.yaml");
const H11_11: &str = include_str!("../benches/data/h11_11.yaml");
const H11_9_CONES: &str = include_str!("data/h11_9_cones.yaml");
const H11_9_FACET: &str = include_str!("data/h11_9_facet.yaml");

/// The text without the top-level fields `keys`, each with the lines below it.
fn without(text: &str, keys: &[&str]) -> String {
    let mut res = String::new();
    let mut skipping = false;
    for line in text.lines() {
        if !line.starts_with([' ', '-']) {
            skipping = keys
                .iter()
                .any(|k| line.strip_prefix(k).is_some_and(|r| r.starts_with(':')));
        }
        if !skipping {
            res += line;
            res.push('\n');
        }
    }
    res
}

/// The input of a bench file at `max_deg`, with `extra` appended: with the generators of the
/// bench file as the Mori cone if `generators`, and the inequalities of `extra` otherwise.
fn input(bench: &str, generators: bool, max_deg: u32, extra: &str) -> Input {
    let bench = without(bench, &["max_deg"]);
    let (extra, bench) = if generators {
        (
            without(extra, &["generators", "mori_hyperplanes", "mori_equations"]),
            bench,
        )
    } else {
        (
            without(extra, &["generators"]),
            without(&bench, &["generators"]),
        )
    };
    let mut inputs = Input::load_all(&format!("{bench}max_deg: {max_deg}\n{extra}")).unwrap();
    assert_eq!(inputs.len(), 1);
    inputs.remove(0)
}

/// The generators of a bench file, or of the facet's file.
fn generators(bench: &str) -> DMatrix<i32> {
    let base = if bench == H11_9_FACET {
        format!(
            "{}{}",
            without(H11_9, &["generators"]),
            without(bench, &["mori_hyperplanes", "mori_equations"])
        )
    } else {
        bench.to_owned()
    };
    match Input::load_all(&base).unwrap().remove(0).problem.cone {
        MoriCone::Generators(g) => g,
        cone => panic!("{cone:?}"),
    }
}

/// The inequalities and the Hilbert bases of an input.
fn cone_data(input: &Input) -> ConeData {
    let MoriCone::Inequalities {
        hyperplanes,
        equations,
    } = input.problem.cone.clone()
    else {
        panic!("no inequalities")
    };
    ConeData {
        hyperplanes,
        equations,
        cones: input.options.hilbert_bases.clone().unwrap(),
    }
}

fn grading_vector(input: &Input) -> Vec<i32> {
    input.problem.grading_vector.clone().unwrap()
}

type Results = HashMap<(Vec<i32>, usize), String>;

fn results(input: &Input, method: Method) -> Result<Results, Error> {
    let options = Options {
        method,
        ..input.options.clone()
    };
    let label = |c: Vec<i32>, s: Option<usize>| (c, s.unwrap_or(0));
    Ok(match compute(&input.problem, input.kind, &options)? {
        Computed::Gv(gv) => gv
            .entries
            .into_iter()
            .map(|e| (label(e.curve, e.surface), e.value.to_string()))
            .collect(),
        Computed::Gw(gw) => gw
            .entries
            .into_iter()
            .map(|e| (label(e.curve, e.surface), e.value.to_string()))
            .collect(),
    })
}

/// The same input without the Hilbert bases.
fn without_bases(input: &Input) -> Input {
    let mut input = input.clone();
    input.options.hilbert_bases = None;
    input
}

#[test]
fn test_h11_9() {
    let input = input(H11_9, true, 8, H11_9_CONES);
    let expected = results(&input, Method::Reference).unwrap();
    assert_eq!(expected.len(), 551);
    assert_eq!(results(&input, Method::Modular).unwrap(), expected);
    assert_eq!(results(&input, Method::Auto).unwrap(), expected);
    assert_eq!(
        results(&without_bases(&input), Method::Modular).unwrap(),
        expected
    );

    // Without the generators both methods take the Mori cone from its inequalities.
    let input = self::input(H11_9, false, 8, H11_9_CONES);
    assert!(matches!(input.problem.cone, MoriCone::Inequalities { .. }));
    assert_eq!(results(&input, Method::Auto).unwrap(), expected);
    assert_eq!(results(&input, Method::Reference).unwrap(), expected);
    assert_eq!(
        results(&without_bases(&input), Method::Auto).unwrap(),
        expected
    );
}

/// Without an element of the Hilbert basis of the Mori cone among the generators, their
/// semigroup has a hole, over which the invariants would not be integers. Both methods use
/// every lattice point of the cone the generators span instead, with or without the Hilbert
/// bases.
#[test]
fn test_saturation() {
    let missing = "[1,0,0,0,0,0,1,0,0]";
    let bench: String = H11_9
        .lines()
        .filter(|l| l.trim().trim_start_matches("- ") != missing)
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(bench.lines().count(), H11_9.lines().count() - 2);
    let input = input(&bench, true, 8, H11_9_CONES);
    let full = self::input(H11_9, true, 8, H11_9_CONES);
    let expected = results(&full, Method::Modular).unwrap();
    assert_eq!(expected[&(vec![1, 0, 0, 0, 0, 0, 1, 0, 0], 0)], "20");
    assert_eq!(results(&input, Method::Modular).unwrap(), expected);
    for method in [Method::Modular, Method::Reference] {
        assert_eq!(results(&without_bases(&input), method).unwrap(), expected);
    }
}

/// The cones of a facet give the invariants of the curves on it.
#[test]
fn test_facet() {
    let input = input(H11_9, false, 10, H11_9_FACET);
    assert_eq!(cone_data(&input).equations.shape(), (9, 1));
    let expected = results(&input, Method::Reference).unwrap();
    assert!(expected.len() > 100);
    let facet = [0, 2, 2, 2, 2, 2, 0, 0, 1];
    assert!(expected
        .keys()
        .all(|(c, _)| { c.iter().zip(facet).map(|(x, f)| x * f).sum::<i32>() == 0 }));
    assert_eq!(results(&input, Method::Modular).unwrap(), expected);
    assert_eq!(
        results(&without_bases(&input), Method::Modular).unwrap(),
        expected
    );
}

/// The columns of a matrix, sorted.
fn sorted_columns(m: &DMatrix<i32>) -> Vec<Vec<i32>> {
    let mut v: Vec<Vec<i32>> = m
        .column_iter()
        .map(|c| c.iter().copied().collect())
        .collect();
    v.sort();
    v
}

/// The Hilbert bases that `cone::cone_data` computes from the hyperplanes are Normaliz's, up to
/// the degree of the computation, on the four benchmark geometries and on the facet of `h11_9`.
#[test]
fn test_cone_data() {
    let cases = [
        (H11_8, include_str!("../benches/data/cones/h11_8.yaml"), 12),
        (H11_9, include_str!("../benches/data/cones/h11_9.yaml"), 12),
        (
            H11_10,
            include_str!("../benches/data/cones/h11_10.yaml"),
            10,
        ),
        (
            H11_11,
            include_str!("../benches/data/cones/h11_11.yaml"),
            12,
        ),
        (H11_9, H11_9_FACET, 10),
    ];
    for (bench, cones, max_deg) in cases {
        let input = input(bench, false, max_deg, cones);
        let normaliz = cone_data(&input);
        let w = grading_vector(&input);
        let found = cygv::cone::cone_data(
            &normaliz.hyperplanes,
            &normaliz.equations,
            &input.problem.geometry.q,
            &[],
            &w,
            max_deg,
        )
        .unwrap();
        // Each distinct cone once (see `cone::cone_data`), so fewer than Normaliz's.
        assert!(found.cones.len() < normaliz.cones.len());
        let degree = |c: &Vec<i32>| -> i32 { c.iter().zip(&w).map(|(x, y)| x * y).sum() };
        let truncated = |c: &CandidateCone| {
            let mut basis = sorted_columns(&c.hilbert_basis);
            basis.retain(|c| degree(c) <= max_deg as i32);
            basis
        };
        for cone in &found.cones {
            let theirs = normaliz
                .cones
                .iter()
                .find(|c| c.divisors == cone.divisors)
                .unwrap();
            assert_eq!(
                sorted_columns(&cone.hilbert_basis),
                truncated(theirs),
                "{:?}",
                cone.divisors
            );
        }
        // And every one of Normaliz's with a divisor negative on it is one of them.
        let q = &input.problem.geometry.q;
        let found: Vec<Vec<Vec<i32>>> = found
            .cones
            .iter()
            .map(|c| sorted_columns(&c.hilbert_basis))
            .collect();
        for theirs in &normaliz.cones {
            let basis = truncated(theirs);
            let negative = theirs.divisors.iter().any(|&r| {
                basis.iter().any(|c| {
                    c.iter()
                        .enumerate()
                        .map(|(t, &x)| x * q[(r, t)])
                        .sum::<i32>()
                        < 0
                })
            });
            assert_eq!(found.contains(&basis), negative, "{:?}", theirs.divisors);
        }
    }
}

/// The facets that `cone::facets` finds from the generators are Normaliz's, on the Mori cones
/// of the four benchmark geometries and on the facet of `h11_9`.
#[test]
fn test_facets() {
    let benches = [
        (H11_8, include_str!("../benches/data/cones/h11_8.yaml")),
        (H11_9, include_str!("../benches/data/cones/h11_9.yaml")),
        (H11_10, include_str!("../benches/data/cones/h11_10.yaml")),
        (H11_11, include_str!("../benches/data/cones/h11_11.yaml")),
    ];
    for (bench, cones) in benches {
        let normaliz = cone_data(&input(bench, false, 8, cones));
        let found = cygv::cone::facets(&generators(bench)).unwrap();
        assert_eq!(found.equations.ncols(), 0);
        assert_eq!(
            sorted_columns(&found.hyperplanes),
            sorted_columns(&normaliz.hyperplanes)
        );
    }

    // On the facet, a hyperplane is only defined up to multiples of the equation, so
    // Normaliz's are moved along it to vanish where `facets` makes them vanish, at the
    // equation's first nonzero entry.
    let normaliz = cone_data(&input(H11_9, false, 10, H11_9_FACET));
    let found = cygv::cone::facets(&generators(H11_9_FACET)).unwrap();
    assert_eq!(found.equations, normaliz.equations);
    let e: Vec<i32> = normaliz.equations.iter().copied().collect();
    let p = e.iter().position(|&x| x != 0).unwrap();
    let reduced: Vec<i32> = normaliz
        .hyperplanes
        .column_iter()
        .flat_map(|h| {
            let v: Vec<i32> = h.iter().zip(&e).map(|(x, y)| e[p] * x - h[p] * y).collect();
            let g = v.iter().fold(0, |g, &x| gcd(g, x.abs()));
            v.into_iter().map(move |x| x / g)
        })
        .collect();
    let reduced = DMatrix::from_vec(9, normaliz.hyperplanes.ncols(), reduced);
    assert_eq!(sorted_columns(&found.hyperplanes), sorted_columns(&reduced));
}

fn gcd(a: i32, b: i32) -> i32 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

/// The curve classes of `results` below one of the columns of `targets`, i.e. those that leave
/// it in the cone of `cone_data` when subtracted from it.
fn below(results: &Results, cone_data: &ConeData, targets: &DMatrix<i32>) -> Results {
    let dot =
        |h: &[i32], c: &[i32]| -> i64 { h.iter().zip(c).map(|(&x, &y)| x as i64 * y as i64).sum() };
    let h11 = cone_data.hyperplanes.nrows();
    let in_cone = |c: &[i32]| {
        cone_data
            .hyperplanes
            .as_slice()
            .chunks(h11)
            .all(|h| dot(h, c) >= 0)
            && cone_data
                .equations
                .as_slice()
                .chunks(h11)
                .all(|e| dot(e, c) == 0)
    };
    results
        .iter()
        .filter(|((c, _), _)| {
            targets.column_iter().any(|p| {
                let d: Vec<i32> = p.iter().zip(c).map(|(x, y)| x - y).collect();
                in_cone(&d)
            })
        })
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect()
}

/// The curves of the two highest degrees in `results`, the last of each.
fn two_targets(input: &Input, results: &Results) -> DMatrix<i32> {
    let w = grading_vector(input);
    let degree = |c: &[i32]| -> i32 { w.iter().zip(c).map(|(w, x)| w * x).sum() };
    let mut curves: Vec<&Vec<i32>> = results.keys().map(|(c, _)| c).collect();
    curves.sort_unstable_by_key(|c| (degree(c), c.to_vec()));
    let top = degree(curves[curves.len() - 1]);
    let lower = curves.iter().rev().find(|c| degree(c) < top).unwrap();
    DMatrix::from_iterator(
        w.len(),
        2,
        curves[curves.len() - 1].iter().chain(lower.iter()).copied(),
    )
}

/// The invariants below the targets, with and without the generators, are those of the full
/// computation there, the same as from the semigroup of the generators, which has every lattice
/// point of the Mori cone in these degrees.
#[test]
fn test_target_points() {
    let full = input(H11_9, true, 8, H11_9_CONES);
    let expected_full = results(&full, Method::Modular).unwrap();
    let cone_data = cone_data(&input(H11_9, false, 8, H11_9_CONES));
    let targets = two_targets(&full, &expected_full);
    let expected = below(&expected_full, &cone_data, &targets);
    assert!(expected.len() < expected_full.len() / 2);
    let w = grading_vector(&full);
    for generators in [true, false] {
        let mut input = input(H11_9, generators, 8, H11_9_CONES);
        input.problem.selection = Selection::TargetPoints {
            points: targets.clone(),
            max_deg: None,
        };
        assert_eq!(results(&input, Method::Auto).unwrap(), expected);
        assert_eq!(results(&input, Method::Modular).unwrap(), expected);
        assert_eq!(
            results(&without_bases(&input), Method::Modular).unwrap(),
            expected
        );
        // With a maximum degree below that of the first target.
        input.set_max_deg(7).unwrap();
        let lower: Results = expected
            .iter()
            .filter(|((c, _), _)| c.iter().zip(&w).map(|(x, w)| x * w).sum::<i32>() <= 7)
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        assert!(lower.len() < expected.len());
        assert_eq!(results(&input, Method::Modular).unwrap(), lower);
    }
}

/// Targets on a facet, with its Hilbert bases, and off it, where nothing is below them.
#[test]
fn test_facet_target_points() {
    let mut input = input(H11_9, false, 10, H11_9_FACET);
    let expected_full = results(&input, Method::Modular).unwrap();
    let targets = two_targets(&input, &expected_full);
    let expected = below(&expected_full, &cone_data(&input), &targets);
    assert!(expected.len() < expected_full.len() / 2);
    input.problem.selection = Selection::TargetPoints {
        points: targets,
        max_deg: None,
    };
    assert_eq!(results(&input, Method::Modular).unwrap(), expected);
    assert_eq!(results(&input, Method::Reference).unwrap(), expected);
    // (0, 1, 0, ..., 0) has degree 1 and is off the facet's span, so nothing is below it,
    // which is taken for a mistake.
    let mut off = DMatrix::zeros(9, 1);
    off[(1, 0)] = 1;
    input.problem.selection = Selection::TargetPoints {
        points: off,
        max_deg: None,
    };
    assert!(matches!(
        results(&input, Method::Modular),
        Err(Error::Input {
            field: "target_points",
            ..
        })
    ));
}

/// On the facet the curve classes span a lattice of rank 8, whose coordinates key them, and
/// that gives the same invariants as the input coordinates: from the cones and from the
/// semigroup of the generators, below targets, one of them off the facet, certified, and with
/// keys of 512 bits.
#[test]
fn test_facet_lattice() {
    use cygv::modular::{self, CandidateSource, Tuning};
    let input = input(H11_9, false, 10, H11_9_FACET);
    let cone_data = cone_data(&input);
    let gens = generators(H11_9_FACET);
    let gens: Vec<&[i32]> = gens.as_slice().chunks(9).collect();
    let w = grading_vector(&input);
    let geometry = &input.problem.geometry;
    let geometry = modular::Geometry::from_toric(
        &geometry.q,
        &[],
        &geometry.intnums,
        cygv::CYKind::Threefold,
        &w,
    )
    .unwrap();
    let full = results(&input, Method::Modular).unwrap();
    let targets = two_targets(&input, &full);
    let mut off = [0; 9];
    off[1] = 1;
    let targets: Vec<&[i32]> = targets.as_slice().chunks(9).chain([&off[..]]).collect();
    let options = |lattice, certify, min_key_bits| modular::Options {
        certify,
        tuning: Tuning {
            lattice,
            min_key_bits,
            ..Tuning::default()
        },
        ..modular::Options::default()
    };
    for source in [
        CandidateSource::Cones(&cone_data),
        CandidateSource::Generators(&gens),
    ] {
        for targets in [&[][..], &targets[..2], &targets[..], &targets[2..]] {
            // Within a degree they come in the order of the keys, which depends on the basis.
            let run = |options| {
                let mut gv =
                    modular::compute_gv(source, &w, 10, targets, &geometry, &options).unwrap();
                gv.sort_by(|a, b| (a.degree, &a.curve).cmp(&(b.degree, &b.curve)));
                gv
            };
            let expected = run(options(false, false, 0));
            assert_eq!(expected.is_empty(), targets.len() == 1);
            for (certify, min_key_bits) in [(false, 0), (true, 0), (false, 512)] {
                assert_eq!(run(options(true, certify, min_key_bits)), expected);
            }
        }
    }
}

/// The Hilbert bases of the cones $K_T$ of P(1,1,2,2,2)[8], with the Mori cone spanned by
/// (1, 0) and (0, 1).
fn p11222_bases() -> Vec<CandidateCone> {
    let mut cones = Vec::new();
    for s in 0..6 {
        for t in s..6 {
            let divisors = if s == t { vec![s] } else { vec![s, t] };
            // Only the last divisor can be negative, on the curves (x, y) with x < 2y.
            let hilbert_basis = if divisors.contains(&5) {
                dmatrix![1, 0; 0, 1]
            } else {
                dmatrix![1, 2; 0, 1]
            };
            cones.push(CandidateCone {
                divisors,
                hilbert_basis,
            });
        }
    }
    cones
}

/// The GV invariants of P(1,1,2,2,2)[8], up to degree 6 or for `min_points`.
fn p11222(
    cone: MoriCone,
    hilbert_bases: Option<Vec<CandidateCone>>,
    min_points: Option<u32>,
    method: Method,
) -> Result<Results, Error> {
    let input = Input {
        name: None,
        problem: Problem {
            geometry: Geometry {
                q: DMatrix::from_row_slice(6, 2, &[0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, -2]),
                intnums: HashMap::from([((0, 0, 0), 8), ((0, 0, 1), 4)]),
                ..Geometry::default()
            },
            cone,
            selection: min_points.map_or(Selection::MaxDeg(6), Selection::MinPoints),
            grading_vector: Some(vec![1, 1]),
        },
        kind: InvariantKind::GV,
        options: Options {
            hilbert_bases,
            ..Options::default()
        },
    };
    results(&input, method)
}

#[test]
fn test_entry_points() {
    let gens = MoriCone::Generators(dmatrix![1, 0; 0, 1]);
    let ineqs = MoriCone::Inequalities {
        hyperplanes: dmatrix![1, 0; 0, 1],
        equations: DMatrix::zeros(2, 0),
    };
    let bases = p11222_bases();
    let expected = p11222(gens.clone(), None, None, Method::Reference).unwrap();
    assert_eq!(expected[&(vec![0, 1], 0)], "4");
    for method in [Method::Auto, Method::Modular, Method::Reference] {
        for cone in [&gens, &ineqs] {
            assert_eq!(p11222(cone.clone(), None, None, method).unwrap(), expected);
            assert_eq!(
                p11222(cone.clone(), Some(bases.clone()), None, method).unwrap(),
                expected
            );
        }
    }
    // The reference method takes the Mori cone from its inequalities too.
    assert_eq!(
        p11222(ineqs.clone(), None, Some(10), Method::Auto).unwrap(),
        p11222(gens.clone(), None, Some(10), Method::Auto).unwrap()
    );

    // The Hilbert bases are checked, against the inequalities or the generators.
    let field = |res: Result<Results, Error>| match res {
        Err(Error::Input { field, .. }) => field,
        res => panic!("{res:?}"),
    };
    let mut bad = bases.clone();
    bad.pop();
    for cone in [&gens, &ineqs] {
        assert_eq!(
            field(p11222(cone.clone(), Some(bad.clone()), None, Method::Auto)),
            "hilbert_bases"
        );
    }
    let mut bad = bases.clone();
    bad[0].hilbert_basis = dmatrix![0; 1];
    assert_eq!(
        field(p11222(ineqs.clone(), Some(bad), None, Method::Reference)),
        "hilbert_bases"
    );
    // The generators span the Mori cone, which (0, 1), in some of the Hilbert bases, is not in
    // with (2, -1) as a generator.
    assert_eq!(
        field(p11222(
            MoriCone::Generators(dmatrix![1, 2; 0, -1]),
            Some(bases),
            None,
            Method::Auto
        )),
        "hilbert_bases"
    );
    // Without a fan, there is no Mori cone to derive.
    assert_eq!(
        field(p11222(MoriCone::FromFan, None, None, Method::Auto)),
        "generators"
    );
}
