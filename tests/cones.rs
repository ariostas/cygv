//! Tests of the cone data, with which the modular pipeline takes its candidates from the
//! Hilbert bases of the cones $K_T$.
//!
//! The cone data in `tests/data` was computed with Normaliz for `benches/data/h11_9.yaml` (see
//! the headers of the files), and is appended to its input.

#![cfg(feature = "cli")]

use cygv::dispatch::{compute_gvgw_strings, ComputeError};
use cygv::io::Input;
use cygv::{CYKind, CandidateCone, ConeData, InvariantKind, Method};
use nalgebra::{dmatrix, DMatrix, DVector, RowDVector};
use std::collections::HashMap;

const H11_8: &str = include_str!("../benches/data/h11_8.yaml");
const H11_9: &str = include_str!("../benches/data/h11_9.yaml");
const H11_10: &str = include_str!("../benches/data/h11_10.yaml");
const H11_11: &str = include_str!("../benches/data/h11_11.yaml");
const H11_9_CONES: &str = include_str!("data/h11_9_cones.yaml");
const H11_9_FACET: &str = include_str!("data/h11_9_facet.yaml");

/// The input of a bench file at `max_deg`, without its generators if `generators` is false,
/// with `extra` appended.
fn input(bench: &str, generators: bool, max_deg: u32, extra: &str) -> Input {
    let (head, gens) = bench.split_at(bench.find("\ngenerators:").unwrap() + 1);
    let head: String = head
        .lines()
        .filter(|l| !l.starts_with("max_deg:"))
        .map(|l| format!("{l}\n"))
        .collect();
    let gens = if generators { gens } else { "" };
    let mut inputs =
        Input::load_all(&format!("{head}max_deg: {max_deg}\n{gens}\n{extra}")).unwrap();
    assert_eq!(inputs.len(), 1);
    inputs.remove(0)
}

type Results = HashMap<(Vec<i32>, usize), String>;

fn results(input: &Input, method: Method) -> Result<Results, ComputeError> {
    let input = Input {
        method,
        ..input.clone()
    };
    Ok(input
        .compute()?
        .into_iter()
        .map(|((c, k), v)| ((c.as_slice().to_vec(), k), v))
        .collect())
}

/// The same input without the cone data.
fn without_cones(input: &Input) -> Input {
    Input {
        cone_data: None,
        ..input.clone()
    }
}

#[test]
fn test_h11_9() {
    let input = input(H11_9, true, 8, H11_9_CONES);
    let expected = results(&input, Method::Reference).unwrap();
    assert_eq!(expected.len(), 551);
    assert_eq!(results(&input, Method::Modular).unwrap(), expected);
    assert_eq!(results(&input, Method::Auto).unwrap(), expected);
    assert_eq!(
        results(&without_cones(&input), Method::Modular).unwrap(),
        expected
    );

    // Without the generators both methods take the Mori cone from the cone data.
    let input = self::input(H11_9, false, 8, H11_9_CONES);
    assert_eq!(input.generators.shape(), (9, 0));
    assert_eq!(results(&input, Method::Auto).unwrap(), expected);
    assert_eq!(results(&input, Method::Reference).unwrap(), expected);
}

/// Without an element of the Hilbert basis of the Mori cone among the generators, their
/// semigroup has a hole, over which the invariants would not be integers. Both methods use
/// every lattice point of the cone the generators span instead, with or without cone data.
#[test]
fn test_saturation() {
    let missing = "[1,0,0,0,0,0,1,0,0]";
    let bench: String = H11_9
        .lines()
        .filter(|l| l.trim().trim_start_matches("- ") != missing)
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(bench.len(), H11_9.len() - 2 * (missing.len() + 5));
    let input = input(&bench, true, 8, H11_9_CONES);
    let full = self::input(H11_9, true, 8, H11_9_CONES);
    let expected = results(&full, Method::Modular).unwrap();
    assert_eq!(expected[&(vec![1, 0, 0, 0, 0, 0, 1, 0, 0], 0)], "20");
    assert_eq!(results(&input, Method::Modular).unwrap(), expected);
    for method in [Method::Modular, Method::Reference] {
        assert_eq!(results(&without_cones(&input), method).unwrap(), expected);
    }
}

/// The cones of a facet give the invariants of the curves on it.
#[test]
fn test_facet() {
    let input = input(H11_9, false, 10, H11_9_FACET);
    let cone_data = input.cone_data.as_ref().unwrap();
    assert_eq!(cone_data.equations.shape(), (9, 1));
    let expected = results(&input, Method::Reference).unwrap();
    assert!(expected.len() > 100);
    let facet = [0, 2, 2, 2, 2, 2, 0, 0, 1];
    assert!(expected
        .keys()
        .all(|(c, _)| { c.iter().zip(facet).map(|(x, f)| x * f).sum::<i32>() == 0 }));
    assert_eq!(results(&input, Method::Modular).unwrap(), expected);
    assert_eq!(
        results(&without_cones(&input), Method::Modular).unwrap(),
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
        let input = input(bench, cones != H11_9_FACET, max_deg, cones);
        let normaliz = input.cone_data.unwrap();
        let w = input.grading_vector.as_slice();
        let found = cygv::cone::cone_data(
            &normaliz.hyperplanes,
            &normaliz.equations,
            &input.q,
            &[],
            w,
            max_deg,
        )
        .unwrap();
        assert_eq!(found.cones.len(), normaliz.cones.len());
        let degree = |c: &Vec<i32>| -> i32 { c.iter().zip(w).map(|(x, y)| x * y).sum() };
        for cone in &found.cones {
            let theirs = normaliz
                .cones
                .iter()
                .find(|c| c.divisors == cone.divisors)
                .unwrap();
            let mut expected = sorted_columns(&theirs.hilbert_basis);
            expected.retain(|c| degree(c) <= max_deg as i32);
            assert_eq!(
                sorted_columns(&cone.hilbert_basis),
                expected,
                "{:?}",
                cone.divisors
            );
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
        let input = input(bench, true, 8, cones);
        let normaliz = input.cone_data.unwrap();
        let found = cygv::cone::facets(&input.generators).unwrap();
        assert_eq!(found.equations.ncols(), 0);
        assert_eq!(
            sorted_columns(&found.hyperplanes),
            sorted_columns(&normaliz.hyperplanes)
        );
    }

    // On the facet, a hyperplane is only defined up to multiples of the equation, so
    // Normaliz's are moved along it to vanish where `facets` makes them vanish, at the
    // equation's first nonzero entry.
    let input = input(H11_9, false, 10, H11_9_FACET);
    let normaliz = input.cone_data.unwrap();
    let found = cygv::cone::facets(&input.generators).unwrap();
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
    let degree =
        |c: &[i32]| -> i32 { input.grading_vector.iter().zip(c).map(|(w, x)| w * x).sum() };
    let mut curves: Vec<&Vec<i32>> = results.keys().map(|(c, _)| c).collect();
    curves.sort_unstable_by_key(|c| (degree(c), c.to_vec()));
    let top = degree(curves[curves.len() - 1]);
    let lower = curves.iter().rev().find(|c| degree(c) < top).unwrap();
    let h11 = input.grading_vector.len();
    DMatrix::from_iterator(
        h11,
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
    let cone_data = full.cone_data.clone().unwrap();
    let targets = two_targets(&full, &expected_full);
    let expected = below(&expected_full, &cone_data, &targets);
    assert!(expected.len() < expected_full.len() / 2);
    for generators in [true, false] {
        let mut input = input(H11_9, generators, 8, H11_9_CONES);
        input.max_deg = None;
        input.target_points = Some(targets.clone());
        assert_eq!(results(&input, Method::Auto).unwrap(), expected);
        assert_eq!(results(&input, Method::Modular).unwrap(), expected);
        if generators {
            assert_eq!(
                results(&without_cones(&input), Method::Modular).unwrap(),
                expected
            );
        }
        // With a maximum degree below that of the first target.
        input.max_deg = Some(7);
        let lower: Results = expected
            .iter()
            .filter(|((c, _), _)| {
                c.iter()
                    .zip(full.grading_vector.iter())
                    .map(|(x, w)| x * w)
                    .sum::<i32>()
                    <= 7
            })
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        assert!(lower.len() < expected.len());
        assert_eq!(results(&input, Method::Modular).unwrap(), lower);
    }
}

/// Targets on a facet, with its cone data, and off it, where nothing is below them.
#[test]
fn test_facet_target_points() {
    let mut input = input(H11_9, false, 10, H11_9_FACET);
    let expected_full = results(&input, Method::Modular).unwrap();
    let cone_data = input.cone_data.clone().unwrap();
    let targets = two_targets(&input, &expected_full);
    let expected = below(&expected_full, &cone_data, &targets);
    assert!(expected.len() < expected_full.len() / 2);
    input.max_deg = None;
    input.target_points = Some(targets);
    assert_eq!(results(&input, Method::Modular).unwrap(), expected);
    assert_eq!(results(&input, Method::Reference).unwrap(), expected);
    // (0, 1, 0, ..., 0) has degree 1 and is off the facet's span.
    let mut off = DMatrix::zeros(9, 1);
    off[(1, 0)] = 1;
    input.target_points = Some(off);
    assert_eq!(results(&input, Method::Modular).unwrap(), Results::new());
}

/// The cone data of P(1,1,2,2,2)[8], with the Mori cone spanned by (1, 0) and (0, 1).
fn p11222_cone_data() -> ConeData {
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
    ConeData {
        hyperplanes: dmatrix![1, 0; 0, 1],
        equations: DMatrix::zeros(2, 0),
        cones,
    }
}

fn p11222(
    generators: DMatrix<i32>,
    cone_data: Option<ConeData>,
    min_points: Option<u32>,
    method: Method,
) -> Result<Vec<((DVector<i32>, usize), String)>, ComputeError> {
    let max_deg = if min_points.is_some() { None } else { Some(6) };
    let mut res = compute_gvgw_strings(
        generators,
        RowDVector::from_row_slice(&[1, 1]),
        DMatrix::from_row_slice(6, 2, &[0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, -2]),
        Vec::new(),
        HashMap::from([((0, 0, 0), 8), ((0, 0, 1), 4)]),
        Vec::new(),
        None,
        cone_data,
        InvariantKind::GV,
        CYKind::Threefold,
        max_deg,
        min_points,
        None,
        None,
        method,
        false,
    )?;
    res.sort_unstable_by(|a, b| a.0 .0.as_slice().cmp(b.0 .0.as_slice()));
    Ok(res)
}

#[test]
fn test_entry_points() {
    let gens = dmatrix![1, 0; 0, 1];
    let none = DMatrix::zeros(2, 0);
    let cd = p11222_cone_data();
    let expected = p11222(gens.clone(), None, None, Method::Reference).unwrap();
    assert!(expected.contains(&((DVector::from_row_slice(&[0, 1]), 0), "4".to_owned())));
    for method in [Method::Auto, Method::Modular] {
        assert_eq!(
            p11222(gens.clone(), Some(cd.clone()), None, method).unwrap(),
            expected
        );
        assert_eq!(
            p11222(none.clone(), Some(cd.clone()), None, method).unwrap(),
            expected
        );
    }

    // Without generators or cone data there is no Mori cone. With cone data, the requests
    // that need the reference method take it from there.
    assert_eq!(
        p11222(none.clone(), None, None, Method::Auto),
        Err(ComputeError::MissingGenerators)
    );
    assert_eq!(
        p11222(none.clone(), Some(cd.clone()), None, Method::Reference).unwrap(),
        expected
    );
    assert_eq!(
        p11222(none.clone(), Some(cd.clone()), Some(10), Method::Auto).unwrap(),
        p11222(gens.clone(), None, Some(10), Method::Auto).unwrap()
    );

    // The cone data is checked.
    let mut bad = cd.clone();
    bad.cones.pop();
    assert!(matches!(
        p11222(gens.clone(), Some(bad), None, Method::Auto),
        Err(ComputeError::InvalidConeData(_))
    ));
    let mut bad = cd.clone();
    bad.cones[0].hilbert_basis = dmatrix![0; 1];
    assert!(matches!(
        p11222(gens.clone(), Some(bad), None, Method::Reference),
        Err(ComputeError::InvalidConeData(_))
    ));
    assert!(matches!(
        p11222(dmatrix![1, -1; 0, 1], Some(cd), None, Method::Auto),
        Err(ComputeError::InvalidConeData(_))
    ));
}
