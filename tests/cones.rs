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

const H11_9: &str = include_str!("../benches/data/h11_9.yaml");
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

    // Without the generators only the modular method can run.
    let input = self::input(H11_9, false, 8, H11_9_CONES);
    assert_eq!(input.generators.shape(), (9, 0));
    assert_eq!(results(&input, Method::Auto).unwrap(), expected);
    assert_eq!(
        results(&input, Method::Reference),
        Err(ComputeError::MissingGenerators)
    );
}

/// Without an element of the Hilbert basis of the Mori cone among the generators, their
/// semigroup has a hole, and the invariants computed over it are not integers. The cones have
/// every lattice point.
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
    assert_eq!(
        results(&without_cones(&input), Method::Modular),
        Err(ComputeError::Modular(cygv::modular::Error::NotStable))
    );
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
        cone_data,
        InvariantKind::GV,
        CYKind::Threefold,
        max_deg,
        min_points,
        None,
        None,
        method,
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

    // Without generators, the requests that need the reference method fail.
    assert_eq!(
        p11222(none.clone(), None, None, Method::Auto),
        Err(ComputeError::MissingGenerators)
    );
    assert_eq!(
        p11222(none.clone(), Some(cd.clone()), None, Method::Reference),
        Err(ComputeError::MissingGenerators)
    );
    assert_eq!(
        p11222(none.clone(), Some(cd.clone()), Some(10), Method::Auto),
        Err(ComputeError::MissingGenerators)
    );
    // With them, such a request ignores the cone data.
    assert_eq!(
        p11222(gens.clone(), Some(cd.clone()), Some(10), Method::Auto).unwrap(),
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
