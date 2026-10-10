//! Comparisons with the reference pipeline.
//!
//! Stage by stage: the candidates, period and instanton stages are run on the same inputs as
//! the reference's `Semigroup::with_max_degree`, `compute_omega` and
//! `compute_instanton_data`, and $c_0$, $c_{1,a}$, $S_2^{(k)}$ and $I^{(k)}$ are compared at
//! every curve class, with the reference's rationals reduced modulo each prime. A mismatch is
//! thereby pinned to a stage. These comparisons hold for any input, consistent geometry or
//! not.
//!
//! End to end: [`compute_gv`] and `run_hkty` must give the same GV invariants. This only
//! holds for consistent geometries: on others the reference reads each invariant off the
//! first nonzero coordinate of the curve, and the modular pipeline off the grading-weighted
//! sum, which disagree.

use super::candidates::candidates_by_enumeration;
use super::field::{primes, LaneField, Lanes};
use super::instanton::instanton;
use super::key::{coordinate_bounds, KeyLayout, PackedKey, PackedLayout};
use super::period::period;
use super::table::SparsePoly;
use super::{compute_gv, Geometry, Options};
use crate::reference::fundamental_period::compute_omega;
use crate::reference::instanton::compute_instanton_data;
use crate::reference::misc::process_int_nums;
use crate::reference::polynomial::properties::PolynomialProperties;
use crate::reference::run_hkty;
use crate::reference::semigroup::Semigroup;
use crate::{CYKind, InvariantKind};
use nalgebra::{DMatrix, DVector, RowDVector};
use rug::{Integer, Rational};
use std::collections::HashMap;

const NL: usize = 2;

struct Case {
    name: &'static str,
    generators: DMatrix<i32>,
    grading_vector: Vec<i32>,
    q: DMatrix<i32>,
    nefpart: Vec<DVector<i32>>,
    intnums: HashMap<(usize, usize, usize), i32>,
    cy_kind: CYKind,
    max_deg: u32,
}

/// The residue of a rational modulo `p`.
fn residue(x: &Rational, p: u64) -> u64 {
    let p = Integer::from(p);
    let (num, den) = x.clone().into_numer_denom();
    (num * den.invert(&p).unwrap()).modulo(&p).to_u64().unwrap()
}

/// Checks that `ours` holds, at every element of the semigroup, the residues of the values
/// that `expected` gives for each component, and nothing anywhere else.
fn compare(
    case: &str,
    what: &str,
    lf: &LaneField<NL>,
    keys: &[PackedKey],
    ours: &SparsePoly<PackedKey, Lanes<NL>>,
    expected: impl Fn(usize, usize) -> Rational,
) {
    let index: HashMap<_, _> = ours
        .keys()
        .iter()
        .enumerate()
        .map(|(i, &k)| (k, i))
        .collect();
    let mut found = 0;
    let mut nonzero = 0;
    for (e, key) in keys.iter().enumerate() {
        let i = index.get(key).copied();
        found += i.is_some() as usize;
        for comp in 0..ours.width() {
            let x = expected(e, comp);
            nonzero += (x != 0) as usize;
            let got = i.map_or([0; NL], |i| lf.residues(ours.values(i)[comp]));
            for (l, f) in lf.fields().iter().enumerate() {
                assert_eq!(
                    got[l],
                    residue(&x, f.p()),
                    "{case}: {what}, element {e}, component {comp}, lane {l}"
                );
            }
        }
    }
    assert_eq!(
        found,
        ours.len(),
        "{case}: {what} has terms outside the semigroup"
    );
    assert!(nonzero > 0, "{case}: {what} vanishes, which tests nothing");
}

fn check(case: &Case) {
    let mut ps = primes();
    let lf = LaneField::new([ps.next().unwrap(), ps.next().unwrap()]);
    let h11 = case.grading_vector.len();
    let name = case.name;

    // The reference.
    let sg = Semigroup::with_max_degree(
        case.generators.clone(),
        RowDVector::from_row_slice(&case.grading_vector),
        case.max_deg,
    )
    .unwrap();
    let pp = PolynomialProperties::<Rational>::new(&sg);
    let (intnum_dict, idx_pairs, n_indices) =
        process_int_nums(case.intnums.clone(), case.cy_kind).unwrap();
    let fp = compute_omega(&pp, &sg, &case.q, &case.nefpart, &idx_pairs).unwrap();
    let ref_c0 = fp.c0.clone(&pp.zero);
    let ref_c1: Vec<_> = fp.c1.iter().map(|p| p.clone(&pp.zero)).collect();
    let ref_c2: HashMap<_, _> = fp.c2.iter().map(|(&k, p)| (k, p.clone(&pp.zero))).collect();
    let ref_inst =
        compute_instanton_data(fp, &pp, &idx_pairs, n_indices, &intnum_dict, case.cy_kind)
            .unwrap()
            .inst;

    // The modular stages.
    let geometry = Geometry::from_toric(
        &case.q,
        &case.nefpart,
        &case.intnums,
        case.cy_kind,
        &case.grading_vector,
    )
    .unwrap();
    let gens: Vec<&[i32]> = case.generators.as_slice().chunks(h11).collect();
    let bounds =
        coordinate_bounds(gens.iter().copied(), &case.grading_vector, case.max_deg).unwrap();
    let layout = PackedLayout::new(&bounds).unwrap();
    let cands = candidates_by_enumeration(
        &layout,
        gens.iter().copied(),
        &case.grading_vector,
        case.max_deg,
        &geometry,
    )
    .unwrap();
    let pd = period(&lf, &layout, &geometry, &cands, case.max_deg).unwrap();
    let inst = instanton(&lf, &geometry, &pd, case.max_deg).unwrap();

    // The candidates are the elements with at most two negative intersections.
    let keys: Vec<PackedKey> = sg
        .elements
        .column_iter()
        .map(|c| layout.pack(c.as_slice()))
        .collect();
    let n_cands = sg
        .elements
        .column_iter()
        .filter(|c| geometry.n_negative(c.as_slice()) <= 2)
        .count();
    assert_eq!(cands.len(), n_cands, "{name}: candidates");
    assert!(
        cands.iter().all(|(k, _)| keys.contains(k)),
        "{name}: candidates"
    );

    let zero = Rational::new();
    let at = |p: &crate::Polynomial<Rational>, e: usize| p.coeffs.get(&e).unwrap_or(&zero).clone();
    compare(name, "c0", &lf, &keys, &pd.c0, |e, _| at(&ref_c0, e));
    compare(name, "c1", &lf, &keys, &pd.c1, |e, a| at(&ref_c1[a], e));
    compare(name, "S2", &lf, &keys, &pd.s2, |e, k| {
        let form = geometry.form(k);
        let mut s = Rational::new();
        for a in 0..h11 {
            for b in a..h11 {
                let kab = form[a * h11 + b];
                if kab == 0 {
                    continue;
                }
                let c2 = at(&ref_c2[&(a, b)], e);
                s += c2 * kab / if a == b { 2 } else { 1 };
            }
        }
        s
    });
    compare(name, "I", &lf, &keys, &inst.inst, |e, k| {
        match case.cy_kind {
            CYKind::Threefold => {
                let mut s = Rational::new();
                for (t, p) in ref_inst.iter().enumerate() {
                    s += at(p, e) * case.grading_vector[t];
                }
                s
            }
            CYKind::Nfold => at(&ref_inst[k], e),
        }
    });
}

/// Checks that the modular pipeline gives the GV invariants of the reference pipeline.
fn check_gv(case: &Case, options: &Options) {
    let reference = run_hkty::<Rational>(
        case.generators.clone(),
        RowDVector::from_row_slice(&case.grading_vector),
        InvariantKind::GV,
        case.cy_kind,
        Some(case.max_deg),
        None,
        None,
        case.q.clone(),
        case.nefpart.clone(),
        case.intnums.clone(),
        None,
    );
    let expected: HashMap<(Vec<i32>, usize), String> = reference
        .into_iter()
        .map(|((c, k), gv)| ((c.as_slice().to_vec(), k), gv.to_string()))
        .collect();

    let geometry = Geometry::from_toric(
        &case.q,
        &case.nefpart,
        &case.intnums,
        case.cy_kind,
        &case.grading_vector,
    )
    .unwrap();
    let h11 = case.grading_vector.len();
    let gens: Vec<&[i32]> = case.generators.as_slice().chunks(h11).collect();
    let result = compute_gv(
        &gens,
        &case.grading_vector,
        case.max_deg,
        &geometry,
        options,
    )
    .unwrap();
    let mut got = HashMap::new();
    for inv in result {
        for (k, v) in inv.values.iter().enumerate() {
            if *v != 0.into() {
                got.insert((inv.curve.clone(), k), v.to_string());
            }
        }
    }
    assert!(!expected.is_empty(), "{}: no invariants", case.name);
    assert_eq!(got, expected, "{}", case.name);
}

fn quintic(max_deg: u32) -> Case {
    Case {
        name: "quintic",
        generators: DMatrix::from_column_slice(1, 1, &[1]),
        grading_vector: vec![1],
        q: DMatrix::from_column_slice(5, 1, &[1; 5]),
        nefpart: vec![],
        intnums: HashMap::from([((0, 0, 0), 5)]),
        cy_kind: CYKind::Threefold,
        max_deg,
    }
}

/// $\mathbb{P}(1,1,2,2,2)[8]$, whose second generator meets the exceptional divisor
/// negatively.
fn p11222(max_deg: u32) -> Case {
    Case {
        name: "P11222",
        generators: DMatrix::from_column_slice(2, 2, &[1, 0, 0, 1]),
        grading_vector: vec![1, 1],
        q: DMatrix::from_row_slice(6, 2, &[0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, -2]),
        nefpart: vec![],
        intnums: HashMap::from([((0, 0, 0), 8), ((0, 0, 1), 4)]),
        cy_kind: CYKind::Threefold,
        max_deg,
    }
}

#[test]
fn test_quintic() {
    check(&quintic(10));
}

/// The GV invariants of degree 7 and up no longer fit one prime, so the default options take
/// a second pass.
#[test]
fn test_quintic_gv() {
    check_gv(&quintic(12), &Options::default());
}

/// One prime per pass and at least five primes, so that the lift merges five passes and
/// goes through big integers.
#[test]
fn test_several_passes() {
    let options = Options {
        lanes: 1,
        min_primes: 5,
        ..Options::default()
    };
    check_gv(&quintic(12), &options);
    check_gv(&p11222(8), &options);
    // Two check primes, with three lanes.
    let options = Options {
        lanes: 3,
        check_primes: 2,
        ..Options::default()
    };
    check_gv(&quintic(12), &options);
}

#[test]
fn test_p11222_gv() {
    check_gv(&p11222(12), &Options::default());
}

#[test]
fn test_invalid_options() {
    let case = p11222(4);
    let geometry =
        Geometry::from_toric(&case.q, &[], &case.intnums, case.cy_kind, &[1, 1]).unwrap();
    let gens: Vec<&[i32]> = case.generators.as_slice().chunks(2).collect();
    for options in [
        Options {
            check_primes: 0,
            ..Options::default()
        },
        Options {
            lanes: 5,
            ..Options::default()
        },
    ] {
        assert_eq!(
            compute_gv(&gens, &[1, 1], 4, &geometry, &options),
            Err(super::Error::InvalidOptions)
        );
    }
    // Never stable with a single pass of one prime.
    let options = Options {
        lanes: 1,
        max_passes: 1,
        ..Options::default()
    };
    assert_eq!(
        compute_gv(&gens, &[1, 1], 4, &geometry, &options),
        Err(super::Error::NotStable)
    );
}

#[test]
fn test_p11222() {
    check(&p11222(10));
}

/// The toy threefold of `examples/threefold.yaml`, which is not a consistent geometry but
/// has curves with one and two negative intersections.
#[test]
fn test_toy_threefold() {
    check(&Case {
        name: "toy threefold",
        generators: DMatrix::from_column_slice(2, 2, &[0, -1, 1, 2]),
        grading_vector: vec![3, -1],
        q: DMatrix::from_column_slice(6, 2, &[1, 1, 1, 0, 1, 2, 0, 0, -1, 1, 1, -1]),
        nefpart: vec![],
        intnums: HashMap::from([
            ((0, 0, 0), 2),
            ((0, 0, 1), 1),
            ((0, 1, 1), -1),
            ((1, 1, 1), 5),
        ]),
        cy_kind: CYKind::Threefold,
        max_deg: 12,
    });
}

/// Loads one of the inputs that the command line interface reads.
#[cfg(feature = "cli")]
fn load(name: &'static str, yaml: &str, max_deg: u32) -> Case {
    let input = crate::io::Input::load_all(yaml).unwrap().remove(0);
    Case {
        name,
        generators: input.generators,
        grading_vector: input.grading_vector.iter().copied().collect(),
        q: input.q,
        nefpart: input.nefpart,
        intnums: input.intnums,
        cy_kind: input.cy_kind,
        max_deg,
    }
}

/// A realistic threefold with thousands of redundant generators.
#[test]
#[cfg(feature = "cli")]
fn test_h11_9() {
    check(&load(
        "h11_9",
        include_str!("../../benches/data/h11_9.yaml"),
        5,
    ));
}

/// The fourfold of the examples, a complete intersection with two anticanonical rows and
/// fourteen reference surfaces. The stages are also compared on the same data without the
/// nef partition, read as a fivefold hypersurface: that is not a consistent geometry, but it
/// covers an n-fold with a single row.
#[test]
#[cfg(feature = "cli")]
fn test_fourfold() {
    let mut case = load("fourfold", include_str!("../../examples/fourfold.yaml"), 5);
    assert_eq!(case.nefpart.len(), 2);
    check(&case);
    case.name = "fivefold hypersurface";
    case.nefpart.clear();
    check(&case);
}

/// The realistic threefolds of the benchmarks, against the reference end to end, one test
/// each so that they run in parallel.
#[cfg(feature = "cli")]
macro_rules! bench_gv_tests {
    ($($test:ident: $name:literal, $max_deg:literal;)*) => {$(
        #[test]
        fn $test() {
            let yaml = include_str!(concat!("../../benches/data/", $name, ".yaml"));
            check_gv(&load($name, yaml, $max_deg), &Options::default());
        }
    )*};
}

#[cfg(feature = "cli")]
bench_gv_tests! {
    test_h11_8_gv: "h11_8", 6;
    test_h11_9_gv: "h11_9", 6;
    test_h11_10_gv: "h11_10", 5;
    test_h11_11_gv: "h11_11", 5;
}

/// The fourfold of the examples, a complete intersection. (Without its nef partition the
/// reference finds GV invariants that are not integers.)
#[test]
#[cfg(feature = "cli")]
fn test_fourfold_gv() {
    let case = load("fourfold", include_str!("../../examples/fourfold.yaml"), 10);
    check_gv(&case, &Options::default());
}
