//! Tests of vex phases of the ambient toric variety, which need strata.
//!
//! The geometries come from the vex references of cgv, the C program the modular pipeline
//! reimplements (`tests/data/refs_vex.yaml`): the charges and intersection numbers there, and
//! the Hilbert basis of the Mori cone and the strata that its `cgv_phase.py` derived from the
//! fans. The tests check properties of the results; `tests/fixtures.rs` compares them with
//! the stored tables.

use cygv::malachite::{Integer, Rational};
use cygv::modular;
use cygv::{
    compute_gv, Ambient, Error, Geometry, Method, MoriCone, Options, Problem, Selection, Stratum,
};
use nalgebra::{DMatrix, DVector};
use std::collections::HashMap;

/// A class of integers.
fn class(entries: &[i64]) -> Vec<Rational> {
    entries.iter().map(|&x| Rational::from(x)).collect()
}

/// The data of a threefold hypersurface, with the vectors as rows.
struct Phase {
    generators: &'static [&'static [i32]],
    grading_vector: &'static [i32],
    q: &'static [&'static [i32]],
    intnums: &'static [[i32; 4]],
    strata: Vec<Stratum>,
}

impl Phase {
    fn gv(&self, max_deg: u32, method: Method) -> Result<HashMap<Vec<i32>, Integer>, Error> {
        let columns = |rows: &[&[i32]]| {
            DMatrix::from_iterator(
                rows[0].len(),
                rows.len(),
                rows.iter().flat_map(|r| r.iter().copied()),
            )
        };
        let problem = Problem {
            geometry: Geometry {
                q: columns(self.q),
                intnums: self
                    .intnums
                    .iter()
                    .map(|&[a, b, c, v]| ((a as usize, b as usize, c as usize), v))
                    .collect(),
                nefpart: Vec::new(),
                ambient: if self.strata.is_empty() {
                    Ambient::Plain
                } else {
                    Ambient::Strata(self.strata.clone())
                },
            },
            cone: MoriCone::Generators(columns(self.generators)),
            selection: Selection::MaxDeg(max_deg),
            grading_vector: Some(self.grading_vector.to_vec()),
        };
        let options = Options {
            method,
            ..Options::default()
        };
        let res = compute_gv(&problem, &options)?;
        Ok(res
            .entries
            .into_iter()
            .map(|e| (e.curve, e.value))
            .collect())
    }
}

/// A vex phase with $h^{1,1} = 2$, a flop of the resolved $\mathbb{P}(1,1,2,2,2)[8]$, whose
/// fiber class carries the familiar invariants (`liam_h2_fan1`).
fn h2_vex() -> Phase {
    Phase {
        generators: &[&[0, 1], &[1, -1]],
        grading_vector: &[2, 1],
        q: &[&[1, 0, -1, 1, 0, 1], &[0, 1, 1, 1, 1, 0]],
        intnums: &[[0, 1, 1, 4], [1, 1, 1, 2]],
        strata: vec![Stratum {
            divisors: [1, 2, 4],
            class: class(&[1, -1]),
        }],
    }
}

#[test]
fn test_h2_vex() {
    let gv = h2_vex().gv(6, Method::Auto).unwrap();
    for (c, expected) in [
        ([0, 1], 640),
        ([0, 2], 10032),
        ([0, 3], 288384),
        ([1, -1], 1),
        ([1, 0], 20),
        ([1, 1], 6474),
    ] {
        assert_eq!(gv[&c[..]], expected, "{c:?}");
    }

    // Without its strata, the modular pipeline finds invariants that are not integers.
    let mut phase = h2_vex();
    phase.strata.clear();
    assert!(matches!(
        phase.gv(6, Method::Auto),
        Err(Error::Modular(modular::Error::NotInteger { .. }))
    ));

    // Both methods give the same invariants, and the strata are checked.
    assert_eq!(h2_vex().gv(6, Method::Modular).unwrap(), gv);
    assert_eq!(h2_vex().gv(6, Method::Reference).unwrap(), gv);
    let mut phase = h2_vex();
    phase.strata[0].divisors = [1, 2, 6];
    assert!(matches!(
        phase.gv(6, Method::Auto),
        Err(Error::Input {
            field: "strata",
            ..
        })
    ));
}

/// A vex phase with $h^{1,1} = 4$ and an FRST of a different polytope with the same CY
/// (`pair_h4_1`): the curve classes are related by `T`, and the invariants must agree.
#[test]
fn test_vex_and_frst_agree() {
    let vex = Phase {
        generators: &[
            &[0, 0, -1, 0],
            &[0, 0, 1, 1],
            &[0, 1, 0, -1],
            &[1, -1, 0, -1],
            &[1, 0, -2, -4],
        ],
        grading_vector: &[7, 3, -1, 2],
        q: &[
            &[1, 0, 0, 0, 2, 0, 0, 1],
            &[0, 1, 0, -2, 1, -1, 0, 1],
            &[0, 0, 1, 2, -1, 1, 0, -3],
            &[0, 0, 0, -1, 1, 0, 1, 0],
        ],
        intnums: &[
            [0, 1, 3, 1],
            [0, 2, 3, 1],
            [0, 0, 0, 8],
            [0, 0, 1, 1],
            [0, 0, 2, 3],
            [0, 0, 3, 2],
            [0, 1, 1, -1],
            [0, 2, 2, 1],
            [0, 3, 3, -4],
            [1, 1, 1, 1],
            [1, 1, 3, -1],
            [1, 3, 3, -1],
            [2, 3, 3, -2],
            [3, 3, 3, 7],
        ],
        strata: vec![Stratum {
            divisors: [3, 5, 6],
            class: class(&[0, 1, 0, -1]),
        }],
    };
    let frst = Phase {
        generators: &[&[0, 0, 0, 1], &[0, 0, 1, -3], &[0, 1, -1, 3], &[1, 0, 0, 0]],
        // The vex grading vector in the basis of the FRST, T^{-1} w.
        grading_vector: &[2, 2, 4, 1],
        q: &[
            &[1, 0, -1, -1, 1, 0, 4, 0],
            &[0, 1, 0, 1, 0, 0, -2, 0],
            &[0, 0, 3, -2, -2, 1, 0, 0],
            &[0, 0, 1, -1, -1, 0, 0, 1],
        ],
        intnums: &[
            [0, 1, 2, 6],
            [0, 0, 0, 8],
            [0, 0, 1, 16],
            [0, 0, 2, 3],
            [0, 1, 1, 28],
            [0, 2, 2, 1],
            [1, 1, 1, 48],
            [1, 1, 2, 10],
            [1, 2, 2, 2],
            [2, 2, 3, 1],
            [2, 3, 3, -3],
            [3, 3, 3, 9],
        ],
        strata: Vec::new(),
    };
    let t = DMatrix::from_row_slice(4, 4, &[1, 2, 0, 1, 0, 1, 0, 1, 0, -1, 1, -3, 0, 1, 0, 0]);

    let max_deg = 10;
    let gv_vex = vex.gv(max_deg, Method::Modular).unwrap();
    assert_eq!(vex.gv(max_deg, Method::Reference).unwrap(), gv_vex);
    let gv_frst = frst.gv(max_deg, Method::Auto).unwrap();
    let mapped: HashMap<Vec<i32>, Integer> = gv_vex
        .into_iter()
        .map(|(c, gv)| {
            (
                (t.transpose() * DVector::from_vec(c)).as_slice().to_vec(),
                gv,
            )
        })
        .collect();
    assert!(mapped.len() > 50);
    assert_eq!(mapped, gv_frst);
}
