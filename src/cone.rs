//! Descriptions of the Mori cone and of the cones of curve classes that carry period data.
//!
//! The fundamental period vanishes at every curve class with more than two negative
//! intersections with the toric divisors. The others, the *candidates*, are the lattice points
//! of the cones
//! $$K_T = \text{Mori cone} \cap \{C : Q_r \cdot C \ge 0 \text{ for } r \notin T\}$$
//! over the sets $T$ of one or two divisors (and, in a vex phase, the divisors of each stratum,
//! whose curves of negative anticanonical degree carry the stratum's term), so the modular
//! pipeline can enumerate them from
//! the Hilbert bases of the $K_T$ without the rest of the Mori cone, which is much faster than
//! enumerating the semigroup of the generators. This module holds that data, [`ConeData`],
//! which callers compute with an external program such as normaliz for now, and the
//! conversion between the two descriptions of a cone, by generators and by inequalities
//! ([`facets`] and [`rays`]), which gives the hyperplanes of the Mori cone from its generators.

mod hull;

pub use hull::{facets, rays, Generators, Inequalities, TooLarge};

use crate::Stratum;
use nalgebra::DMatrix;
use std::collections::HashSet;

/// The Mori cone (or a face of it) by inequalities and equations, with the Hilbert bases of the
/// cones $K_T$.
///
/// In a vex phase, each stratum $S$ (see [`Stratum`]) also needs its cone $K_S$, of three
/// divisors, whose curves of negative anticanonical degree carry the stratum's term.
///
/// The data is trusted to be complete: a Hilbert basis that misses an element silently drops
/// curve classes. What can be checked cheaply is checked, i.e. that every $T$ with one or two
/// divisors and every stratum has a cone, that no other $T$ does, and that every Hilbert basis
/// element and generator lies in its cone.
///
/// The lattice points of the cones are taken in $\mathbb{Z}^{h^{1,1}}$, so the candidates are
/// those of the cone the hyperplanes describe, even when the generators do not generate all of
/// its lattice points. The results can then differ from those of the reference pipeline, which
/// only uses the semigroup the generators generate.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConeData {
    /// The inequalities $h \cdot C \ge 0$ of the Mori cone, one per column. Redundant ones are
    /// allowed.
    pub hyperplanes: DMatrix<i32>,
    /// The equations $e \cdot C = 0$ of the linear span of the cone, one per column, when it is
    /// a face of the Mori cone rather than all of it. Otherwise there are none.
    pub equations: DMatrix<i32>,
    /// The cone $K_T$ of every set $T$ of one or two divisors, and of the divisors of every
    /// stratum, in any order.
    pub cones: Vec<CandidateCone>,
}

/// The cone $K_T$ of the curve classes in the Mori cone whose negative intersections with the
/// toric divisors are all with divisors in $T$.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CandidateCone {
    /// The set $T$, as indices of rows of the charge matrix.
    pub divisors: Vec<usize>,
    /// The Hilbert basis of $K_T \cap \mathbb{Z}^{h^{1,1}}$, one element per column. Only the
    /// elements of degree at most the maximum degree are used, so the others can be left out.
    pub hilbert_basis: DMatrix<i32>,
}

/// The columns of a matrix, which is column-major.
pub(crate) fn columns(m: &DMatrix<i32>) -> impl Iterator<Item = &[i32]> {
    m.as_slice().chunks_exact(m.nrows().max(1)).take(m.ncols())
}

/// Whether `c` lies in the cone cut out by the columns of `hyperplanes` and `equations`.
fn in_cone(c: &[i32], hyperplanes: &DMatrix<i32>, equations: &DMatrix<i32>) -> bool {
    columns(hyperplanes).all(|h| dot(h, c) >= 0) && columns(equations).all(|e| dot(e, c) == 0)
}

fn dot(a: &[i32], b: &[i32]) -> i64 {
    a.iter().zip(b).map(|(&x, &y)| x as i64 * y as i64).sum()
}

/// Checks the cone data against the charge matrix `q` (one row per divisor), the generators,
/// which may be empty, and the strata, which must have been checked, returning the input field
/// at fault and what is wrong with it otherwise.
pub(crate) fn check_cone_data(
    cone_data: &ConeData,
    q: &DMatrix<i32>,
    generators: &DMatrix<i32>,
    strata: &[Stratum],
) -> Result<(), (&'static str, &'static str)> {
    let cones_err = |reason| Err(("cones", reason));
    let (n_divisors, h11) = q.shape();
    let ConeData {
        hyperplanes,
        equations,
        cones,
    } = cone_data;
    if hyperplanes.ncols() == 0 {
        return Err((
            "mori_hyperplanes",
            "the Mori cone needs at least one hyperplane",
        ));
    }
    if hyperplanes.nrows() != h11 {
        return Err((
            "mori_hyperplanes",
            "the hyperplanes must have one entry per Kähler parameter",
        ));
    }
    if equations.ncols() > 0 && equations.nrows() != h11 {
        return Err((
            "mori_equations",
            "the equations must have one entry per Kähler parameter",
        ));
    }
    if generators.ncols() > 0 && generators.nrows() == h11 {
        for g in columns(generators) {
            if !in_cone(g, hyperplanes, equations) {
                return Err((
                    "generators",
                    "a generator lies outside the Mori cone given by its hyperplanes and \
                     equations",
                ));
            }
        }
    }

    let strata_sets: Vec<[usize; 3]> = strata
        .iter()
        .map(|s| {
            let mut d = s.divisors;
            d.sort_unstable();
            d
        })
        .collect();
    let mut seen = HashSet::new();
    let q_rows: Vec<Vec<i32>> = q.row_iter().map(|r| r.iter().copied().collect()).collect();
    for cone in cones {
        let mut t = cone.divisors.clone();
        t.sort_unstable();
        if !(1..=3).contains(&t.len()) {
            return cones_err("a cone K_T must have one or two divisors, or three for a stratum");
        }
        if t[t.len() - 1] >= n_divisors || t.windows(2).any(|w| w[0] == w[1]) {
            return cones_err("the divisors of a cone K_T must be distinct rows of the charges");
        }
        if t.len() == 3 && !strata_sets.iter().any(|s| s[..] == t[..]) {
            return cones_err("a cone K_T of three divisors must be on a stratum");
        }
        if !seen.insert(t.clone()) {
            return cones_err("two cones K_T have the same divisors");
        }
        let hb = &cone.hilbert_basis;
        if hb.ncols() > 0 && hb.nrows() != h11 {
            return cones_err(
                "the Hilbert basis elements must have one entry per Kähler parameter",
            );
        }
        for g in columns(hb) {
            let outside = !in_cone(g, hyperplanes, equations)
                || (0..n_divisors).any(|r| !t.contains(&r) && dot(&q_rows[r], g) < 0);
            if outside {
                return cones_err("a Hilbert basis element lies outside its cone K_T");
            }
        }
    }
    if seen.len() != n_divisors + n_divisors * n_divisors.saturating_sub(1) / 2 + strata.len() {
        return cones_err(
            "there must be a cone K_T for every set of one or two divisors, and for every stratum",
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::dmatrix;

    /// The charges of P(1,1,2,2,2)[8], one row per divisor, and its cone data, with the Mori
    /// cone spanned by (1, 0) and (0, 1).
    fn p11222() -> (DMatrix<i32>, ConeData) {
        let q = DMatrix::from_row_slice(6, 2, &[0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, -2]);
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
        let cone_data = ConeData {
            hyperplanes: dmatrix![1, 0; 0, 1],
            equations: DMatrix::zeros(2, 0),
            cones,
        };
        (q, cone_data)
    }

    #[test]
    fn test_check_cone_data() {
        let (q, cd) = p11222();
        let gens = dmatrix![1, 0; 0, 1];
        assert_eq!(check_cone_data(&cd, &q, &gens, &[]), Ok(()));
        assert_eq!(check_cone_data(&cd, &q, &DMatrix::zeros(0, 0), &[]), Ok(()));

        let check = |f: &dyn Fn(&mut ConeData)| {
            let mut cd = cd.clone();
            f(&mut cd);
            check_cone_data(&cd, &q, &gens, &[])
        };
        assert!(check(&|cd| cd.hyperplanes = DMatrix::zeros(2, 0)).is_err());
        assert!(check(&|cd| cd.hyperplanes = dmatrix![1; 0; 0]).is_err());
        assert!(check(&|cd| cd.equations = dmatrix![1; 0; 0]).is_err());
        // The generator (1, 0) is not on the face x = 0.
        assert!(check(&|cd| cd.equations = dmatrix![1; 0]).is_err());
        assert!(check(&|cd| cd.cones.pop().map(|_| ()).unwrap()).is_err());
        assert!(check(&|cd| cd.cones[0].divisors = vec![]).is_err());
        assert!(check(&|cd| cd.cones[0].divisors = vec![0, 1, 2]).is_err());
        assert!(check(&|cd| cd.cones[0].divisors = vec![6]).is_err());
        // The order of the divisors does not matter, but repeating a set does.
        assert_eq!(check(&|cd| cd.cones[1].divisors = vec![1, 0]), Ok(()));
        assert!(check(&|cd| cd.cones[1].divisors = vec![0]).is_err());
        assert!(check(&|cd| cd.cones[1].divisors = vec![1, 1]).is_err());
        assert!(check(&|cd| cd.cones[0].hilbert_basis = dmatrix![1; 0; 0]).is_err());
        // (0, 1) meets divisor 5 negatively, which cone {0} does not allow.
        assert!(check(&|cd| cd.cones[0].hilbert_basis = dmatrix![0; 1]).is_err());
        assert!(check(&|cd| cd.cones[0].hilbert_basis = dmatrix![-1; 1]).is_err());
        assert_eq!(
            check(&|cd| cd.cones[0].hilbert_basis = DMatrix::zeros(0, 0)),
            Ok(())
        );
        assert!(check_cone_data(&cd, &q, &dmatrix![1, -1; 0, 1], &[]).is_err());

        // A stratum needs its cone, of three divisors, which other sets cannot have.
        let strata = [Stratum {
            divisors: [5, 0, 1],
            class: vec![(1, 1), (0, 1)],
        }];
        assert!(check_cone_data(&cd, &q, &gens, &strata).is_err());
        let mut with_stratum = cd.clone();
        with_stratum.cones.push(CandidateCone {
            divisors: vec![0, 1, 5],
            hilbert_basis: dmatrix![1, 0; 0, 1],
        });
        assert_eq!(check_cone_data(&with_stratum, &q, &gens, &strata), Ok(()));
        assert!(check_cone_data(&with_stratum, &q, &gens, &[]).is_err());
        with_stratum.cones.last_mut().unwrap().divisors = vec![0, 2, 5];
        assert!(check_cone_data(&with_stratum, &q, &gens, &strata).is_err());
    }
}
