//! The invariants from the coefficients of the instanton series.
//!
//! The extraction ([`extraction`](super::extraction)) gives the coefficients $n_C$ of the
//! instanton series in the basis of the $q^C$. For a threefold $n_C = \deg(C)\, A_C$ with
//! $A_C = \sum_{m \mid C} \mathrm{GV}_{C/m} / m^3$, and for an n-fold
//! $n^{(k)}_C = A_{C,k} = \sum_{m \mid C} \mathrm{GV}_{C/m,k} / m^2$, where $m \mid C$ means
//! that $C/m$ is a lattice point. (For a threefold $A_C$ is the genus-zero GW invariant.)
//! Solving for the GV invariants,
//! $\mathrm{GV}_C = A_C - \sum_{m \ge 2,\, m \mid C} \mathrm{GV}_{C/m} / m^e$ with $e = 3$ or
//! $2$, where $m$ runs over the divisors of the gcd of the coordinates of $C$. The curves are
//! taken in increasing order of degree, so that the GV invariants of the $C/m$ are final
//! when $C$ is reached.
//!
//! A multiple of a curve with a nonzero GV invariant can have one of its own even when it has
//! no coefficient in the instanton series, so those multiples are added on the way.
//!
//! The GW invariants are the $A_C$ (for an n-fold, the $A_{C,k}$). They are rational, so they
//! are not computed modulo the primes but from the lifted GV invariants ([`gw_from_gv`]).

use super::field::Ring;
use super::key::KeyLayout;
use super::table::{SlotMap, SparsePoly};
use super::{CurveInvariants, Error};
use crate::CYKind;
use malachite::{Integer, Rational};
use rayon::prelude::*;

/// The greatest common divisor of the absolute values of the coordinates.
pub(crate) fn gcd(c: &[i32]) -> u32 {
    c.iter().fold(0u32, |mut a, &x| {
        let mut b = x.unsigned_abs();
        while b != 0 {
            (a, b) = (b, a % b);
        }
        a
    })
}

/// The GV invariants, from the coefficients $n_C$ of the instanton series at the curve
/// classes of degree at most `max_deg`, as a series with as many values per term as `coeffs`.
///
/// `inv` holds the inverses of the integers up to `max_deg` at least. The values are solved
/// for in the place of the coefficients, in order of degree, each looking up the $C/m$ by
/// binary search among the terms of its degree; only the multiples that have no coefficient
/// take a table of their own.
pub fn multicover<R: Ring, L: KeyLayout>(
    ring: &R,
    layout: &L,
    kind: CYKind,
    mut coeffs: SparsePoly<L::Key, R::Elem>,
    inv: &[R::Elem],
    max_deg: u32,
) -> Result<SparsePoly<L::Key, R::Elem>, Error> {
    if inv.len() <= max_deg as usize {
        return Err(Error::DimensionMismatch);
    }
    let width = coeffs.width();
    let exponent = if kind.is_threefold() { 3 } else { 2 };
    let inv_pow: Vec<R::Elem> = inv
        .iter()
        .map(|&x| (1..exponent).fold(x, |p, _| ring.mul(p, x)))
        .collect();

    // A at every curve of positive degree.
    coeffs.retain_terms(|_, d, _| d != 0 && d <= max_deg);
    if kind.is_threefold() {
        for i in 0..coeffs.len() {
            let d = coeffs.degs()[i] as usize;
            for v in coeffs.values_mut(i) {
                *v = ring.mul(*v, inv[d]);
            }
        }
    }
    // The terms of degree e are those from starts[e] to starts[e + 1].
    let starts: Vec<usize> = (0..=max_deg + 1)
        .map(|e| coeffs.degs().partition_point(|&d| d < e))
        .collect();
    let find = |coeffs: &SparsePoly<L::Key, R::Elem>, key: L::Key, deg: u32| {
        let (a, b) = (starts[deg as usize], starts[deg as usize + 1]);
        coeffs.keys()[a..b].binary_search(&key).ok().map(|i| a + i)
    };
    // The multiples without a coefficient, by degree.
    let mut extra: SlotMap<L::Key, R::Elem> = SlotMap::new(width, 16);
    let mut extra_by_deg: Vec<Vec<u32>> = vec![Vec::new(); max_deg as usize + 1];

    let (mut c, mut cn) = (vec![0; layout.dim()], vec![0; layout.dim()]);
    let mut gv = vec![R::Elem::default(); width];
    for d in 1..=max_deg {
        let main = (starts[d as usize]..starts[d as usize + 1]).map(Ok);
        let more = core::mem::take(&mut extra_by_deg[d as usize])
            .into_iter()
            .map(Err);
        for term in main.chain(more) {
            let key = match term {
                Ok(i) => {
                    gv.copy_from_slice(coeffs.values(i));
                    coeffs.keys()[i]
                }
                Err(s) => {
                    gv.copy_from_slice(extra.values(s));
                    extra.key(s)
                }
            };
            layout.unpack(key, &mut c);
            let g = gcd(&c);
            if g > 1 {
                for m in (2..=g).filter(|m| g.is_multiple_of(*m)) {
                    for (y, &x) in cn.iter_mut().zip(&c) {
                        *y = x / m as i32;
                    }
                    let low = layout.pack(&cn);
                    let vals = match find(&coeffs, low, d / m) {
                        Some(j) => coeffs.values(j),
                        None => match extra.get(low) {
                            Some(t) => extra.values(t),
                            None => continue,
                        },
                    };
                    for (v, &w) in gv.iter_mut().zip(vals) {
                        *v = ring.sub(*v, ring.mul(w, inv_pow[m as usize]));
                    }
                }
                match term {
                    Ok(i) => coeffs.values_mut(i).copy_from_slice(&gv),
                    Err(s) => extra.values_mut(s).copy_from_slice(&gv),
                }
            }
            if gv.iter().all(|&x| ring.is_zero(x)) {
                continue;
            }
            let mut multiple = key;
            for j in 2..=max_deg / d {
                multiple = multiple + key;
                if find(&coeffs, multiple, j * d).is_none() {
                    let (t, created) = extra.insert(multiple, j * d);
                    if created {
                        extra_by_deg[(j * d) as usize].push(t);
                    }
                }
            }
        }
    }
    coeffs.retain_terms(|_, _, v| v.iter().any(|&x| !ring.is_zero(x)));
    coeffs.merge_from(extra.take());
    Ok(coeffs)
}

/// The GW invariants $A_C = \sum_{m \mid C} \mathrm{GV}_{C/m} / m^e$ of the curve classes of
/// degree at most `max_deg`, from the GV invariants of those classes, with $e = 3$ for a
/// threefold and $2$ for an n-fold. The result is in increasing order of degree and then
/// curve class, and leaves out the classes where all of them vanish.
///
/// Every multiple is a term (its degree, the index of its GV invariant and the factor), the
/// terms are sorted in parallel by the class they give, which is computed in the comparisons
/// rather than held, and each run of equal classes is summed, the runs in parallel.
pub fn gw_from_gv(
    kind: CYKind,
    gv: &[CurveInvariants],
    max_deg: u32,
) -> Result<Vec<CurveInvariants<Rational>>, Error> {
    let exponent = if kind.is_threefold() { 3 } else { 2 };
    // The multiples must fit the curves' integers, as the classes up to the degree do when
    // the keys could hold them.
    let fits = gv.par_iter().all(|inv| {
        let m = (max_deg / inv.degree.max(1)) as i64;
        inv.curve
            .iter()
            .all(|&x| i32::try_from(x as i64 * m).is_ok())
    });
    if !fits {
        return Err(Error::Key(super::key::KeyError::BoundTooLarge));
    }
    let mut terms: Vec<(u32, u32, u32)> = gv
        .par_iter()
        .enumerate()
        .filter(|(_, inv)| inv.degree != 0 && inv.degree <= max_deg)
        .flat_map_iter(|(i, inv)| {
            (1..=max_deg / inv.degree).map(move |m| (m * inv.degree, i as u32, m))
        })
        .collect();
    let class =
        |&(_, i, m): &(u32, u32, u32)| gv[i as usize].curve.iter().map(move |&x| x * m as i32);
    terms.par_sort_unstable_by(|a, b| a.0.cmp(&b.0).then_with(|| class(a).cmp(class(b))));
    let same = |a: &(u32, u32, u32), b: &(u32, u32, u32)| a.0 == b.0 && class(a).eq(class(b));
    Ok(terms
        .par_chunk_by(same)
        .filter_map(|run| {
            let first = &gv[run[0].1 as usize];
            let mut values = vec![Rational::default(); first.values.len()];
            for &(_, i, m) in run {
                let denom = Integer::from(u64::from(m).pow(exponent));
                for (x, v) in values.iter_mut().zip(&gv[i as usize].values) {
                    *x += Rational::from_integers_ref(v, &denom);
                }
            }
            values.iter().any(|x| *x != 0).then(|| CurveInvariants {
                curve: class(&run[0]).collect(),
                degree: run[0].0,
                values,
            })
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modular::field::{primes, LaneField, Lanes};
    use crate::modular::key::{PackedKey, PackedLayout};
    use std::collections::HashMap;

    type E = Lanes<2>;

    fn ring() -> LaneField<2> {
        let mut ps = primes();
        LaneField::new([ps.next().unwrap(), ps.next().unwrap()])
    }

    #[test]
    fn test_gcd() {
        assert_eq!(gcd(&[0, 0]), 0);
        assert_eq!(gcd(&[0, -6, 4]), 2);
        assert_eq!(gcd(&[9, 6, -15]), 3);
        assert_eq!(gcd(&[7]), 7);
    }

    /// GV invariants on the quadrant of two coordinates, graded by (1, 2), turned into
    /// coefficients by the definition and back.
    #[test]
    fn test_round_trip() {
        let lf = ring();
        let layout = PackedLayout::<u128>::new(&[30, 30]).unwrap();
        let max_deg = 12;
        let inv = lf.inverses(max_deg as usize + 1).unwrap();
        let deg = |c: [i32; 2]| (c[0] + 2 * c[1]) as u32;
        let int = |v: [i64; 2]| v.map(|x| lf.int(x)).to_vec();
        for kind in [CYKind::Threefold, CYKind::Nfold] {
            let e = if kind.is_threefold() { 3 } else { 2 };
            let mut gvs: HashMap<[i32; 2], Vec<E>> = HashMap::from([
                ([1, 0], int([3, -1])),
                ([0, 1], int([5, 0])),
                ([1, 1], int([-7, 2])),
                ([2, 2], int([1, 1])),
                ([3, 0], int([0, 13])),
                ([1, 3], int([8, 0])),
                ([4, 2], int([2, 9])),
            ]);
            // Cancels the coefficient at (0, 2), which is then only reached as a multiple of
            // (0, 1).
            let cancel = lf.neg(lf.inv(lf.int(1 << e)).unwrap());
            gvs.insert(
                [0, 2],
                gvs[&[0, 1]].iter().map(|&x| lf.mul(x, cancel)).collect(),
            );

            // n_C = (deg C for a threefold) sum_{m | C} GV_{C/m} / m^e
            let mut map = SlotMap::new(2, 16);
            for (c, v) in &gvs {
                for m in 1..=max_deg as i32 {
                    let cm = [c[0] * m, c[1] * m];
                    if deg(cm) > max_deg {
                        break;
                    }
                    let f = if kind.is_threefold() {
                        lf.int(deg(cm) as i64)
                    } else {
                        lf.one()
                    };
                    let scale = lf.mul(f, lf.inv(lf.int((m as i64).pow(e))).unwrap());
                    let entry = map.get_or_insert(layout.pack(&cm), deg(cm));
                    for (x, &y) in entry.iter_mut().zip(v) {
                        lf.fma(x, y, scale);
                    }
                }
            }
            let coeffs = map.take();
            assert!(!coeffs.keys().contains(&layout.pack(&[0, 2])));

            let got = multicover(&lf, &layout, kind, coeffs, &inv, max_deg).unwrap();
            let got: HashMap<PackedKey, Vec<E>> = (0..got.len())
                .map(|i| (got.keys()[i], got.values(i).to_vec()))
                .collect();
            let expected: HashMap<PackedKey, Vec<E>> = gvs
                .iter()
                .map(|(c, v)| (layout.pack(c), v.clone()))
                .collect();
            assert_eq!(got, expected, "{kind:?}");
        }
    }

    #[test]
    fn test_short_table() {
        let lf = ring();
        let layout = PackedLayout::<u128>::new(&[30, 30]).unwrap();
        let inv = lf.inverses(5).unwrap();
        let coeffs = SparsePoly::new(1);
        let kind = CYKind::Threefold;
        assert_eq!(
            multicover(&lf, &layout, kind, coeffs.clone(), &inv, 5),
            Err(Error::DimensionMismatch)
        );
        assert!(multicover(&lf, &layout, kind, coeffs, &inv, 4)
            .unwrap()
            .is_empty());
    }

    #[test]
    fn test_gw_from_gv() {
        let inv = |curve: Vec<i32>, degree, values: Vec<i64>| CurveInvariants {
            curve,
            degree,
            values: values.into_iter().map(Integer::from).collect(),
        };
        let gv = [
            inv(vec![1, 0], 1, vec![8, 2]),
            inv(vec![0, 1], 2, vec![0, 3]),
            inv(vec![2, 0], 2, vec![-1, 0]),
            inv(vec![1, 1], 3, vec![5, 0]),
            // The degree that is too high is ignored.
            inv(vec![3, 1], 5, vec![1, 1]),
        ];
        let rat = |v: &[(i64, i64)]| -> Vec<Rational> {
            v.iter()
                .map(|&(n, d)| Rational::from_signeds(n, d))
                .collect()
        };
        let got = gw_from_gv(CYKind::Threefold, &gv, 4).unwrap();
        let expected = vec![
            (vec![1, 0], 1, rat(&[(8, 1), (2, 1)])),
            (vec![0, 1], 2, rat(&[(0, 1), (3, 1)])),
            // GV(2, 0) + GV(1, 0) / 8
            (vec![2, 0], 2, rat(&[(0, 1), (1, 4)])),
            (vec![1, 1], 3, rat(&[(5, 1), (0, 1)])),
            (vec![3, 0], 3, rat(&[(8, 27), (2, 27)])),
            (vec![0, 2], 4, rat(&[(0, 1), (3, 8)])),
            // GV(1, 0) / 64 + GV(2, 0) / 8
            (vec![4, 0], 4, rat(&[(0, 1), (1, 32)])),
        ];
        let got: Vec<_> = got
            .into_iter()
            .map(|inv| (inv.curve, inv.degree, inv.values))
            .collect();
        assert_eq!(got, expected);

        // Exponent 2 for an n-fold.
        let got = gw_from_gv(CYKind::Nfold, &gv[..1], 2).unwrap();
        assert_eq!(got[1].curve, [2, 0]);
        assert_eq!(got[1].values, rat(&[(2, 1), (1, 2)]));
    }
}
