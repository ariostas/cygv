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
//! taken in increasing order of that gcd, so that the GV invariants of the $C/m$ are final
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
use std::collections::HashMap;

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
/// `inv` holds the inverses of the integers up to `max_deg` at least.
pub fn multicover<R: Ring, L: KeyLayout>(
    ring: &R,
    layout: &L,
    kind: CYKind,
    coeffs: &SparsePoly<L::Key, R::Elem>,
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

    // A at every curve of positive degree, and the curves by the gcd of their coordinates.
    let mut map = SlotMap::new(width, 2 * coeffs.len());
    let mut by_gcd: Vec<Vec<u32>> = vec![Vec::new(); max_deg as usize + 1];
    let mut c = vec![0; layout.dim()];
    for i in 0..coeffs.len() {
        let d = coeffs.degs()[i];
        if d == 0 || d > max_deg {
            continue;
        }
        let (s, _) = map.insert(coeffs.keys()[i], d);
        let vals = map.values_mut(s);
        vals.copy_from_slice(coeffs.values(i));
        if kind.is_threefold() {
            for v in vals {
                *v = ring.mul(*v, inv[d as usize]);
            }
        }
        layout.unpack(coeffs.keys()[i], &mut c);
        by_gcd[gcd(&c) as usize].push(s);
    }

    let mut cn = vec![0; layout.dim()];
    let mut gv = vec![R::Elem::default(); width];
    for g in 1..=max_deg as usize {
        for s in core::mem::take(&mut by_gcd[g]) {
            gv.copy_from_slice(map.values(s));
            let key = map.key(s);
            if g > 1 {
                layout.unpack(key, &mut c);
                for m in (2..=g).filter(|m| g.is_multiple_of(*m)) {
                    for (y, &x) in cn.iter_mut().zip(&c) {
                        *y = x / m as i32;
                    }
                    if let Some(t) = map.get(layout.pack(&cn)) {
                        for (v, &w) in gv.iter_mut().zip(map.values(t)) {
                            *v = ring.sub(*v, ring.mul(w, inv_pow[m]));
                        }
                    }
                }
                map.values_mut(s).copy_from_slice(&gv);
            }
            if gv.iter().all(|&x| ring.is_zero(x)) {
                continue;
            }
            let d = map.deg(s);
            let mut multiple = key;
            for j in 2..=max_deg / d {
                multiple = multiple + key;
                let (t, created) = map.insert(multiple, j * d);
                if created {
                    by_gcd[j as usize * g].push(t);
                }
            }
        }
    }
    Ok(map.take())
}

/// The GW invariants $A_C = \sum_{m \mid C} \mathrm{GV}_{C/m} / m^e$ of the curve classes of
/// degree at most `max_deg`, from the GV invariants of those classes, with $e = 3$ for a
/// threefold and $2$ for an n-fold. The result is in increasing order of degree and then
/// curve class, and leaves out the classes where all of them vanish.
pub fn gw_from_gv(
    kind: CYKind,
    gv: &[CurveInvariants],
    max_deg: u32,
) -> Vec<CurveInvariants<Rational>> {
    let exponent = if kind.is_threefold() { 3 } else { 2 };
    let mut gw: HashMap<Vec<i32>, (u32, Vec<Rational>)> = HashMap::new();
    for inv in gv.iter().filter(|inv| inv.degree != 0) {
        for m in 1..=max_deg / inv.degree {
            let denom = Integer::from(u64::from(m).pow(exponent));
            let curve = inv.curve.iter().map(|&x| x * m as i32).collect();
            let (_, values) = gw
                .entry(curve)
                .or_insert_with(|| (m * inv.degree, vec![Rational::default(); inv.values.len()]));
            for (x, v) in values.iter_mut().zip(&inv.values) {
                *x += Rational::from_integers_ref(v, &denom);
            }
        }
    }
    let mut gw: Vec<CurveInvariants<Rational>> = gw
        .into_iter()
        .filter(|(_, (_, values))| values.iter().any(|x| *x != 0))
        .map(|(curve, (degree, values))| CurveInvariants {
            curve,
            degree,
            values,
        })
        .collect();
    gw.sort_unstable_by(|a, b| (a.degree, &a.curve).cmp(&(b.degree, &b.curve)));
    gw
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modular::field::{primes, LaneField, Lanes};
    use crate::modular::key::{PackedKey, PackedLayout};

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
        let layout = PackedLayout::new(&[30, 30]).unwrap();
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

            let got = multicover(&lf, &layout, kind, &coeffs, &inv, max_deg).unwrap();
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
        let layout = PackedLayout::new(&[30, 30]).unwrap();
        let inv = lf.inverses(5).unwrap();
        let coeffs = SparsePoly::new(1);
        let kind = CYKind::Threefold;
        assert_eq!(
            multicover(&lf, &layout, kind, &coeffs, &inv, 5),
            Err(Error::DimensionMismatch)
        );
        assert!(multicover(&lf, &layout, kind, &coeffs, &inv, 4)
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
        let got = gw_from_gv(CYKind::Threefold, &gv, 4);
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
        let got = gw_from_gv(CYKind::Nfold, &gv[..1], 2);
        assert_eq!(got[1].curve, [2, 0]);
        assert_eq!(got[1].values, rat(&[(2, 1), (1, 2)]));
    }
}
