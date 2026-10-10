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

use super::field::Ring;
use super::key::KeyLayout;
use super::table::{SlotMap, SparsePoly};
use super::Error;
use crate::CYKind;

/// The greatest common divisor of the absolute values of the coordinates.
fn gcd(c: &[i32]) -> u32 {
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
}
