//! The series that the invariants are extracted from.
//!
//! From the period data on the candidates this computes
//!
//! - $\alpha_a = c_{1,a} / c_0$, as one series with $h^{1,1}$ values per term, so that its
//!   support is the union of the supports of the $\alpha_a$;
//! - the instanton series $I^{(k)} = S_2^{(k)} / c_0 - \sum_a \alpha_a V^{(k)}_a$ with
//!   $V^{(k)}_a = \frac12 \sum_b K^{(k)}_{ab}\, \alpha_b$, one value per bilinear form.
//!
//! For a threefold this is $I = \sum_t w_t\, \mathrm{inst}_t$ in terms of the reference
//! pipeline's instanton corrections, and for an n-fold $I^{(k)} = \mathrm{inst}_k$; in both
//! cases it replaces the $\beta$ and $F$ series and the per-index products of the reference.
//! All series are truncated at the maximum degree.

use super::field::Ring;
use super::key::Key;
use super::period::PeriodData;
use super::table::{DegreeHeap, SlotMap, SparsePoly};
use super::{Error, Geometry};

/// The output of the instanton stage.
#[derive(Debug, Clone)]
pub struct InstantonData<K, E> {
    /// $\alpha_a$, $h^{1,1}$ values per term.
    pub alpha: SparsePoly<K, E>,
    /// $I^{(k)}$, one value per bilinear form.
    pub inst: SparsePoly<K, E>,
}

/// Adds $f(a_x, b_y)$ to the slot of $x + y$ in `acc` for every pair of terms with total
/// degree at most `max_deg`.
///
/// `f` gets the values of the two terms and those of the slot.
pub fn mul_into<K: Key, E: Copy + Default + PartialEq>(
    a: &SparsePoly<K, E>,
    b: &SparsePoly<K, E>,
    max_deg: u32,
    acc: &mut SlotMap<K, E>,
    mut f: impl FnMut(&[E], &[E], &mut [E]),
) {
    for x in 0..a.len() {
        let dx = a.degs()[x];
        if dx > max_deg {
            break;
        }
        let (kx, vx) = (a.keys()[x], a.values(x));
        for y in 0..b.len_upto(max_deg - dx) {
            let out = acc.get_or_insert(kx + b.keys()[y], dx + b.degs()[y]);
            f(vx, b.values(y), out);
        }
    }
}

/// The product of a series with one value per term and one with any number, truncated at
/// `max_deg`, with `acc` as the working space (it is left empty).
pub fn mul_scalar<R: Ring, K: Key>(
    ring: &R,
    a: &SparsePoly<K, R::Elem>,
    b: &SparsePoly<K, R::Elem>,
    max_deg: u32,
    acc: &mut SlotMap<K, R::Elem>,
) -> SparsePoly<K, R::Elem> {
    debug_assert_eq!((a.width(), acc.width()), (1, b.width()));
    acc.clear();
    mul_into(a, b, max_deg, acc, |va, vb, out| {
        for (o, &v) in out.iter_mut().zip(vb) {
            ring.fma(o, va[0], v);
        }
    });
    acc.take()
}

/// The reciprocal of a series with one value per term, truncated at `max_deg`, with `acc`
/// as the working space (it is left empty).
///
/// The series must start with an invertible constant term, and its other terms must have
/// positive degree. The terms of the result are finalized in order of degree: once every term
/// of lower degree has added its share, $r_C = -c_0^{-1} \sum_{0 < D \le C} c_D\, r_{C - D}$.
pub fn recip<R: Ring, K: Key>(
    ring: &R,
    c: &SparsePoly<K, R::Elem>,
    max_deg: u32,
    acc: &mut SlotMap<K, R::Elem>,
) -> Result<SparsePoly<K, R::Elem>, Error> {
    debug_assert_eq!((c.width(), acc.width()), (1, 1));
    if c.is_empty() || c.keys()[0] != K::default() {
        return Err(Error::NotInvertible);
    }
    let inv0 = ring.inv(c.values(0)[0]).ok_or(Error::NotInvertible)?;
    acc.clear();
    let mut heap = DegreeHeap::new();
    heap.push(0, acc.insert(K::default(), 0).0);
    while let Some(s) = heap.pop() {
        let d = acc.deg(s);
        let r = if d == 0 {
            inv0
        } else {
            ring.neg(ring.mul(acc.values(s)[0], inv0))
        };
        acc.values_mut(s)[0] = r;
        if ring.is_zero(r) {
            continue;
        }
        let k = acc.key(s);
        for y in 1..c.len_upto(max_deg - d) {
            let dy = d + c.degs()[y];
            let (t, created) = acc.insert(k + c.keys()[y], dy);
            if created {
                heap.push(dy, t);
            }
            ring.fma(&mut acc.values_mut(t)[0], c.values(y)[0], r);
        }
    }
    Ok(acc.take())
}

/// Computes $\alpha$ and the instanton series from the period data.
pub fn instanton<R: Ring, K: Key>(
    ring: &R,
    geometry: &Geometry,
    period: &PeriodData<K, R::Elem>,
    max_deg: u32,
) -> Result<InstantonData<K, R::Elem>, Error> {
    let h11 = geometry.h11();
    let m = geometry.n_forms();
    if period.c1.width() != h11 || period.s2.width() != m {
        return Err(Error::DimensionMismatch);
    }
    let zero = R::Elem::default();

    let mut acc = SlotMap::new(1, period.c0.len());
    let c0_inv = recip(ring, &period.c0, max_deg, &mut acc)?;
    let mut acc = SlotMap::new(h11, period.c1.len());
    let alpha = mul_scalar(ring, &c0_inv, &period.c1, max_deg, &mut acc);
    drop(acc);

    // V, with the values of all forms side by side: v[k * h11 + a] = V^(k)_a.
    let half = period.tables.inv[2];
    let forms: Vec<R::Elem> = (0..m)
        .flat_map(|k| {
            geometry
                .form(k)
                .iter()
                .map(|&x| ring.mul(ring.int(x), half))
        })
        .collect();
    let mut v = SparsePoly::new(m * h11);
    let mut vals = vec![zero; m * h11];
    for x in 0..alpha.len() {
        let al = alpha.values(x);
        for (row, out) in forms.chunks_exact(h11).zip(vals.iter_mut()) {
            *out = zero;
            for (&kab, &ab) in row.iter().zip(al) {
                ring.fma(out, kab, ab);
            }
        }
        if vals.iter().any(|&x| !ring.is_zero(x)) {
            v.push(alpha.keys()[x], alpha.degs()[x], &vals);
        }
    }

    let mut acc = SlotMap::new(m, period.s2.len() + alpha.len());
    mul_into(&c0_inv, &period.s2, max_deg, &mut acc, |a, b, out| {
        for (o, &s) in out.iter_mut().zip(b) {
            ring.fma(o, a[0], s);
        }
    });
    mul_into(&alpha, &v, max_deg, &mut acc, |al, vv, out| {
        for (o, vk) in out.iter_mut().zip(vv.chunks_exact(h11)) {
            let mut d = zero;
            for (&x, &y) in al.iter().zip(vk) {
                ring.fma(&mut d, x, y);
            }
            *o = ring.sub(*o, d);
        }
    });
    let inst = acc.take();
    Ok(InstantonData { alpha, inst })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modular::field::{primes, LaneField, Lanes};
    use crate::modular::key::{KeyLayout, PackedKey, PackedLayout};
    use crate::modular::TestRng;
    use std::collections::HashMap;

    type E = Lanes<2>;

    fn ring() -> LaneField<2> {
        let mut ps = primes();
        LaneField::new([ps.next().unwrap(), ps.next().unwrap()])
    }

    /// A random series in the quadrant of two coordinates under the grading (1, 1), with
    /// `width` values per term, and the constant term `c` if given.
    fn random_series(
        rng: &mut TestRng,
        lf: &LaneField<2>,
        layout: &PackedLayout,
        width: usize,
        max_deg: u32,
        c: Option<i64>,
    ) -> SparsePoly<PackedKey, E> {
        let mut map = SlotMap::new(width, 16);
        for _ in 0..40 {
            let x = rng.range(0, max_deg as i64) as i32;
            let y = rng.range(0, (max_deg as i32 - x) as i64) as i32;
            if x + y == 0 {
                continue;
            }
            for v in map.get_or_insert(layout.pack(&[x, y]), (x + y) as u32) {
                *v = lf.int(rng.range(-50, 50));
            }
        }
        if let Some(c) = c {
            map.get_or_insert(PackedKey(0), 0)[0] = lf.int(c);
        }
        map.take()
    }

    /// The product by the definition, as a map from keys to values.
    fn naive_mul(
        lf: &LaneField<2>,
        a: &SparsePoly<PackedKey, E>,
        b: &SparsePoly<PackedKey, E>,
        max_deg: u32,
    ) -> HashMap<PackedKey, Vec<E>> {
        let mut res: HashMap<PackedKey, Vec<E>> = HashMap::new();
        for x in 0..a.len() {
            for y in 0..b.len() {
                if a.degs()[x] + b.degs()[y] > max_deg {
                    continue;
                }
                let e = res
                    .entry(a.keys()[x] + b.keys()[y])
                    .or_insert_with(|| vec![E::default(); b.width()]);
                for (o, &v) in e.iter_mut().zip(b.values(y)) {
                    lf.fma(o, a.values(x)[0], v);
                }
            }
        }
        res.retain(|_, v| v.iter().any(|&x| !lf.is_zero(x)));
        res
    }

    fn as_map(p: &SparsePoly<PackedKey, E>) -> HashMap<PackedKey, Vec<E>> {
        (0..p.len())
            .map(|i| (p.keys()[i], p.values(i).to_vec()))
            .collect()
    }

    #[test]
    fn test_mul_scalar() {
        let lf = ring();
        let layout = PackedLayout::new(&[30, 30]).unwrap();
        let mut rng = TestRng::new(11);
        for (max_deg, width) in [(10, 1), (10, 3), (6, 2), (0, 1)] {
            let a = random_series(&mut rng, &lf, &layout, 1, 10, Some(3));
            let b = random_series(&mut rng, &lf, &layout, width, 10, None);
            let mut acc = SlotMap::new(width, 4);
            let p = mul_scalar(&lf, &a, &b, max_deg, &mut acc);
            assert!(acc.is_empty());
            assert_eq!(as_map(&p), naive_mul(&lf, &a, &b, max_deg));
            assert!(p.degs().iter().all(|&d| d <= max_deg));
        }
    }

    #[test]
    fn test_recip() {
        let lf = ring();
        let layout = PackedLayout::new(&[30, 30]).unwrap();
        let mut rng = TestRng::new(12);
        let max_deg = 12;
        for c in [1, -7, 123_456] {
            let s = random_series(&mut rng, &lf, &layout, 1, max_deg, Some(c));
            let mut acc = SlotMap::new(1, 4);
            let r = recip(&lf, &s, max_deg, &mut acc).unwrap();
            // s * r = 1 up to the maximum degree, and the same the other way around.
            let one = HashMap::from([(PackedKey(0), vec![lf.one()])]);
            assert_eq!(naive_mul(&lf, &s, &r, max_deg), one);
            assert_eq!(naive_mul(&lf, &r, &s, max_deg), one);
            assert_eq!(r.keys()[0], PackedKey(0));
            assert_eq!(r.values(0)[0], lf.inv(lf.int(c)).unwrap());
        }
        // No constant term, or one that is not invertible.
        let mut acc = SlotMap::new(1, 4);
        let s = random_series(&mut rng, &lf, &layout, 1, max_deg, None);
        assert_eq!(recip(&lf, &s, max_deg, &mut acc), Err(Error::NotInvertible));
        let lf7 = LaneField::new([7, 11]);
        let s = random_series(&mut rng, &lf7, &layout, 1, max_deg, Some(14));
        assert_eq!(
            recip(&lf7, &s, max_deg, &mut acc),
            Err(Error::NotInvertible)
        );
    }
}
