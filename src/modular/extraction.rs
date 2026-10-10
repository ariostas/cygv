//! The extraction of the coefficients of the instanton series in the basis of curves.
//!
//! In terms of $q^C = z^C \exp(C \cdot \alpha)$, with $C \cdot \alpha = \sum_a C_a \alpha_a$,
//! the instanton series is $I^{(k)} = \sum_C n^{(k)}_C\, q^C$ over the curve classes of
//! positive degree. For a threefold $n_C = \deg(C)\, A_C$, and for an n-fold
//! $n^{(k)}_C = A_{C,k}$, where $A$ are the invariants before the multicover correction
//! ([`multicover`](super::multicover)). Every term of $\exp(C \cdot \alpha)$ but the first has
//! positive degree, since $\alpha$ has no constant term, so once $n_D\, q^D$ has been
//! subtracted from $I$ for every $D$ of lower degree than $C$, the coefficient of $z^C$ left
//! in $I$ is $n_C$. The curves are therefore taken in order of degree, and each one's
//! coefficient is read off and its $n_C\, q^C$ subtracted, truncated at the maximum degree.
//! What is left in $I$ in the end is $n_C$ at every curve.
//!
//! The exponential is built with the Euler recurrence: with $L = C \cdot \alpha$ and
//! $f = \exp(L)$, applying $\sum_a z_a \partial_{z_a}$ gives
//! $\deg(m)\, f_m = \sum_{u + v = m} \deg(u)\, L_u\, f_v$. The terms of $f$ are therefore
//! final degree by degree, and those of each degree are kept in a table of their own
//! ([`LevelTable`]), which stays small enough to be cached while the terms of lower degree
//! are scattered into it.
//!
//! This stage only takes plain series and tables and returns a series, so that other
//! implementations of it, such as a parallel or a GPU one, can take its place (`PLAN.md`,
//! section 4.1).

use super::field::Ring;
use super::key::KeyLayout;
use super::table::{DegreeHeap, LevelTable, SlotMap, SparsePoly};
use super::Error;

/// Level tables larger than this are freed after each curve rather than kept for the next.
const LEVEL_KEEP: usize = 1 << 15;

/// The working space for the exponentials, kept from one curve to the next.
struct ExpWorkspace<K, E> {
    /// The terms of each degree, from zero to the maximum degree.
    levels: Vec<LevelTable<K, E>>,
    /// The nonzero terms $\deg(u)\, L_u$, in order of degree.
    l_keys: Vec<K>,
    l_vals: Vec<E>,
    /// The degrees present in `l_keys`, each with the range of its terms.
    groups: Vec<(u32, usize, usize)>,
    /// The nonzero terms of the current degree.
    sources: Vec<(K, E)>,
    /// The nonzero coordinates of the current curve, with their images.
    coords: Vec<(usize, E)>,
}

impl<K: super::key::Key, E: Copy + Default + PartialEq> ExpWorkspace<K, E> {
    fn new(sentinel: K, max_deg: u32) -> Self {
        Self {
            levels: (0..=max_deg).map(|_| LevelTable::new(sentinel)).collect(),
            l_keys: Vec::new(),
            l_vals: Vec::new(),
            groups: Vec::new(),
            sources: Vec::new(),
            coords: Vec::new(),
        }
    }

    /// Calls `emit(m, deg(m), f_m)` for every nonzero term of positive degree at most `top`
    /// of $f = \exp(C \cdot \alpha)$, in order of degree.
    ///
    /// `ints` and `inv` hold the images of the integers and their inverses up to `top`.
    #[allow(clippy::too_many_arguments)]
    fn exp<R: Ring<Elem = E>>(
        &mut self,
        ring: &R,
        curve: &[i32],
        alpha: &SparsePoly<K, E>,
        ints: &[E],
        inv: &[E],
        top: u32,
        mut emit: impl FnMut(K, u32, E),
    ) {
        self.coords.clear();
        self.coords.extend(
            curve
                .iter()
                .enumerate()
                .filter(|&(_, &c)| c != 0)
                .map(|(a, &c)| (a, ring.int(c as i64))),
        );
        self.l_keys.clear();
        self.l_vals.clear();
        self.groups.clear();
        for u in 0..alpha.len_upto(top) {
            let al = alpha.values(u);
            let mut v = E::default();
            for &(a, c) in &self.coords {
                ring.fma(&mut v, c, al[a]);
            }
            if ring.is_zero(v) {
                continue;
            }
            let deg = alpha.degs()[u];
            match self.groups.last_mut() {
                Some(g) if g.0 == deg => g.2 += 1,
                _ => self
                    .groups
                    .push((deg, self.l_keys.len(), self.l_keys.len() + 1)),
            }
            self.l_keys.push(alpha.keys()[u]);
            self.l_vals.push(ring.mul(v, ints[deg as usize]));
        }

        *self.levels[0].get_or_insert(K::default()) = ring.one();
        for e in 0..=top {
            let level = &mut self.levels[e as usize];
            if level.is_empty() {
                continue;
            }
            self.sources.clear();
            for (key, val) in level.iter() {
                let f = if e == 0 {
                    val
                } else {
                    ring.mul(val, inv[e as usize])
                };
                if ring.is_zero(f) {
                    continue;
                }
                if e > 0 {
                    emit(key, e, f);
                }
                self.sources.push((key, f));
            }
            if level.capacity() > LEVEL_KEEP {
                level.release();
            } else {
                level.clear();
            }
            for &(g, j0, j1) in &self.groups {
                if g > top - e {
                    break;
                }
                let target = &mut self.levels[(e + g) as usize];
                for &(kb, fb) in &self.sources {
                    for j in j0..j1 {
                        ring.fma(
                            target.get_or_insert(kb + self.l_keys[j]),
                            self.l_vals[j],
                            fb,
                        );
                    }
                }
            }
        }
    }
}

/// Extracts the coefficients $n_C$ of the instanton series `inst` in the basis of the
/// $q^C$, at every curve class of degree at most `max_deg`, as a series with as many values
/// per term as `inst`.
///
/// `alpha` holds the $\alpha_a$, one value per coordinate of `layout`, and must have no
/// constant term; `inv` holds the inverses of the integers up to `max_deg` at least.
pub fn extract<R: Ring, L: KeyLayout>(
    ring: &R,
    layout: &L,
    alpha: &SparsePoly<L::Key, R::Elem>,
    inst: &SparsePoly<L::Key, R::Elem>,
    inv: &[R::Elem],
    max_deg: u32,
) -> Result<SparsePoly<L::Key, R::Elem>, Error> {
    let h11 = layout.dim();
    if alpha.width() != h11 || inv.len() <= max_deg as usize {
        return Err(Error::DimensionMismatch);
    }
    assert!(
        alpha.degs().first().is_none_or(|&d| d > 0),
        "alpha has a constant term"
    );
    let width = inst.width();
    let ints: Vec<R::Elem> = (0..=max_deg).map(|n| ring.int(n as i64)).collect();

    let mut map = SlotMap::new(width, 2 * inst.len());
    let mut heap = DegreeHeap::new();
    for i in 0..inst.len() {
        let (s, _) = map.insert(inst.keys()[i], inst.degs()[i]);
        map.values_mut(s).copy_from_slice(inst.values(i));
        heap.push(inst.degs()[i], s);
    }

    let mut ws = ExpWorkspace::new(layout.sentinel(), max_deg);
    let mut curve = vec![0; h11];
    let mut n_c = vec![R::Elem::default(); width];
    while let Some(slot) = heap.pop() {
        let d = map.deg(slot);
        // A curve of the maximum degree has nothing left to change.
        if d == 0 || d >= max_deg {
            continue;
        }
        n_c.copy_from_slice(map.values(slot));
        if n_c.iter().all(|&x| ring.is_zero(x)) {
            continue;
        }
        let kc = map.key(slot);
        layout.unpack(kc, &mut curve);
        ws.exp(ring, &curve, alpha, &ints, inv, max_deg - d, |m, e, f| {
            let (t, created) = map.insert(kc + m, d + e);
            if created {
                heap.push(d + e, t);
            }
            for (x, &n) in map.values_mut(t).iter_mut().zip(&n_c) {
                *x = ring.sub(*x, ring.mul(n, f));
            }
        });
    }
    Ok(map.take())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modular::field::{primes, LaneField, Lanes};
    use crate::modular::key::{PackedKey, PackedLayout};
    use crate::modular::TestRng;
    use std::collections::HashMap;

    type E = Lanes<2>;

    fn ring() -> LaneField<2> {
        let mut ps = primes();
        LaneField::new([ps.next().unwrap(), ps.next().unwrap()])
    }

    /// $\sum_C n_C\, z^C \exp(C \cdot \alpha)$ by the definition, with the exponential as a
    /// power series of $C \cdot \alpha$, over the quadrant of two coordinates graded by
    /// (1, 1).
    fn naive_series(
        lf: &LaneField<2>,
        layout: &PackedLayout,
        alpha: &SparsePoly<PackedKey, E>,
        n: &HashMap<[i32; 2], Vec<E>>,
        width: usize,
        max_deg: u32,
    ) -> HashMap<PackedKey, Vec<E>> {
        let mut res: HashMap<PackedKey, Vec<E>> = HashMap::new();
        let inv = lf.inverses(max_deg as usize + 1).unwrap();
        for (c, nc) in n {
            let dc = (c[0] + c[1]) as u32;
            // L = C . alpha, and its powers L^j / j!.
            let mut l = HashMap::new();
            for u in 0..alpha.len() {
                let v = lf.add(
                    lf.mul(lf.int(c[0] as i64), alpha.values(u)[0]),
                    lf.mul(lf.int(c[1] as i64), alpha.values(u)[1]),
                );
                l.insert((alpha.keys()[u], alpha.degs()[u]), v);
            }
            let mut power = HashMap::from([((PackedKey(0), 0), lf.one())]);
            let mut exp = power.clone();
            for j in 1..=max_deg {
                let mut next = HashMap::new();
                for (&(kp, dp), &vp) in &power {
                    for (&(kl, dl), &vl) in &l {
                        if dc + dp + dl <= max_deg {
                            let e = next.entry((kp + kl, dp + dl)).or_insert(E::default());
                            lf.fma(e, lf.mul(vp, vl), inv[j as usize]);
                        }
                    }
                }
                for (&k, &v) in &next {
                    let e = exp.entry(k).or_insert(E::default());
                    *e = lf.add(*e, v);
                }
                power = next;
            }
            let kc = layout.pack(c);
            for (&(k, _), &f) in &exp {
                let e = res
                    .entry(kc + k)
                    .or_insert_with(|| vec![E::default(); width]);
                for (x, &y) in e.iter_mut().zip(nc) {
                    lf.fma(x, y, f);
                }
            }
        }
        res.retain(|_, v| v.iter().any(|&x| !lf.is_zero(x)));
        res
    }

    /// Building the series from random coefficients and extracting them gives them back.
    #[test]
    fn test_round_trip() {
        let lf = ring();
        let layout = PackedLayout::new(&[30, 30]).unwrap();
        let mut rng = TestRng::new(21);
        let max_deg = 9;
        let inv = lf.inverses(max_deg as usize + 1).unwrap();
        for width in [1, 3] {
            let mut alpha_map = SlotMap::new(2, 16);
            for _ in 0..12 {
                let x = rng.range(0, 4) as i32;
                let y = rng.range(0, 4 - x as i64) as i32;
                if x + y == 0 {
                    continue;
                }
                for v in alpha_map.get_or_insert(layout.pack(&[x, y]), (x + y) as u32) {
                    *v = lf.int(rng.range(-9, 9));
                }
            }
            let alpha = alpha_map.take();
            let mut n = HashMap::new();
            for _ in 0..15 {
                let x = rng.range(0, max_deg as i64) as i32;
                let y = rng.range(0, max_deg as i64 - x as i64) as i32;
                if x + y == 0 {
                    continue;
                }
                let v: Vec<E> = (0..width).map(|_| lf.int(rng.range(-50, 50))).collect();
                n.insert([x, y], v);
            }
            let series = naive_series(&lf, &layout, &alpha, &n, width, max_deg);
            let mut inst_map = SlotMap::new(width, 16);
            for (&k, v) in &series {
                let mut c = [0; 2];
                layout.unpack(k, &mut c);
                inst_map
                    .get_or_insert(k, (c[0] + c[1]) as u32)
                    .copy_from_slice(v);
            }
            let inst = inst_map.take();
            let got = extract(&lf, &layout, &alpha, &inst, &inv, max_deg).unwrap();
            let got: HashMap<PackedKey, Vec<E>> = (0..got.len())
                .map(|i| (got.keys()[i], got.values(i).to_vec()))
                .collect();
            let expected: HashMap<PackedKey, Vec<E>> =
                n.iter().map(|(c, v)| (layout.pack(c), v.clone())).collect();
            assert_eq!(got, expected, "width {width}");
        }
    }

    #[test]
    fn test_errors() {
        let lf = ring();
        let layout = PackedLayout::new(&[30, 30]).unwrap();
        let inv = lf.inverses(5).unwrap();
        let empty = SparsePoly::new(2);
        let inst = SparsePoly::new(1);
        assert_eq!(
            extract(&lf, &layout, &empty, &inst, &inv, 5),
            Err(Error::DimensionMismatch)
        );
        assert_eq!(
            extract(&lf, &layout, &SparsePoly::new(3), &inst, &inv, 4),
            Err(Error::DimensionMismatch)
        );
        assert!(extract(&lf, &layout, &empty, &inst, &inv, 4)
            .unwrap()
            .is_empty());
    }
}
