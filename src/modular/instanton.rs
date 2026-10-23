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
//!
//! The products run in parallel over chunks of their first factor, each accumulating in a map
//! of its own that it then adds to a [`ShardedMap`]. Since each $\alpha_a$ has few terms,
//! $\sum_a \alpha_a V^{(k)}_a = \frac12 \sum_{a, b} K^{(k)}_{ab}\, \alpha_a \alpha_b$ is taken
//! as products of the columns of $\alpha$, only for the pairs $a \le b$ on which some form
//! does not vanish, and by symmetry each pair of terms only once.

use super::field::Ring;
use super::key::Key;
use super::period::PeriodData;
use super::table::{DegreeHeap, ShardedMap, SlotMap, SparsePoly};
use super::{Error, Geometry};
use core::ops::{Range, RangeInclusive};
use rayon::prelude::*;

/// The output of the instanton stage.
#[derive(Debug, Clone)]
pub struct InstantonData<K, E> {
    /// $\alpha_a$, $h^{1,1}$ values per term.
    pub alpha: SparsePoly<K, E>,
    /// $I^{(k)}$, one value per bilinear form.
    pub inst: SparsePoly<K, E>,
    /// $\beta = c_1 / c_0$ along the rows of the charges, if the period has $c_1$ along them,
    /// so that $\alpha_a = \sum_j Q^0_{ja} \beta_j + \sum_r Q_{ra} \beta_r$.
    pub beta: Option<SparsePoly<K, E>>,
}

/// Splits the indices `0..n` into at most `n_chunks` consecutive ranges of similar total
/// `cost`.
fn chunks(n: usize, n_chunks: usize, cost: impl Fn(usize) -> usize) -> Vec<Range<usize>> {
    let total: usize = (0..n).map(&cost).sum();
    let per_chunk = total.div_ceil(n_chunks.max(1)).max(1);
    let mut res = Vec::new();
    let (mut start, mut acc) = (0, 0);
    for x in 0..n {
        acc += cost(x);
        if acc >= per_chunk || x + 1 == n {
            res.push(start..x + 1);
            (start, acc) = (x + 1, 0);
        }
    }
    res
}

/// The most terms a chunk of [`par_accumulate`] collects before adding them, when the sums
/// may be taken in any order.
const LOCAL_TERMS: usize = 1 << 16;

/// The number of chunks of [`par_accumulate`] when the sums must be taken in a fixed order.
const ORDERED_CHUNKS: usize = 128;

/// Calls `body(x, map)` for every `x` in `0..n`, in parallel over chunks of similar total
/// `cost` on the threads of the current rayon pool, where `map` is a map with `width` values
/// per slot of the chunk's own. Each chunk then adds its map to `acc` with `add`.
///
/// Unless `exact`, the sums are taken in a fixed order: the chunks do not depend on the
/// number of threads, and their maps are added in order, a round of one chunk per thread at a
/// time.
fn par_accumulate<K: Key, E: Copy + Default + PartialEq + Send + Sync>(
    acc: &ShardedMap<K, E>,
    width: usize,
    n: usize,
    cost: impl Fn(usize) -> usize,
    body: impl Fn(usize, &mut SlotMap<K, E>) + Sync,
    add: impl Fn(&mut [E], &[E]) + Sync,
    exact: bool,
) {
    let threads = rayon::current_num_threads();
    if !exact {
        let chunks = chunks(n, ORDERED_CHUNKS, cost);
        for round in chunks.chunks(threads) {
            let mut maps: Vec<SlotMap<K, E>> = round
                .par_iter()
                .map(|xs| {
                    let mut map = SlotMap::new(width, 1024);
                    for x in xs.clone() {
                        body(x, &mut map);
                    }
                    map
                })
                .collect();
            acc.add_in_order(&mut maps, &add);
        }
        return;
    }
    let n_chunks = if threads == 1 { 1 } else { 4 * threads };
    chunks(n, n_chunks, cost).into_par_iter().for_each_init(
        || SlotMap::new(width, 1024),
        |local, xs| {
            for x in xs {
                body(x, local);
                // A chunk's terms can be most of the result: added in batches, they take
                // a bounded share of the memory.
                if local.len() >= LOCAL_TERMS {
                    acc.add_from(local, &add);
                }
            }
            acc.add_from(local, &add);
        },
    );
}

/// The terms of `p` of degree below `deg`, i.e. the start of the others.
#[inline]
fn len_below<K: Key, E: Copy + Default + PartialEq>(p: &SparsePoly<K, E>, deg: u32) -> usize {
    deg.checked_sub(1).map_or(0, |d| p.len_upto(d))
}

/// Adds the terms of $a \cdot b$ of degree in `degs` to `acc`, where `a` has one value per
/// term and `b` as many as `acc`.
fn mul_scalar_into<R: Ring, K: Key>(
    ring: &R,
    a: &SparsePoly<K, R::Elem>,
    b: &SparsePoly<K, R::Elem>,
    degs: RangeInclusive<u32>,
    acc: &ShardedMap<K, R::Elem>,
) {
    debug_assert_eq!(a.width(), 1);
    let (lo, hi) = (*degs.start(), *degs.end());
    let inner = |x: usize| {
        let dx = a.degs()[x];
        len_below(b, lo.saturating_sub(dx))..b.len_upto(hi - dx)
    };
    par_accumulate(
        acc,
        b.width(),
        a.len_upto(hi),
        |x| inner(x).len(),
        |x, local| {
            let (kx, dx, ax) = (a.keys()[x], a.degs()[x], a.values(x)[0]);
            for y in inner(x) {
                let out = local.get_or_insert(kx + b.keys()[y], dx + b.degs()[y]);
                for (o, &v) in out.iter_mut().zip(b.values(y)) {
                    ring.fma(o, ax, v);
                }
            }
        },
        |out, v| add_into(ring, out, v),
        R::EXACT,
    );
}

/// The number of shards for a table that the threads of the current rayon pool add to.
fn n_shards() -> usize {
    match rayon::current_num_threads() {
        1 => 1,
        n => 4 * n,
    }
}

/// The product of a series with one value per term and one with any number, truncated at
/// `max_deg`, computed on the threads of the current rayon pool.
pub fn mul_scalar<R: Ring, K: Key>(
    ring: &R,
    a: &SparsePoly<K, R::Elem>,
    b: &SparsePoly<K, R::Elem>,
    max_deg: u32,
) -> SparsePoly<K, R::Elem> {
    let acc = ShardedMap::new(b.width(), n_shards(), b.len());
    mul_scalar_into(ring, a, b, 0..=max_deg, &acc);
    acc.take()
}

/// Adds `v` to `out`, value by value.
#[inline]
fn add_into<R: Ring>(ring: &R, out: &mut [R::Elem], v: &[R::Elem]) {
    for (o, &x) in out.iter_mut().zip(v) {
        *o = ring.add(*o, x);
    }
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

/// About the number of blocks of degrees that the instanton series is taken in.
const INSTANTON_BLOCKS: usize = 8;

/// Splits the degrees `0..=max_deg` into ranges `(lo, hi)` in which `p` has about
/// $1/n$ of its terms each, or the terms of one degree when they are more.
fn degree_blocks<K: Key, E: Copy + Default + PartialEq>(
    p: &SparsePoly<K, E>,
    max_deg: u32,
    n: usize,
) -> Vec<(u32, u32)> {
    let per_block = p.len_upto(max_deg).div_ceil(n.max(1)).max(1);
    let mut blocks = Vec::new();
    let mut lo = 0;
    for d in 0..=max_deg {
        if d == max_deg || p.len_upto(d) - len_below(p, lo) >= per_block {
            blocks.push((lo, d));
            lo = d + 1;
        }
    }
    blocks
}

/// Computes $\alpha$ and the instanton series from the period data.
pub fn instanton<R: Ring, K: Key>(
    ring: &R,
    geometry: &Geometry,
    period: &mut PeriodData<K, R::Elem>,
    max_deg: u32,
) -> Result<InstantonData<K, R::Elem>, Error> {
    let h11 = geometry.h11();
    let m = geometry.n_forms();
    if period.c1.width() != h11 || period.s2.width() != m {
        return Err(Error::DimensionMismatch);
    }
    let mut acc = SlotMap::new(1, period.c0.len());
    let c0_inv = recip(ring, &period.c0, max_deg, &mut acc)?;
    drop(acc);
    // Each part of the period is freed once it is used, which lowers the peak.
    let alpha = mul_scalar(ring, &c0_inv, &period.c1, max_deg);
    period.c1 = SparsePoly::new(h11);
    let beta = period
        .c1_rows
        .take()
        .map(|c1| mul_scalar(ring, &c0_inv, &c1, max_deg));

    // sum_a alpha_a V^(k)_a = 1/2 sum_{a, b} K^(k)_ab alpha_a alpha_b, by pairs of columns
    // a <= b with some K_ab nonzero, since each alpha_a has few terms: the product of the two
    // columns times K_ab for a < b, and for a = b half of it, i.e. the products of the terms
    // x < y and half of those of x = y.
    // The terms of each column, as indices into alpha.
    let mut columns: Vec<Vec<u32>> = vec![Vec::new(); h11];
    for x in 0..alpha.len() {
        for (a, &v) in alpha.values(x).iter().enumerate() {
            if !ring.is_zero(v) {
                columns[a].push(x as u32);
            }
        }
    }
    let term = |a: usize, i: usize| {
        let x = columns[a][i] as usize;
        (alpha.keys()[x], alpha.degs()[x], alpha.values(x)[a])
    };
    let half = period.tables.inv[2];
    // For each pair of columns, the coefficients K^(k)_ab of the forms.
    let mut pairs: Vec<(usize, usize, Vec<R::Elem>)> = Vec::new();
    for a in 0..h11 {
        for b in a..h11 {
            let kab: Vec<i64> = (0..m).map(|k| geometry.form(k)[a * h11 + b]).collect();
            if kab.iter().any(|&x| x != 0) && !columns[a].is_empty() && !columns[b].is_empty() {
                pairs.push((a, b, kab.iter().map(|&x| ring.int(x)).collect()));
            }
        }
    }
    // The outer terms: each pair with each term of its first column.
    let outer: Vec<(u32, u32)> = pairs
        .iter()
        .enumerate()
        .flat_map(|(p, &(a, _, _))| (0..columns[a].len() as u32).map(move |i| (p as u32, i)))
        .collect();
    // The inner terms of an outer one whose products have degree in lo..=hi.
    let inner = |o: usize, lo: u32, hi: u32| {
        let (p, i) = (outer[o].0 as usize, outer[o].1 as usize);
        let (a, b, _) = &pairs[p];
        let d = term(*a, i).1;
        if d > hi {
            return 0..0;
        }
        let column = &columns[*b];
        let deg = |y: &u32| alpha.degs()[*y as usize];
        let mut start = column.partition_point(|y| deg(y) < lo.saturating_sub(d));
        if a == b {
            start = start.max(i);
        }
        let end = column.partition_point(|y| deg(y) <= hi - d);
        start..end.max(start)
    };

    // The instanton series is taken a block of degrees at a time, so that only one block is
    // accumulated at once, in tables twice to four times the size of the terms.
    let blocks = degree_blocks(&period.s2, max_deg, INSTANTON_BLOCKS);
    let mut parts = Vec::with_capacity(blocks.len());
    for (lo, hi) in blocks {
        let terms = |p: &SparsePoly<K, R::Elem>| p.len_upto(hi) - len_below(p, lo);
        let acc = ShardedMap::new(m, n_shards(), terms(&period.s2) + terms(&alpha));
        mul_scalar_into(ring, &c0_inv, &period.s2, lo..=hi, &acc);
        par_accumulate(
            &acc,
            m,
            outer.len(),
            |o| inner(o, lo, hi).len(),
            |o, local| {
                let (p, i) = (outer[o].0 as usize, outer[o].1 as usize);
                let (a, b, kab) = &pairs[p];
                let (kx, dx, vx) = term(*a, i);
                for j in inner(o, lo, hi) {
                    let (ky, dy, vy) = term(*b, j);
                    let mut v = ring.mul(vx, vy);
                    if a == b && i == j {
                        v = ring.mul(v, half);
                    }
                    for (out, &k) in local.get_or_insert(kx + ky, dx + dy).iter_mut().zip(kab) {
                        *out = ring.sub(*out, ring.mul(k, v));
                    }
                }
            },
            |out, v| add_into(ring, out, v),
            R::EXACT,
        );
        parts.push(acc.take());
    }
    period.s2 = SparsePoly::new(m);
    let inst = SparsePoly::concat(m, parts);
    Ok(InstantonData { alpha, inst, beta })
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
        let layout = PackedLayout::<u128>::new(&[30, 30]).unwrap();
        let mut rng = TestRng::new(11);
        for (max_deg, width) in [(10, 1), (10, 3), (6, 2), (0, 1)] {
            let a = random_series(&mut rng, &lf, &layout, 1, 10, Some(3));
            let b = random_series(&mut rng, &lf, &layout, width, 10, None);
            let p = mul_scalar(&lf, &a, &b, max_deg);
            assert_eq!(as_map(&p), naive_mul(&lf, &a, &b, max_deg));
            assert!(p.degs().iter().all(|&d| d <= max_deg));
        }
    }

    #[test]
    fn test_recip() {
        let lf = ring();
        let layout = PackedLayout::<u128>::new(&[30, 30]).unwrap();
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
