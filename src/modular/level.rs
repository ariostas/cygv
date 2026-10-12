//! The level tables of the extraction's exponentials, keyed by linear hashes of offsets.
//!
//! The terms of $\exp(C \cdot \alpha)$ sit at offsets $m$ from $C$ that are sums of terms of
//! $\alpha$, and the extraction builds them one degree at a time, each degree in a
//! [`LevelTable`] of its own. Those tables are where the pipeline spends most of its time, so
//! they do not key the offsets by their curve class keys. They key each one by
//! $H(m) = \varphi(m)\, M \bmod 2^b$ instead, with $\varphi$ packing the coordinates of $m$
//! into a $b$-bit [`LevelWord`] (see [`LevelLayout`]) and $M$ an odd constant. Both maps are
//! linear, so the entry that a term $u$ times a source $v$ adds to is
//! $H(u + v) = H(u) + H(v)$, one addition with no hash to compute; $H$ is a bijection, so
//! comparing hashes compares offsets; and its high bits are a multiplicative hash, which
//! picks the bucket. An offset is only decoded, by multiplying by $M^{-1}$ and unpacking,
//! for the terms the extraction emits.
//!
//! $\varphi$ packs a box that holds every offset, and leaves out a coordinate that the degree
//! determines, which usually makes it fit 64 bits. The layout and the tables follow the C
//! program's (`PLAN.md`, PR #102).

use super::key::Key;
use core::cmp::Reverse;
use core::fmt::Debug;

/// A word that offsets are packed into: `u64`, or `u128` when the box needs more bits.
pub trait LevelWord: Copy + Eq + Default + Debug + Send + Sync + 'static {
    /// The number of bits.
    const BITS: u32;

    /// The odd multiplier of the hash.
    const MUL: Self;

    /// Its inverse modulo $2^b$.
    const MUL_INV: Self;

    fn wadd(self, other: Self) -> Self;

    fn wsub(self, other: Self) -> Self;

    fn wmul(self, other: Self) -> Self;

    fn xor(self, other: Self) -> Self;

    /// $x \bmod 2^b$, i.e. `x` sign-extended.
    fn from_i64(x: i64) -> Self;

    /// $2^s$, for $s < b$.
    fn pow2(s: u32) -> Self;

    /// The `width` bits from bit `offset` on, for `width` at most 64.
    fn field(self, offset: u32, width: u32) -> u64;

    /// The top `bits` bits, for `bits` between 1 and 32.
    fn top(self, bits: u32) -> usize;
}

/// The inverse of an odd number modulo $2^{64}$ or $2^{128}$, by Newton's iteration, which
/// doubles the number of correct low bits each time, starting from three ($m m \equiv 1
/// \pmod 8$).
macro_rules! level_word {
    ($t:ty, $mul:expr) => {
        impl LevelWord for $t {
            const BITS: u32 = <$t>::BITS;
            const MUL: Self = $mul;
            const MUL_INV: Self = {
                let mut x: $t = $mul;
                let mut i = 0;
                while i < 7 {
                    x = x.wrapping_mul((2 as $t).wrapping_sub(($mul as $t).wrapping_mul(x)));
                    i += 1;
                }
                x
            };

            #[inline]
            fn wadd(self, other: Self) -> Self {
                self.wrapping_add(other)
            }

            #[inline]
            fn wsub(self, other: Self) -> Self {
                self.wrapping_sub(other)
            }

            #[inline]
            fn wmul(self, other: Self) -> Self {
                self.wrapping_mul(other)
            }

            #[inline]
            fn xor(self, other: Self) -> Self {
                self ^ other
            }

            #[inline]
            fn from_i64(x: i64) -> Self {
                x as $t
            }

            #[inline]
            fn pow2(s: u32) -> Self {
                1 << s
            }

            #[inline]
            fn field(self, offset: u32, width: u32) -> u64 {
                let x = (self >> offset) as u64;
                if width >= 64 {
                    x
                } else {
                    x & ((1 << width) - 1)
                }
            }

            #[inline]
            fn top(self, bits: u32) -> usize {
                (self >> (Self::BITS - bits)) as usize
            }
        }
    };
}

level_word!(u64, 0x9E37_79B9_7F4A_7C15);
level_word!(u128, 0x9E37_79B9_7F4A_7C15_D6E8_FEB8_6659_FD93);

/// The box of [`level_box`]: per coordinate the lowest value and the width of its field, and
/// the coordinate that the degree determines, which gets no field.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LevelBox {
    lo: Vec<i64>,
    widths: Vec<u32>,
    dropped: Option<usize>,
}

impl LevelBox {
    /// The number of bits of all fields together.
    pub fn bits(&self) -> u32 {
        self.widths.iter().sum()
    }
}

/// The box that holds every sum of `terms` of degree at most `max_deg`, the terms given by
/// their coordinates and degree, which must be positive, under `grading`.
///
/// A sum $x = \sum_u n_u u$ of degree $e = \sum_u n_u \deg(u) \le D$ has
/// $x_t = \sum_u n_u \deg(u)\, (u_t / \deg(u))$, which lies between $e$ times the smallest
/// and the largest of the ratios $u_t / \deg(u)$, and therefore between $D$ times them, or
/// zero, rounded outwards. One coordinate with a nonzero weight is left out, the widest:
/// $x_{t_0} = (e - \sum_{t \ne t_0} w_t x_t) / w_{t_0}$.
pub fn level_box<'a>(
    terms: impl IntoIterator<Item = (&'a [i32], u32)>,
    grading: &[i32],
    max_deg: u32,
) -> LevelBox {
    let dim = grading.len();
    let (mut lo, mut hi) = (vec![0i64; dim], vec![0i64; dim]);
    let max_deg = max_deg as i64;
    for (u, deg) in terms {
        assert!(deg > 0, "the terms of a level box need a positive degree");
        let deg = deg as i64;
        for t in 0..dim {
            let x = max_deg * u[t] as i64;
            lo[t] = lo[t].min(x.div_euclid(deg));
            hi[t] = hi[t].max(-((-x).div_euclid(deg)));
        }
    }
    let mut widths: Vec<u32> = lo
        .iter()
        .zip(&hi)
        .map(|(&l, &h)| u64::BITS - ((h - l) as u64).leading_zeros())
        .collect();
    let dropped = (0..dim)
        .filter(|&t| grading[t] != 0)
        .max_by_key(|&t| (widths[t], Reverse(t)));
    if let Some(t) = dropped {
        widths[t] = 0;
    }
    LevelBox {
        lo,
        widths,
        dropped,
    }
}

/// How offsets are packed into a [`LevelWord`]: the coordinates of a [`LevelBox`] but the
/// dropped one, each $x_t - \ell_t$ in a field of its own, so that
/// $\varphi(x) = \sum_t x_t 2^{o_t} \bmod 2^b$, which is linear, and
/// $\varphi(x) - \varphi(\ell)$ holds the fields as digits.
#[derive(Debug, Clone)]
pub struct LevelLayout<W> {
    lo: Vec<i64>,
    offsets: Vec<u32>,
    widths: Vec<u32>,
    grading: Vec<i64>,
    dropped: Option<usize>,
    /// $\varphi(\ell)$.
    base: W,
    /// The hash of $\varphi(\ell) + 2^B$, with $B$ the bits of the fields, which is no
    /// offset's: the tables store hashes XORed with it, so that zero marks an empty entry.
    marker: W,
}

impl<W: LevelWord> LevelLayout<W> {
    /// The layout of `level_box` under `grading`, or `None` if its fields do not fit in
    /// $b - 1$ bits, which leaves room for the marker.
    pub fn new(level_box: &LevelBox, grading: &[i32]) -> Option<Self> {
        let bits = level_box.bits();
        if bits >= W::BITS {
            return None;
        }
        let mut offsets = Vec::with_capacity(grading.len());
        let mut o = 0;
        for &w in &level_box.widths {
            offsets.push(o);
            o += w;
        }
        let mut layout = Self {
            lo: level_box.lo.clone(),
            offsets,
            widths: level_box.widths.clone(),
            grading: grading.iter().map(|&w| w as i64).collect(),
            dropped: level_box.dropped,
            base: W::default(),
            marker: W::default(),
        };
        layout.base = layout.phi(level_box.lo.iter().copied());
        layout.marker = layout.base.wadd(W::pow2(bits)).wmul(W::MUL);
        Some(layout)
    }

    fn phi(&self, coords: impl IntoIterator<Item = i64>) -> W {
        let mut p = W::default();
        for (t, x) in coords.into_iter().enumerate() {
            if self.widths[t] > 0 {
                p = p.wadd(W::from_i64(x).wmul(W::pow2(self.offsets[t])));
            }
        }
        p
    }

    /// The hash of the offset with coordinates `coords`, which must lie in the box.
    #[inline]
    pub fn hash(&self, coords: &[i32]) -> W {
        self.phi(coords.iter().map(|&x| x as i64)).wmul(W::MUL)
    }

    /// The marker, see [`LevelTable::new`].
    #[inline]
    pub fn marker(&self) -> W {
        self.marker
    }

    /// Writes the coordinates of the offset of degree `deg` whose hash is `h` into `coords`.
    #[inline]
    pub fn decode(&self, h: W, deg: u32, coords: &mut [i32]) {
        let y = h.wmul(W::MUL_INV).wsub(self.base);
        let mut rest = deg as i64;
        for (t, c) in coords.iter_mut().enumerate() {
            let x = self.lo[t] + y.field(self.offsets[t], self.widths[t]) as i64;
            *c = x as i32;
            rest -= self.grading[t] * x;
        }
        if let Some(t) = self.dropped {
            // Its field is empty, so the loop counted it as its lowest value.
            rest += self.grading[t] * self.lo[t];
            debug_assert_eq!(rest % self.grading[t], 0, "the degree does not fit");
            coords[t] = (rest / self.grading[t]) as i32;
        }
    }
}

/// The inverse of an odd number modulo $2^{64}$, by Newton's iteration.
fn inverse_mod_2_64(m: i64) -> i64 {
    let mut x = m;
    for _ in 0..6 {
        x = x.wrapping_mul(2i64.wrapping_sub(m.wrapping_mul(x)));
    }
    x
}

impl<W: LevelWord> LevelLayout<W> {
    /// What [`OffsetKeys::key`] needs to turn hashes into keys, with `units` the keys of the
    /// unit vectors, or anything for a coordinate that is always zero.
    pub fn offset_keys<K: Key>(&self, units: &[K]) -> OffsetKeys<W, K> {
        let mut fields = Vec::new();
        let (mut key_lo, mut weight_lo) = (K::default(), 0i64);
        for (t, &unit) in units.iter().enumerate() {
            if Some(t) == self.dropped {
                continue;
            }
            key_lo = key_lo + unit.times(self.lo[t]);
            weight_lo += self.grading[t] * self.lo[t];
            if self.widths[t] > 0 {
                fields.push((self.offsets[t], self.widths[t], unit, self.grading[t]));
            }
        }
        let dropped = self.dropped.map(|t| {
            let w = self.grading[t];
            let shift = w.trailing_zeros();
            (units[t], shift, inverse_mod_2_64(w >> shift))
        });
        OffsetKeys {
            fields,
            base: self.base,
            key_lo,
            weight_lo,
            dropped,
        }
    }
}

/// Turns the hashes of a [`LevelLayout`] into the keys of the offsets, without unpacking the
/// coordinates: a key is linear in them, and so is the degree, which gives the dropped one.
#[derive(Debug, Clone)]
pub struct OffsetKeys<W, K> {
    /// Per field: its offset and width, the key of the unit vector and the weight.
    fields: Vec<(u32, u32, K, i64)>,
    /// $\varphi(\ell)$.
    base: W,
    /// The key and degree of $\ell$ without the dropped coordinate.
    key_lo: K,
    weight_lo: i64,
    /// For the dropped coordinate: the key of its unit vector, and its weight $2^s r$ with $r$
    /// odd, as $s$ and $r^{-1} \bmod 2^{64}$, which divide exactly.
    dropped: Option<(K, u32, i64)>,
}

impl<W: LevelWord, K: Key> OffsetKeys<W, K> {
    /// The key of the offset of degree `deg` whose hash is `h`.
    #[inline]
    pub fn key(&self, h: W, deg: u32) -> K {
        let y = h.wmul(W::MUL_INV).wsub(self.base);
        let (mut key, mut weight) = (self.key_lo, self.weight_lo);
        for &(offset, width, unit, w) in &self.fields {
            let d = y.field(offset, width) as i64;
            key = key + unit.times(d);
            weight += w * d;
        }
        if let Some((unit, shift, inv)) = self.dropped {
            let x = ((deg as i64 - weight) >> shift).wrapping_mul(inv);
            key = key + unit.times(x);
        }
        key
    }
}

/// The capacity of a [`LevelTable`] when it is first used.
const LEVEL_MIN_CAPACITY: usize = 256;

/// An open-addressing table of the terms of one degree of an exponential being built, keyed
/// by the hashes of their offsets.
///
/// The value is stored next to the hash, so that a probe touches both, and the entries keep
/// no degree, since all of them have the same one. Each entry stores its hash XORed with the
/// layout's marker, which no offset hashes to, so that zero marks an empty entry and the
/// table can be allocated zeroed. The bucket comes from the high bits of the hash. The table
/// remembers which entries are in use, in order of insertion, so that iterating over and
/// clearing it cost the number of entries rather than the capacity. Memory is only
/// allocated on first use, and the table grows when it is three quarters full.
#[derive(Debug, Clone)]
pub struct LevelTable<W, V> {
    entries: Vec<(W, V)>,
    /// The indices of the entries in use, in order of insertion.
    used: Vec<u32>,
    /// The base-2 logarithm of the capacity, when there is one.
    log2: u32,
    marker: W,
}

impl<W: LevelWord, V: Copy + Default> LevelTable<W, V> {
    /// An empty table for the hashes of a layout with this `marker`
    /// ([`LevelLayout::marker`]).
    pub fn new(marker: W) -> Self {
        Self {
            entries: Vec::new(),
            used: Vec::new(),
            log2: 0,
            marker,
        }
    }

    /// The number of entries in use.
    #[inline]
    pub fn len(&self) -> usize {
        self.used.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.used.is_empty()
    }

    /// The number of entries allocated.
    #[inline]
    pub fn capacity(&self) -> usize {
        self.entries.len()
    }

    /// The entry holding hash `h`, or else the empty entry where it would go.
    #[inline]
    fn probe(&self, h: W) -> Result<usize, usize> {
        let mask = self.entries.len() - 1;
        let x = h.xor(self.marker);
        let mut i = h.top(self.log2);
        loop {
            let k = self.entries[i].0;
            if k == x {
                return Ok(i);
            }
            if k == W::default() {
                return Err(i);
            }
            i = (i + 1) & mask;
        }
    }

    /// The value of hash `h`, if it has an entry.
    #[inline]
    pub fn get(&self, h: W) -> Option<&V> {
        if self.entries.is_empty() {
            return None;
        }
        self.probe(h).ok().map(|i| &self.entries[i].1)
    }

    /// The value of hash `h`, inserted as zero if it has no entry.
    #[inline]
    pub fn get_or_insert(&mut self, h: W) -> &mut V {
        debug_assert!(h != self.marker, "the marker was inserted");
        if self.entries.is_empty() {
            self.entries = vec![(W::default(), V::default()); LEVEL_MIN_CAPACITY];
            self.log2 = LEVEL_MIN_CAPACITY.trailing_zeros();
        }
        let i = match self.probe(h) {
            Ok(i) => i,
            Err(i) => {
                let i = if 4 * (self.used.len() + 1) > 3 * self.entries.len() {
                    self.grow();
                    self.probe(h).unwrap_err()
                } else {
                    i
                };
                self.entries[i] = (h.xor(self.marker), V::default());
                self.used.push(i as u32);
                i
            }
        };
        &mut self.entries[i].1
    }

    /// Doubles the capacity.
    fn grow(&mut self) {
        let capacity = 2 * self.entries.len();
        assert!(capacity <= 1 << 32, "too many entries for a level table");
        let old = core::mem::replace(
            &mut self.entries,
            vec![(W::default(), V::default()); capacity],
        );
        self.log2 += 1;
        for i in self.used.iter_mut() {
            let entry = old[*i as usize];
            let mut j = entry.0.xor(self.marker).top(self.log2);
            while self.entries[j].0 != W::default() {
                j = (j + 1) & (capacity - 1);
            }
            self.entries[j] = entry;
            *i = j as u32;
        }
    }

    /// Asks the processor to fetch the entry where hash `h` would be probed first into the
    /// cache, ahead of a lookup a few iterations later. It does nothing on other
    /// architectures than x86-64, or before the table has memory.
    #[inline]
    pub fn prefetch(&self, h: W) {
        #[cfg(target_arch = "x86_64")]
        if !self.entries.is_empty() {
            let ptr = self.entries[h.top(self.log2)..].as_ptr().cast::<i8>();
            // SAFETY: a prefetch has no effect on memory, and SSE is part of x86-64.
            #[allow(unused_unsafe)]
            unsafe {
                core::arch::x86_64::_mm_prefetch(ptr, core::arch::x86_64::_MM_HINT_T0);
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        let _ = h;
    }

    /// The entries in use, as hashes and values, in order of insertion.
    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = (W, V)> + '_ {
        self.used.iter().map(|&i| {
            let (k, v) = self.entries[i as usize];
            (k.xor(self.marker), v)
        })
    }

    /// Removes all entries, keeping the memory.
    pub fn clear(&mut self) {
        for &i in &self.used {
            self.entries[i as usize].0 = W::default();
        }
        self.used.clear();
    }

    /// Removes all entries and frees the memory.
    pub fn release(&mut self) {
        self.entries = Vec::new();
        self.used = Vec::new();
        self.log2 = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modular::key::{KeyLayout, PackedKey, PackedLayout};
    use crate::modular::TestRng;
    use std::collections::HashMap;

    #[test]
    fn test_inverses() {
        assert_eq!(u64::MUL.wmul(u64::MUL_INV), 1);
        assert_eq!(u128::MUL.wmul(u128::MUL_INV), 1);
        for m in [1i64, 3, -5, 12345677, i64::MAX, -1] {
            assert_eq!(m.wrapping_mul(inverse_mod_2_64(m)), 1);
        }
    }

    /// Random terms of positive degree under `grading`, in the quadrant-like region where
    /// the grading is positive.
    fn terms(rng: &mut TestRng, grading: &[i32], n: usize) -> Vec<(Vec<i32>, u32)> {
        let mut res = Vec::new();
        while res.len() < n {
            let u: Vec<i32> = grading.iter().map(|_| rng.range(-3, 6) as i32).collect();
            let d: i64 = u
                .iter()
                .zip(grading)
                .map(|(&x, &w)| x as i64 * w as i64)
                .sum();
            if d > 0 {
                res.push((u, d as u32));
            }
        }
        res
    }

    /// Hashes add like the offsets, and decode to them, for sums of terms up to the degree.
    fn check_layout<W: LevelWord>(grading: &[i32], seed: u64) {
        let mut rng = TestRng::new(seed);
        let ts = terms(&mut rng, grading, 12);
        let max_deg = 40;
        let b = level_box(ts.iter().map(|(u, d)| (u.as_slice(), *d)), grading, max_deg);
        let Some(layout) = LevelLayout::<W>::new(&b, grading) else {
            assert!(b.bits() >= W::BITS);
            return;
        };
        let dim = grading.len();
        let mut coords = vec![0; dim];
        // Keys of a packed layout wide enough for the box.
        let packed = PackedLayout::new(&vec![1 << 12; dim]).unwrap();
        let units: Vec<PackedKey> = (0..dim)
            .map(|t| {
                let mut e = vec![0; dim];
                e[t] = 1;
                packed.pack(&e)
            })
            .collect();
        let keys = layout.offset_keys(&units);
        for _ in 0..300 {
            // A random sum of terms of degree at most max_deg, built term by term.
            let (mut x, mut h, mut deg) = (vec![0i32; dim], W::default(), 0);
            loop {
                let (u, d) = &ts[rng.range(0, ts.len() as i64 - 1) as usize];
                if deg + d > max_deg {
                    break;
                }
                deg += d;
                for (a, b) in x.iter_mut().zip(u) {
                    *a += b;
                }
                h = h.wadd(layout.hash(u));
                assert_eq!(h, layout.hash(&x));
                assert_ne!(h, layout.marker());
                layout.decode(h, deg, &mut coords);
                assert_eq!(coords, x);
                assert_eq!(keys.key(h, deg), packed.pack(&x));
            }
        }
    }

    #[test]
    fn test_layout() {
        for (i, grading) in [
            vec![1, 1],
            vec![1, 2, 3],
            vec![6, -12, -14, -17, -18, -4, -5, -12],
            vec![0, 1, 0, 2, 5],
            vec![-1, 3],
        ]
        .iter()
        .enumerate()
        {
            check_layout::<u64>(grading, i as u64);
            check_layout::<u128>(grading, i as u64);
        }
    }

    #[test]
    fn test_box() {
        // Terms (1, 0), (0, 1), (1, -1) under (1, 2): the ratios are 1, 0, 1 and 0, 1/2,
        // -1, so at degree 4 the box is [0, 4] x [-4, 2].
        let ts = [(vec![1, 0], 1), (vec![0, 1], 2), (vec![1, -1], 1)];
        let b = level_box(ts.iter().map(|(u, d)| (u.as_slice(), *d)), &[1, 2], 4);
        assert_eq!(b.lo, [0, -4]);
        // Widths 3 (0..=4) and 3 (-4..=2), the first dropped.
        assert_eq!(b.widths, [0, 3]);
        assert_eq!(b.dropped, Some(0));
        assert_eq!(b.bits(), 3);
        // Nothing to drop when the grading is zero, and no terms give an empty box.
        let b = level_box(std::iter::empty(), &[0, 0], 4);
        assert_eq!((b.bits(), b.dropped), (0, None));
        let layout = LevelLayout::<u64>::new(&b, &[0, 0]).unwrap();
        let mut c = [7, 7];
        layout.decode(layout.hash(&[0, 0]), 0, &mut c);
        assert_eq!(c, [0, 0]);
    }

    #[test]
    fn test_level_table() {
        let mut rng = TestRng::new(8);
        let grading = [1, 1, 1, 1];
        let ts = terms(&mut rng, &grading, 10);
        let b = level_box(ts.iter().map(|(u, d)| (u.as_slice(), *d)), &grading, 30);
        let layout = LevelLayout::<u64>::new(&b, &grading).unwrap();
        let mut table: LevelTable<u64, u64> = LevelTable::new(layout.marker());
        assert_eq!(table.get(0), None);
        assert_eq!(table.capacity(), 0);
        for n in [10, 3000, 200, 0, 5000] {
            let mut model: HashMap<u64, u64> = HashMap::new();
            let mut order = Vec::new();
            for _ in 0..n {
                let mut x = vec![0i32; 4];
                for _ in 0..rng.range(0, 6) {
                    let (u, _) = &ts[rng.range(0, ts.len() as i64 - 1) as usize];
                    for (a, b) in x.iter_mut().zip(u) {
                        *a += b;
                    }
                }
                let h = layout.hash(&x);
                let v = rng.next_u64();
                let entry = table.get_or_insert(h);
                *entry = entry.wrapping_add(v);
                let m = model.entry(h).or_insert_with(|| {
                    order.push(h);
                    0
                });
                *m = m.wrapping_add(v);
            }
            assert_eq!(table.len(), model.len());
            assert!(4 * table.len() <= 3 * table.capacity());
            let entries: Vec<(u64, u64)> = table.iter().collect();
            let expected: Vec<(u64, u64)> = order.iter().map(|h| (*h, model[h])).collect();
            assert_eq!(entries, expected);
            for (h, v) in &model {
                assert_eq!(table.get(*h), Some(v));
            }
            table.clear();
            assert!(table.is_empty());
            assert_eq!(table.get(order.first().copied().unwrap_or_default()), None);
        }
        table.release();
        assert_eq!(table.capacity(), 0);
        *table.get_or_insert(0) = 3;
        assert_eq!(table.iter().collect::<Vec<_>>(), [(0, 3)]);
    }
}
