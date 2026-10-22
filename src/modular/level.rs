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

use super::key::{Key, Wide};
use super::table::prefetch;
use core::cmp::Reverse;
use core::fmt::Debug;
use core::ops::Range;

/// A word that offsets are packed into: `u64`, or `u128` when the box needs more bits, or for
/// keys wider than 128 bits, a [`Wide`] word of theirs.
pub trait LevelWord: Copy + Ord + Default + Debug + Send + Sync + 'static {
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

/// $a b \bmod 2^{64 N}$.
const fn wide_mul<const N: usize>(a: &[u64; N], b: &[u64; N]) -> [u64; N] {
    let mut r = [0u64; N];
    let mut i = 0;
    while i < N {
        let mut carry = 0u128;
        let mut j = 0;
        while i + j < N {
            let t = r[i + j] as u128 + a[i] as u128 * b[j] as u128 + carry;
            r[i + j] = t as u64;
            carry = t >> 64;
            j += 1;
        }
        i += 1;
    }
    r
}

/// $a - b \bmod 2^{64 N}$.
const fn wide_sub<const N: usize>(a: &[u64; N], b: &[u64; N]) -> [u64; N] {
    let mut r = [0u64; N];
    let mut borrow = 0u64;
    let mut i = 0;
    while i < N {
        let (d, b1) = a[i].overflowing_sub(b[i]);
        let (d, b2) = d.overflowing_sub(borrow);
        r[i] = d;
        borrow = (b1 || b2) as u64;
        i += 1;
    }
    r
}

/// An odd multiplier of $64 N$ bits.
const fn wide_multiplier<const N: usize>() -> [u64; N] {
    const WORDS: [u64; 4] = [
        0x9E37_79B9_7F4A_7C15,
        0xD6E8_FEB8_6659_FD93,
        0x94D0_49BB_1331_11EB,
        0xBF58_476D_1CE4_E5B9,
    ];
    let mut r = [0u64; N];
    let mut i = 0;
    while i < N {
        r[i] = WORDS[i % 4].rotate_left(8 * (i / 4) as u32);
        i += 1;
    }
    r
}

/// The inverse of the odd `m` modulo $2^{64 N}$, by the same Newton iteration.
const fn wide_inverse<const N: usize>(m: &[u64; N]) -> [u64; N] {
    let mut two = [0u64; N];
    two[0] = 2;
    let mut x = *m;
    let mut i = 0;
    while i < 10 {
        x = wide_mul(&x, &wide_sub(&two, &wide_mul(m, &x)));
        i += 1;
    }
    x
}

impl<const N: usize> LevelWord for Wide<N> {
    const BITS: u32 = 64 * N as u32;
    const MUL: Self = Wide(wide_multiplier::<N>());
    const MUL_INV: Self = Wide(wide_inverse::<N>(&wide_multiplier::<N>()));

    #[inline]
    fn wadd(self, other: Self) -> Self {
        super::key::Word::wrapping_add(self, other)
    }

    #[inline]
    fn wsub(self, other: Self) -> Self {
        super::key::Word::wrapping_sub(self, other)
    }

    #[inline]
    fn wmul(self, other: Self) -> Self {
        Wide(wide_mul(&self.0, &other.0))
    }

    #[inline]
    fn xor(self, other: Self) -> Self {
        let mut r = self.0;
        for (x, y) in r.iter_mut().zip(other.0) {
            *x ^= y;
        }
        Wide(r)
    }

    #[inline]
    fn from_i64(x: i64) -> Self {
        super::key::Word::from_i64(x)
    }

    #[inline]
    fn pow2(s: u32) -> Self {
        let mut r = [0u64; N];
        r[(s / 64) as usize] = 1 << (s % 64);
        Wide(r)
    }

    #[inline]
    fn field(self, offset: u32, width: u32) -> u64 {
        let (i, b) = ((offset / 64) as usize, offset % 64);
        let mut x = self.0[i] >> b;
        if b > 0 && i + 1 < N {
            x |= self.0[i + 1] << (64 - b);
        }
        if width >= 64 {
            x
        } else {
            x & ((1 << width) - 1)
        }
    }

    #[inline]
    fn top(self, bits: u32) -> usize {
        (self.0[N - 1] >> (64 - bits)) as usize
    }
}

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
    let mut bounds = BoxBounds::new(grading.len(), max_deg);
    for (u, deg) in terms {
        bounds.add(u, deg);
    }
    bounds.finish(grading)
}

/// [`level_box`] one term at a time, for terms that are not stored as coordinates.
pub struct BoxBounds {
    lo: Vec<i64>,
    hi: Vec<i64>,
    max_deg: i64,
}

impl BoxBounds {
    pub fn new(dim: usize, max_deg: u32) -> Self {
        Self {
            lo: vec![0; dim],
            hi: vec![0; dim],
            max_deg: max_deg as i64,
        }
    }

    /// Adds the term of coordinates `u` and degree `deg`, which must be positive.
    pub fn add(&mut self, u: &[i32], deg: u32) {
        assert!(deg > 0, "the terms of a level box need a positive degree");
        let deg = deg as i64;
        for (t, &x) in u.iter().enumerate() {
            let x = self.max_deg * x as i64;
            self.lo[t] = self.lo[t].min(x.div_euclid(deg));
            self.hi[t] = self.hi[t].max(-((-x).div_euclid(deg)));
        }
    }

    /// The box, with the widest coordinate of nonzero weight under `grading` left out.
    pub fn finish(self, grading: &[i32]) -> LevelBox {
        let Self { lo, hi, .. } = self;
        let mut widths: Vec<u32> = lo
            .iter()
            .zip(&hi)
            .map(|(&l, &h)| u64::BITS - ((h - l) as u64).leading_zeros())
            .collect();
        let dropped = (0..grading.len())
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
    ///
    /// The lookup is inlined into the scatter, and the insertion, much rarer, is not.
    #[inline(always)]
    pub fn get_or_insert(&mut self, h: W) -> &mut V {
        debug_assert!(h != self.marker, "the marker was inserted");
        if self.entries.is_empty() {
            return self.insert(h, None);
        }
        let x = h.xor(self.marker);
        let mask = self.entries.len() - 1;
        let mut i = h.top(self.log2);
        loop {
            let k = self.entries[i].0;
            if k == x {
                return &mut self.entries[i].1;
            }
            if k == W::default() {
                return self.insert(h, Some(i));
            }
            i = (i + 1) & mask;
        }
    }

    /// Inserts hash `h` as zero at the empty entry `at` where its probe ended, or allocates
    /// the table first if there is none, growing it if it is too full.
    #[cold]
    #[inline(never)]
    fn insert(&mut self, h: W, at: Option<usize>) -> &mut V {
        let at = match at {
            None => {
                self.entries = vec![(W::default(), V::default()); LEVEL_MIN_CAPACITY];
                self.log2 = LEVEL_MIN_CAPACITY.trailing_zeros();
                self.probe(h).unwrap_err()
            }
            Some(_) if 4 * (self.used.len() + 1) > 3 * self.entries.len() => {
                self.grow();
                self.probe(h).unwrap_err()
            }
            Some(i) => i,
        };
        self.entries[at] = (h.xor(self.marker), V::default());
        self.used.push(at as u32);
        &mut self.entries[at].1
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
    /// cache, ahead of a lookup a few iterations later. It does nothing before the table has
    /// memory.
    #[inline]
    pub fn prefetch(&self, h: W) {
        if !self.entries.is_empty() {
            prefetch(self.entries[h.top(self.log2)..].as_ptr());
        }
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

/// A compact open-addressing table of terms keyed by the hashes of their classes, for terms
/// that wait for their degree (the nodes of the extraction's tree): any capacity, the bucket
/// from the top 32 bits of the hash scaled to it, filled up to [`COMPACT_LOAD`] before it
/// grows by [`COMPACT_GROWTH`], and no list of the entries in use, so that it holds about
/// 1.4 entries per term rather than [`LevelTable`]'s two to three.
#[derive(Debug, Clone)]
pub struct CompactTable<W, V> {
    entries: Vec<(W, V)>,
    len: usize,
    marker: W,
}

/// The largest fraction of a [`CompactTable`] in use, as a numerator over 16.
pub const COMPACT_LOAD: usize = 13;

/// The factor a [`CompactTable`] grows by, as a numerator over 16.
pub const COMPACT_GROWTH: usize = 20;

/// The capacity of a [`CompactTable`] when it is first used.
const COMPACT_MIN_CAPACITY: usize = 16;

impl<W: LevelWord, V: Copy + Default> CompactTable<W, V> {
    pub fn new(marker: W) -> Self {
        Self {
            entries: Vec::new(),
            len: 0,
            marker,
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.len
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    #[inline]
    pub fn capacity(&self) -> usize {
        self.entries.len()
    }

    #[inline]
    fn bucket(h: W, capacity: usize) -> usize {
        ((h.top(32) as u64 * capacity as u64) >> 32) as usize
    }

    /// The value of hash `h`, if it has an entry.
    #[inline]
    pub fn get_mut(&mut self, h: W) -> Option<&mut V> {
        let n = self.entries.len();
        if n == 0 {
            return None;
        }
        let x = h.xor(self.marker);
        let mut i = Self::bucket(h, n);
        loop {
            let k = self.entries[i].0;
            if k == x {
                return Some(&mut self.entries[i].1);
            }
            if k == W::default() {
                return None;
            }
            i += 1;
            if i == n {
                i = 0;
            }
        }
    }

    /// The value of hash `h`, inserted as zero if it has no entry.
    #[inline(always)]
    pub fn get_or_insert(&mut self, h: W) -> &mut V {
        debug_assert!(h != self.marker, "the marker was inserted");
        let n = self.entries.len();
        if n == 0 {
            return self.insert(h);
        }
        let x = h.xor(self.marker);
        let mut i = Self::bucket(h, n);
        loop {
            let k = self.entries[i].0;
            if k == x {
                return &mut self.entries[i].1;
            }
            if k == W::default() {
                break;
            }
            i += 1;
            if i == n {
                i = 0;
            }
        }
        if 16 * (self.len + 1) > COMPACT_LOAD * n {
            return self.insert(h);
        }
        self.len += 1;
        self.entries[i] = (x, V::default());
        &mut self.entries[i].1
    }

    /// Inserts hash `h` as zero after allocating or growing the table.
    #[cold]
    #[inline(never)]
    fn insert(&mut self, h: W) -> &mut V {
        let n = self.entries.len();
        let capacity = if n == 0 {
            COMPACT_MIN_CAPACITY
        } else {
            (n * COMPACT_GROWTH / 16).max(n + 1)
        };
        assert!(capacity <= 1 << 32, "too many entries for a compact table");
        let old = core::mem::replace(
            &mut self.entries,
            vec![(W::default(), V::default()); capacity],
        );
        for entry in old {
            if entry.0 != W::default() {
                let mut j = Self::bucket(entry.0.xor(self.marker), capacity);
                while self.entries[j].0 != W::default() {
                    j += 1;
                    if j == capacity {
                        j = 0;
                    }
                }
                self.entries[j] = entry;
            }
        }
        let x = h.xor(self.marker);
        let mut i = Self::bucket(h, capacity);
        while self.entries[i].0 != W::default() {
            i += 1;
            if i == capacity {
                i = 0;
            }
        }
        self.len += 1;
        self.entries[i] = (x, V::default());
        &mut self.entries[i].1
    }

    /// Prefetches the entry where hash `h` would be probed first, as
    /// [`LevelTable::prefetch`] does.
    #[inline]
    pub fn prefetch(&self, h: W) {
        if !self.entries.is_empty() {
            prefetch(self.entries[Self::bucket(h, self.entries.len())..].as_ptr());
        }
    }

    /// The entries in use, as hashes and values.
    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = (W, V)> + '_ {
        let marker = self.marker;
        self.entries
            .iter()
            .filter(|e| e.0 != W::default())
            .map(move |&(k, v)| (k.xor(marker), v))
    }

    /// The terms of the entries in `range`, a range of `0..capacity()`, as hashes and values.
    pub fn iter_range(&self, range: Range<usize>) -> impl Iterator<Item = (W, V)> + '_ {
        let marker = self.marker;
        self.entries[range]
            .iter()
            .filter(|e| e.0 != W::default())
            .map(move |&(k, v)| (k.xor(marker), v))
    }

    /// Removes all entries and frees the memory.
    pub fn release(&mut self) {
        self.entries = Vec::new();
        self.len = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modular::key::{KeyLayout, PackedKey, PackedLayout, Word};
    use crate::modular::TestRng;
    use std::collections::HashMap;

    #[test]
    fn test_inverses() {
        assert_eq!(u64::MUL.wmul(u64::MUL_INV), 1);
        assert_eq!(u128::MUL.wmul(u128::MUL_INV), 1);
        let one = <Wide<4> as Word>::from_i64(1);
        assert_eq!(Wide::<4>::MUL.wmul(Wide::<4>::MUL_INV), one);
        let one = <Wide<8> as Word>::from_i64(1);
        assert_eq!(Wide::<8>::MUL.wmul(Wide::<8>::MUL_INV), one);
        // The wide product agrees with u128's.
        let mut rng = TestRng::new(2);
        for _ in 0..1000 {
            let (a, b) = (rng.next_u64(), rng.next_u64());
            let (c, d) = (rng.next_u64(), rng.next_u64());
            let x = (a as u128) << 64 | b as u128;
            let y = (c as u128) << 64 | d as u128;
            let p = Wide([b, a]).wmul(Wide([d, c]));
            assert_eq!(
                p,
                Wide([x.wrapping_mul(y) as u64, (x.wrapping_mul(y) >> 64) as u64])
            );
            for (o, w) in [(0, 64), (5, 64), (64, 3), (100, 28), (60, 10)] {
                assert_eq!(Wide([b, a]).field(o, w), x.field(o, w));
            }
        }
        for m in [1i64, 3, -5, 12345677, i64::MAX, -1] {
            assert_eq!(m.wrapping_mul(inverse_mod_2_64(m)), 1);
        }
    }

    /// Random terms of positive degree under `grading`, in the quadrant-like region where
    /// the grading is positive.
    /// With more than ten coordinates, the terms have about three nonzero ones each, which
    /// makes their degrees small and the box wide.
    fn terms(rng: &mut TestRng, grading: &[i32], n: usize) -> Vec<(Vec<i32>, u32)> {
        let mut res = Vec::new();
        let dim = grading.len() as i64;
        while res.len() < n {
            let u: Vec<i32> = grading
                .iter()
                .map(|_| {
                    if dim > 10 && rng.range(1, dim) > 3 {
                        0
                    } else {
                        rng.range(-3, 6) as i32
                    }
                })
                .collect();
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

    /// Hashes add like the offsets, and decode to them, for sums of terms up to the degree,
    /// when the box fits `W`; the bits of the box.
    fn check_layout<W: LevelWord>(grading: &[i32], seed: u64) -> u32 {
        let mut rng = TestRng::new(seed);
        let ts = terms(&mut rng, grading, 12.max(3 * grading.len()));
        let max_deg = 40;
        let b = level_box(ts.iter().map(|(u, d)| (u.as_slice(), *d)), grading, max_deg);
        let Some(layout) = LevelLayout::<W>::new(&b, grading) else {
            assert!(b.bits() >= W::BITS);
            return b.bits();
        };
        let dim = grading.len();
        let mut coords = vec![0; dim];
        // Keys of a packed layout wide enough for the box: a coordinate is at most the degree
        // times six.
        let packed = PackedLayout::<Wide<8>>::new(&vec![1 << 8; dim]).unwrap();
        let units: Vec<PackedKey<Wide<8>>> = (0..dim)
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
        b.bits()
    }

    #[test]
    fn test_layout() {
        for (i, grading) in [
            vec![1, 1],
            vec![1, 2, 3],
            vec![6, -12, -14, -17, -18, -4, -5, -12],
            vec![0, 1, 0, 2, 5],
            vec![-1, 3],
            // Boxes of more than 128 and 256 bits.
            vec![1; 24],
            vec![1; 48],
        ]
        .iter()
        .enumerate()
        {
            check_layout::<u64>(grading, i as u64);
            check_layout::<u128>(grading, i as u64);
            check_layout::<Wide<4>>(grading, i as u64);
            let bits = check_layout::<Wide<8>>(grading, i as u64);
            assert_eq!(bits > 128, i >= 5, "{bits}");
            assert_eq!(bits > 256, i >= 6, "{bits}");
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
