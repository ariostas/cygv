//! Tables keyed by curve class.
//!
//! - [`SlotMap`] maps keys to slots, each holding the key, its degree and a fixed number of
//!   values. It is where series are accumulated, and the slot numbers stay valid while it
//!   grows, so a [`DegreeHeap`] can refer to them.
//! - [`SparsePoly`] is a finished series: its nonzero terms in order of degree, then key.
//! - [`DegreeHeap`] hands out slots in order of degree, for the loops that finalize the terms
//!   of a series degree by degree while adding new ones above.
//!
//! A slot or term carries `width` values rather than one, fixed when the table is made. A
//! threefold needs one value per term, but an n-fold needs one per reference surface, a
//! number only known at run time.

use super::key::Key;
use core::cmp::Reverse;
use rayon::prelude::*;
use std::collections::BinaryHeap;
use std::sync::Mutex;

/// Marks an unused bucket of a [`SlotMap`].
const EMPTY: u64 = 0;

/// Whether all values are zero, i.e. equal to [`Default`].
#[inline]
fn all_zero<V: Default + PartialEq>(vals: &[V]) -> bool {
    vals.iter().all(|v| *v == V::default())
}

/// A hash map from curve classes to slots, each with the class, its degree and `width` values.
///
/// Slots are numbered from zero in order of insertion. Lookups go through an open-addressing
/// index of 8-byte buckets, each holding a slot number and a 32-bit tag from the high bits
/// of the key's hash (the low ones pick the bucket), so that a probe only follows the slot
/// number to compare keys when the tags match. The index is at most half full.
///
/// Slot numbers are `u32`, which limits a table to about four billion slots. That is far
/// more than fits in memory in any case: each slot takes at least 36 bytes with packed keys.
#[derive(Debug, Clone)]
pub struct SlotMap<K, V> {
    width: usize,
    /// The index: per bucket the tag and one more than the slot, or [`EMPTY`].
    index: Vec<u64>,
    keys: Vec<K>,
    degs: Vec<u32>,
    /// `width` values per slot.
    vals: Vec<V>,
    /// Scratch space for [`SlotMap::clear`].
    buckets: Vec<usize>,
}

impl<K: Key, V: Copy + Default + PartialEq> SlotMap<K, V> {
    /// An empty map with `width` values per slot and room for `capacity` slots before it
    /// grows.
    ///
    /// # Panics
    ///
    /// Panics if `width` is zero.
    pub fn new(width: usize, capacity: usize) -> Self {
        assert!(width > 0, "a slot map needs at least one value per slot");
        let n_buckets = (2 * capacity + 1).max(16).next_power_of_two();
        Self {
            width,
            index: vec![EMPTY; n_buckets],
            keys: Vec::with_capacity(capacity),
            degs: Vec::with_capacity(capacity),
            vals: Vec::with_capacity(capacity * width),
            buckets: Vec::new(),
        }
    }

    /// The number of values per slot.
    #[inline]
    pub fn width(&self) -> usize {
        self.width
    }

    /// The number of slots.
    #[inline]
    pub fn len(&self) -> usize {
        self.keys.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    /// The bucket that holds `key`, whose hash is `hash`, with its slot, or else the empty
    /// bucket where it would go.
    #[inline]
    fn probe(&self, key: K, hash: u64) -> (usize, Option<u32>) {
        let mask = self.index.len() - 1;
        let tag = hash >> 32;
        let mut h = hash as usize & mask;
        loop {
            let e = self.index[h];
            if e == EMPTY {
                return (h, None);
            }
            if e >> 32 == tag {
                let s = (e as u32).wrapping_sub(1);
                if self.keys[s as usize] == key {
                    return (h, Some(s));
                }
            }
            h = (h + 1) & mask;
        }
    }

    /// The slot of `key`, if it has one.
    #[inline]
    pub fn get(&self, key: K) -> Option<u32> {
        self.probe(key, key.hash64()).1
    }

    /// The slot of `key`, created with zero values if there is none, and whether it was
    /// created.
    ///
    /// The degree is only stored when the slot is created.
    #[inline]
    pub fn insert(&mut self, key: K, deg: u32) -> (u32, bool) {
        self.insert_hashed(key, key.hash64(), deg)
    }

    /// [`SlotMap::insert`] with the hash of the key, [`Key::hash64`], already computed.
    #[inline]
    pub fn insert_hashed(&mut self, key: K, hash: u64, deg: u32) -> (u32, bool) {
        let (h, found) = self.probe(key, hash);
        if let Some(s) = found {
            return (s, false);
        }
        let s = self.keys.len();
        assert!(s < u32::MAX as usize, "too many slots for a slot map");
        self.keys.push(key);
        self.degs.push(deg);
        self.vals.resize(self.vals.len() + self.width, V::default());
        self.index[h] = (hash >> 32) << 32 | (s as u64 + 1);
        if 2 * self.keys.len() >= self.index.len() {
            self.rehash(2 * self.index.len());
        }
        (s as u32, true)
    }

    /// Asks the processor to fetch the bucket where a key with hash `hash` is probed first
    /// into the cache, ahead of an insertion a few iterations later. It does nothing on other
    /// architectures than x86-64.
    #[inline]
    pub fn prefetch(&self, hash: u64) {
        #[cfg(target_arch = "x86_64")]
        {
            let h = hash as usize & (self.index.len() - 1);
            let ptr = self.index[h..].as_ptr().cast::<i8>();
            // SAFETY: a prefetch has no effect on memory, and SSE is part of x86-64.
            #[allow(unused_unsafe)]
            unsafe {
                core::arch::x86_64::_mm_prefetch(ptr, core::arch::x86_64::_MM_HINT_T0);
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        let _ = hash;
    }

    /// The values of `key`'s slot, created with zero values if there is none.
    #[inline]
    pub fn get_or_insert(&mut self, key: K, deg: u32) -> &mut [V] {
        let (s, _) = self.insert(key, deg);
        self.values_mut(s)
    }

    /// Rebuilds the index with `n_buckets` buckets.
    fn rehash(&mut self, n_buckets: usize) {
        self.index = vec![EMPTY; n_buckets];
        for (s, &key) in self.keys.iter().enumerate() {
            let hash = key.hash64();
            let mut h = hash as usize & (n_buckets - 1);
            while self.index[h] != EMPTY {
                h = (h + 1) & (n_buckets - 1);
            }
            self.index[h] = (hash >> 32) << 32 | (s as u64 + 1);
        }
    }

    #[inline]
    pub fn key(&self, slot: u32) -> K {
        self.keys[slot as usize]
    }

    #[inline]
    pub fn deg(&self, slot: u32) -> u32 {
        self.degs[slot as usize]
    }

    #[inline]
    pub fn values(&self, slot: u32) -> &[V] {
        let i = slot as usize * self.width;
        &self.vals[i..i + self.width]
    }

    #[inline]
    pub fn values_mut(&mut self, slot: u32) -> &mut [V] {
        let i = slot as usize * self.width;
        &mut self.vals[i..i + self.width]
    }

    /// The keys of all slots, in slot order.
    #[inline]
    pub fn keys(&self) -> &[K] {
        &self.keys
    }

    /// The degrees of all slots, in slot order.
    #[inline]
    pub fn degs(&self) -> &[u32] {
        &self.degs
    }

    /// Removes all slots, keeping the memory.
    ///
    /// When the map holds few slots compared with the size of its index, only their buckets
    /// are reset, so that clearing a large but mostly empty map stays cheap.
    pub fn clear(&mut self) {
        if 8 * self.keys.len() >= self.index.len() {
            self.index.fill(EMPTY);
        } else {
            // Find every bucket before emptying any, so that the probe sequences stay intact.
            let mut buckets = core::mem::take(&mut self.buckets);
            buckets.clear();
            buckets.extend(self.keys.iter().map(|&key| self.probe(key, key.hash64()).0));
            for &h in &buckets {
                self.index[h] = EMPTY;
            }
            self.buckets = buckets;
        }
        self.keys.clear();
        self.degs.clear();
        self.vals.clear();
    }

    /// Takes the slots with a nonzero value as a [`SparsePoly`], and clears the map.
    pub fn take(&mut self) -> SparsePoly<K, V> {
        let mut order: Vec<u32> = (0..self.len() as u32)
            .filter(|&s| !all_zero(self.values(s)))
            .collect();
        order.sort_unstable_by_key(|&s| (self.deg(s), self.key(s)));
        let poly = SparsePoly {
            width: self.width,
            keys: order.iter().map(|&s| self.key(s)).collect(),
            degs: order.iter().map(|&s| self.deg(s)).collect(),
            vals: order
                .iter()
                .flat_map(|&s| self.values(s).iter().copied())
                .collect(),
        };
        self.clear();
        poly
    }
}

/// A sparse series: its nonzero terms in increasing order of degree, then key, each with a
/// key, a degree and `width` values.
#[derive(Debug, Clone, PartialEq)]
pub struct SparsePoly<K, V> {
    width: usize,
    keys: Vec<K>,
    degs: Vec<u32>,
    vals: Vec<V>,
}

impl<K: Key, V: Copy + Default + PartialEq> SparsePoly<K, V> {
    /// The series without terms.
    pub fn new(width: usize) -> Self {
        Self {
            width,
            keys: Vec::new(),
            degs: Vec::new(),
            vals: Vec::new(),
        }
    }

    /// The number of values per term.
    #[inline]
    pub fn width(&self) -> usize {
        self.width
    }

    /// The number of terms.
    #[inline]
    pub fn len(&self) -> usize {
        self.keys.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    #[inline]
    pub fn keys(&self) -> &[K] {
        &self.keys
    }

    #[inline]
    pub fn degs(&self) -> &[u32] {
        &self.degs
    }

    /// The values of the `i`-th term.
    #[inline]
    pub fn values(&self, i: usize) -> &[V] {
        &self.vals[i * self.width..(i + 1) * self.width]
    }

    /// Appends a term.
    ///
    /// The terms must be pushed in increasing order of degree, then key, and not all of
    /// `values` may be zero; this is only checked in debug builds.
    ///
    /// # Panics
    ///
    /// Panics if `values` does not have `width` entries.
    pub fn push(&mut self, key: K, deg: u32, values: &[V]) {
        assert_eq!(
            values.len(),
            self.width,
            "a term needs one value per component"
        );
        debug_assert!(
            !all_zero(values),
            "the terms of a sparse series are nonzero"
        );
        debug_assert!(
            self.keys
                .last()
                .is_none_or(|&k| (*self.degs.last().unwrap(), k) < (deg, key)),
            "the terms of a sparse series are pushed in order"
        );
        self.keys.push(key);
        self.degs.push(deg);
        self.vals.extend_from_slice(values);
    }

    /// Appends the terms of `other`, which must all come after those of `self` in the order
    /// of degree, then key; this is only checked in debug builds.
    ///
    /// # Panics
    ///
    /// Panics if the widths differ.
    pub fn append(&mut self, mut other: Self) {
        assert_eq!(
            self.width, other.width,
            "appended series need the same width"
        );
        debug_assert!(
            self.keys
                .last()
                .zip(other.keys.first())
                .is_none_or(|(&a, &b)| { (*self.degs.last().unwrap(), a) < (other.degs[0], b) }),
            "the terms of a sparse series are appended in order"
        );
        self.keys.append(&mut other.keys);
        self.degs.append(&mut other.degs);
        self.vals.append(&mut other.vals);
    }

    /// The values of the `i`-th term, to change in place. A term must not be left with only
    /// zeros but by [`SparsePoly::retain_terms`] dropping it.
    #[inline]
    pub fn values_mut(&mut self, i: usize) -> &mut [V] {
        &mut self.vals[i * self.width..(i + 1) * self.width]
    }

    /// Keeps only the terms whose key satisfies `keep`, in their order.
    pub fn retain(&mut self, mut keep: impl FnMut(K) -> bool) {
        self.retain_terms(|key, _, _| keep(key));
    }

    /// Keeps only the terms whose key, degree and values satisfy `keep`, in their order.
    pub fn retain_terms(&mut self, mut keep: impl FnMut(K, u32, &[V]) -> bool) {
        let mut n = 0;
        for i in 0..self.len() {
            if !keep(self.keys[i], self.degs[i], self.values(i)) {
                continue;
            }
            self.keys[n] = self.keys[i];
            self.degs[n] = self.degs[i];
            self.vals
                .copy_within(i * self.width..(i + 1) * self.width, n * self.width);
            n += 1;
        }
        self.keys.truncate(n);
        self.degs.truncate(n);
        self.vals.truncate(n * self.width);
    }

    /// Adds the terms of `other`, which has the same width and no key in common with `self`,
    /// merging them into their places from the end, without a second buffer.
    ///
    /// # Panics
    ///
    /// Panics if the widths differ.
    pub fn merge_from(&mut self, other: Self) {
        assert_eq!(self.width, other.width, "merged series need the same width");
        let (w, n, m) = (self.width, self.len(), other.len());
        self.keys.reserve_exact(m);
        self.degs.reserve_exact(m);
        self.vals.reserve_exact(m * w);
        self.keys.extend_from_slice(&other.keys);
        self.degs.extend_from_slice(&other.degs);
        self.vals.extend_from_slice(&other.vals);
        let (mut i, mut j) = (n, m);
        while j > 0 {
            let k = i + j - 1;
            if i > 0
                && (self.degs[i - 1], self.keys[i - 1]) > (other.degs[j - 1], other.keys[j - 1])
            {
                i -= 1;
                self.keys[k] = self.keys[i];
                self.degs[k] = self.degs[i];
                self.vals.copy_within(i * w..(i + 1) * w, k * w);
            } else {
                j -= 1;
                self.keys[k] = other.keys[j];
                self.degs[k] = other.degs[j];
                self.vals[k * w..(k + 1) * w].copy_from_slice(other.values(j));
            }
        }
    }

    /// The number of terms of degree at most `deg`, i.e. the end of their range.
    #[inline]
    pub fn len_upto(&self, deg: u32) -> usize {
        self.degs.partition_point(|&d| d <= deg)
    }

    /// The series with the terms of all of `parts`, which must have the same width and no key
    /// in common, such as the shards of a table split by key.
    ///
    /// # Panics
    ///
    /// Panics if the widths differ.
    pub fn merge(width: usize, parts: Vec<Self>) -> Self {
        assert!(
            parts.iter().all(|p| p.width == width),
            "merged series need the same width"
        );
        let mut parts: Vec<Self> = parts.into_iter().filter(|p| !p.is_empty()).collect();
        if parts.len() <= 1 {
            return parts.pop().unwrap_or_else(|| Self::new(width));
        }
        // A merge of the parts if each is in order, which they usually are; otherwise a sort.
        let sorted = parts
            .iter()
            .all(|p| (1..p.len()).all(|i| (p.degs[i - 1], p.keys[i - 1]) < (p.degs[i], p.keys[i])));
        let n: usize = parts.iter().map(Self::len).sum();
        let mut res = Self::new(width);
        res.keys.reserve_exact(n);
        res.degs.reserve_exact(n);
        res.vals.reserve_exact(n * width);
        if sorted {
            let mut heap: BinaryHeap<Reverse<(u32, K, u32)>> = parts
                .iter()
                .enumerate()
                .map(|(p, part)| Reverse((part.degs[0], part.keys[0], p as u32)))
                .collect();
            let mut at = vec![0usize; parts.len()];
            while let Some(Reverse((deg, key, p))) = heap.pop() {
                let (p, part) = (p as usize, &parts[p as usize]);
                let i = at[p];
                debug_assert!(
                    res.keys.last().is_none_or(|&k| k != key),
                    "merged series have no key in common"
                );
                res.keys.push(key);
                res.degs.push(deg);
                res.vals.extend_from_slice(part.values(i));
                at[p] += 1;
                if at[p] < part.len() {
                    heap.push(Reverse((part.degs[i + 1], part.keys[i + 1], p as u32)));
                }
            }
            return res;
        }
        let mut order: Vec<(u32, K, u32, u32)> = parts
            .iter()
            .enumerate()
            .flat_map(|(p, part)| {
                (0..part.len()).map(move |i| (part.degs[i], part.keys[i], p as u32, i as u32))
            })
            .collect();
        order.par_sort_unstable();
        debug_assert!(
            order.windows(2).all(|w| w[0].1 != w[1].1),
            "merged series have no key in common"
        );
        for (deg, key, p, i) in order {
            res.keys.push(key);
            res.degs.push(deg);
            res.vals
                .extend_from_slice(parts[p as usize].values(i as usize));
        }
        res
    }
}

/// The shard of a key among `n_shards`, from the high bits of its hash, since tables index by
/// the low ones; at most $2^{32}$ shards. The tags of a [`SlotMap`] overlap these bits, so with a power of two of
/// shards a few bits of the tags are the same throughout a shard, which only makes their
/// collisions a little more likely.
#[inline]
pub fn shard_of<K: Key>(key: K, n_shards: usize) -> usize {
    // The bits from 32 on, scaled to the shards by a multiplication rather than a division.
    (((key.hash64() >> 32) * n_shards as u64) >> 32) as usize
}

/// A [`SlotMap`] split into shards by key, each behind its own lock, which several threads
/// add to at once.
#[derive(Debug)]
pub struct ShardedMap<K, V> {
    shards: Vec<Mutex<SlotMap<K, V>>>,
}

impl<K: Key, V: Copy + Default + PartialEq + Send> ShardedMap<K, V> {
    /// An empty map with `width` values per slot, in `n_shards` shards with room for
    /// `capacity` slots in all.
    pub fn new(width: usize, n_shards: usize, capacity: usize) -> Self {
        let n_shards = n_shards.max(1);
        Self {
            shards: (0..n_shards)
                .map(|_| Mutex::new(SlotMap::new(width, capacity / n_shards)))
                .collect(),
        }
    }

    /// Adds the values of every slot of `map` to those of the same key, with `add`, and
    /// clears `map`. Each shard is locked once.
    pub fn add_from(&self, map: &mut SlotMap<K, V>, add: impl Fn(&mut [V], &[V])) {
        let n = self.shards.len();
        let mut by_shard: Vec<Vec<u32>> = vec![Vec::new(); n];
        for (s, &key) in map.keys().iter().enumerate() {
            by_shard[shard_of(key, n)].push(s as u32);
        }
        for (shard, slots) in self.shards.iter().zip(by_shard) {
            if slots.is_empty() {
                continue;
            }
            let mut shard = shard.lock().unwrap();
            for s in slots {
                add(shard.get_or_insert(map.key(s), map.deg(s)), map.values(s));
            }
        }
        map.clear();
    }

    /// Adds the values of every slot of each of `maps` to those of the same key, with `add`,
    /// the maps one after the other, and clears them. The shards take their slots in
    /// parallel, so that a sum is taken in the order of the maps whatever the threads.
    pub fn add_in_order(&self, maps: &mut [SlotMap<K, V>], add: impl Fn(&mut [V], &[V]) + Sync)
    where
        K: Sync,
        V: Sync,
    {
        let n = self.shards.len();
        let by_shard: Vec<Vec<Vec<u32>>> = maps
            .par_iter()
            .map(|map| {
                let mut by_shard = vec![Vec::new(); n];
                for (s, &key) in map.keys().iter().enumerate() {
                    by_shard[shard_of(key, n)].push(s as u32);
                }
                by_shard
            })
            .collect();
        self.shards.par_iter().enumerate().for_each(|(i, shard)| {
            let mut shard = shard.lock().unwrap();
            for (map, slots) in maps.iter().zip(&by_shard) {
                for &s in &slots[i] {
                    add(shard.get_or_insert(map.key(s), map.deg(s)), map.values(s));
                }
            }
        });
        maps.iter_mut().for_each(SlotMap::clear);
    }

    /// The slots with a nonzero value, as a [`SparsePoly`].
    pub fn take(self) -> SparsePoly<K, V> {
        let parts: Vec<SparsePoly<K, V>> = self
            .shards
            .into_par_iter()
            .map(|s| s.into_inner().unwrap().take())
            .collect();
        let width = parts[0].width();
        SparsePoly::merge(width, parts)
    }
}

/// A min-heap of slots ordered by degree. Slots of equal degree come out in increasing order.
#[derive(Debug, Clone, Default)]
pub struct DegreeHeap {
    heap: BinaryHeap<Reverse<(u32, u32)>>,
}

impl DegreeHeap {
    pub fn new() -> Self {
        Self::default()
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.heap.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.heap.is_empty()
    }

    #[inline]
    pub fn push(&mut self, deg: u32, slot: u32) {
        self.heap.push(Reverse((deg, slot)));
    }

    /// The lowest degree of a slot in the heap.
    #[inline]
    pub fn min_deg(&self) -> Option<u32> {
        self.heap.peek().map(|Reverse((deg, _))| *deg)
    }

    /// Removes and returns a slot of the lowest degree.
    #[inline]
    pub fn pop(&mut self) -> Option<u32> {
        self.heap.pop().map(|Reverse((_, slot))| slot)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modular::key::{KeyLayout, PackedKey, PackedLayout};
    use crate::modular::TestRng;
    use std::collections::HashMap;

    /// A key from a pool of `n` keys spread over the whole range, so that keys repeat.
    fn pool_key(rng: &mut TestRng, n: u64) -> PackedKey {
        let i = rng.next_u64() % n;
        PackedKey((i as u128).wrapping_mul(0x9E37_79B9_7F4A_7C15_F39C_C060_5CED_C835))
    }

    #[test]
    fn test_slot_map() {
        let mut rng = TestRng::new(5);
        // Starting small, so that the map grows several times.
        let mut map: SlotMap<PackedKey, u64> = SlotMap::new(2, 4);
        let mut model: HashMap<PackedKey, (u32, u32, [u64; 2])> = HashMap::new();
        for round in 0..3 {
            for _ in 0..5000 {
                let key = pool_key(&mut rng, 3000);
                let deg = (key.0 % 31) as u32;
                let (slot, created) = map.insert(key, deg);
                assert_eq!(created, !model.contains_key(&key));
                let entry = model.entry(key).or_insert((slot, deg, [0; 2]));
                assert_eq!(entry.0, slot);
                let (a, b) = (rng.next_u64(), rng.next_u64());
                let vals = map.values_mut(slot);
                vals[0] = vals[0].wrapping_add(a);
                vals[1] ^= b;
                entry.2[0] = entry.2[0].wrapping_add(a);
                entry.2[1] ^= b;
            }
            assert_eq!(map.len(), model.len());
            for (&key, &(slot, deg, vals)) in &model {
                assert_eq!(map.get(key), Some(slot));
                assert_eq!(map.key(slot), key);
                assert_eq!(map.deg(slot), deg);
                assert_eq!(map.values(slot), vals);
            }
            for _ in 0..1000 {
                let key = pool_key(&mut rng, 1 << 40);
                assert_eq!(map.get(key), model.get(&key).map(|e| e.0));
            }
            // Clear in both modes: after the first round the map is full, and in the later
            // ones only a few slots go in before the clear.
            map.clear();
            model.clear();
            assert!(map.is_empty());
            for _ in 0..(10 * round) {
                let key = pool_key(&mut rng, 3000);
                map.get_or_insert(key, 1).copy_from_slice(&[7, 7]);
            }
            map.clear();
            for _ in 0..1000 {
                assert_eq!(map.get(pool_key(&mut rng, 3000)), None);
            }
        }
    }

    /// Adding maps to a sharded map from several threads, merging its shards, and appending
    /// series give what one map gives.
    #[test]
    fn test_sharded() {
        let mut rng = TestRng::new(8);
        let mut whole: SlotMap<PackedKey, u64> = SlotMap::new(2, 16);
        let mut parts: Vec<SlotMap<PackedKey, u64>> = Vec::new();
        for _ in 0..12 {
            let mut part = SlotMap::new(2, 16);
            for _ in 0..200 {
                let key = pool_key(&mut rng, 300);
                let deg = (key.hash64() % 7) as u32;
                let v = [rng.range(0, 5) as u64, rng.range(0, 5) as u64];
                for (map, v) in [(&mut whole, v), (&mut part, v)] {
                    for (x, y) in map.get_or_insert(key, deg).iter_mut().zip(v) {
                        *x = (*x + y) % 5;
                    }
                }
            }
            parts.push(part);
        }
        let expected = whole.take();
        let add = |out: &mut [u64], v: &[u64]| {
            for (x, &y) in out.iter_mut().zip(v) {
                *x = (*x + y) % 5;
            }
        };
        for n_shards in [1, 3, 16] {
            let sharded = ShardedMap::new(2, n_shards, 64);
            parts.clone().into_par_iter().for_each(|mut part| {
                sharded.add_from(&mut part, add);
                assert!(part.is_empty());
            });
            assert_eq!(sharded.take(), expected, "{n_shards} shards");
        }
        // Split by degree and appended in order.
        let mut appended = SparsePoly::new(2);
        for d in 0..7 {
            let mut part = SparsePoly::new(2);
            for i in 0..expected.len() {
                if expected.degs()[i] == d {
                    part.push(expected.keys()[i], d, expected.values(i));
                }
            }
            appended.append(part);
        }
        assert_eq!(appended, expected);
    }

    #[test]
    fn test_take() {
        let mut rng = TestRng::new(6);
        let layout = PackedLayout::<u128>::new(&[20, 20, 20]).unwrap();
        let mut map: SlotMap<PackedKey, u64> = SlotMap::new(2, 16);
        let mut model = Vec::new();
        for _ in 0..500 {
            let c: Vec<i32> = (0..3).map(|_| rng.range(-20, 20) as i32).collect();
            let key = layout.pack(&c);
            let deg = rng.range(0, 10) as u32;
            let vals = [rng.range(0, 2) as u64, rng.range(0, 2) as u64];
            let (slot, created) = map.insert(key, deg);
            if !created {
                continue;
            }
            map.values_mut(slot).copy_from_slice(&vals);
            if vals != [0, 0] {
                model.push((deg, key, vals));
            }
        }
        model.sort_unstable_by_key(|&(d, k, _)| (d, k));
        let poly = map.take();
        assert!(map.is_empty());
        assert_eq!(poly.width(), 2);
        assert_eq!(poly.len(), model.len());
        for (i, (d, k, v)) in model.iter().enumerate() {
            assert_eq!(
                (poly.degs()[i], poly.keys()[i], poly.values(i)),
                (*d, *k, &v[..])
            );
        }
        for deg in 0..=11 {
            let n = model.iter().filter(|e| e.0 <= deg).count();
            assert_eq!(poly.len_upto(deg), n);
        }
        // Pushing the same terms gives the same series.
        let mut pushed = SparsePoly::new(2);
        for (d, k, v) in &model {
            pushed.push(*k, *d, v);
        }
        assert_eq!(pushed, poly);
        // Retaining keeps the order and the values of the kept terms.
        let mut kept = SparsePoly::new(2);
        for (d, k, v) in model.iter().filter(|e| e.1 .0 % 3 != 0) {
            kept.push(*k, *d, v);
        }
        pushed.retain(|k| k.0 % 3 != 0);
        assert_eq!(pushed, kept);
        assert!(kept.len() < poly.len());
        // The map is usable again after taking.
        assert_eq!(map.insert(PackedKey(5), 0), (0, true));
        assert!(SparsePoly::<PackedKey, u64>::new(1).is_empty());
    }

    #[test]
    fn test_degree_heap() {
        let mut rng = TestRng::new(7);
        let mut heap = DegreeHeap::new();
        let mut degs = Vec::new();
        for slot in 0..1000 {
            let deg = rng.range(0, 50) as u32;
            heap.push(deg, slot);
            degs.push(deg);
        }
        assert_eq!(heap.len(), 1000);
        let mut last = None;
        while let Some(deg) = heap.min_deg() {
            let slot = heap.pop().unwrap();
            let cur = (degs[slot as usize], slot);
            assert_eq!(cur.0, deg);
            assert!(last < Some(cur));
            last = Some(cur);
        }
        assert!(heap.is_empty());
        assert_eq!(heap.pop(), None);
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "pushed in order")]
    fn test_push_out_of_order() {
        let mut poly = SparsePoly::new(1);
        poly.push(PackedKey(1u128), 1, &[1u64]);
        poly.push(PackedKey(0), 1, &[1u64]);
    }
}
