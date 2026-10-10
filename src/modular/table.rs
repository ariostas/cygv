//! Tables keyed by curve class.
//!
//! - [`SlotMap`] maps keys to slots, each holding the key, its degree and a fixed number of
//!   values. It is where series are accumulated, and the slot numbers stay valid while it
//!   grows, so a [`DegreeHeap`] can refer to them.
//! - [`SparsePoly`] is a finished series: its nonzero terms in order of degree, then key.
//! - [`DegreeHeap`] hands out slots in order of degree, for the loops that finalize the terms
//!   of a series degree by degree while adding new ones above.
//! - [`LevelTable`] holds one degree level of the terms of a single exponential, with the
//!   value stored next to the key, so that one probe touches both.
//!
//! A slot or term carries `width` values rather than one, fixed when the table is made. A
//! threefold needs one value per term, but an n-fold needs one per reference surface, a
//! number only known at run time.

use super::key::Key;
use core::cmp::Reverse;
use std::collections::BinaryHeap;

/// Marks an unused bucket of a [`SlotMap`].
const EMPTY: u32 = u32::MAX;

/// Whether all values are zero, i.e. equal to [`Default`].
#[inline]
fn all_zero<V: Default + PartialEq>(vals: &[V]) -> bool {
    vals.iter().all(|v| *v == V::default())
}

/// A hash map from curve classes to slots, each with the class, its degree and `width` values.
///
/// Slots are numbered from zero in order of insertion. Lookups go through an open-addressing
/// index that keeps a copy of each key next to its slot number, so that a probe does not
/// have to follow the slot number to compare keys. The index is at most half full.
///
/// Slot numbers are `u32`, which limits a table to about four billion slots. That is far
/// more than fits in memory in any case: each slot takes at least 36 bytes with packed keys.
#[derive(Debug, Clone)]
pub struct SlotMap<K, V> {
    width: usize,
    /// The index: per bucket a key and its slot, or [`EMPTY`].
    hkeys: Vec<K>,
    hslots: Vec<u32>,
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
            hkeys: vec![K::default(); n_buckets],
            hslots: vec![EMPTY; n_buckets],
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

    /// The bucket that holds `key`, with its slot, or else the empty bucket where it would go.
    #[inline]
    fn probe(&self, key: K) -> (usize, Option<u32>) {
        let mask = self.hslots.len() - 1;
        let mut h = key.hash64() as usize & mask;
        loop {
            let s = self.hslots[h];
            if s == EMPTY {
                return (h, None);
            }
            if self.hkeys[h] == key {
                return (h, Some(s));
            }
            h = (h + 1) & mask;
        }
    }

    /// The slot of `key`, if it has one.
    #[inline]
    pub fn get(&self, key: K) -> Option<u32> {
        self.probe(key).1
    }

    /// The slot of `key`, created with zero values if there is none, and whether it was
    /// created.
    ///
    /// The degree is only stored when the slot is created.
    #[inline]
    pub fn insert(&mut self, key: K, deg: u32) -> (u32, bool) {
        let (h, found) = self.probe(key);
        if let Some(s) = found {
            return (s, false);
        }
        let s = self.keys.len();
        assert!(s < EMPTY as usize, "too many slots for a slot map");
        self.keys.push(key);
        self.degs.push(deg);
        self.vals.resize(self.vals.len() + self.width, V::default());
        self.hkeys[h] = key;
        self.hslots[h] = s as u32;
        if 2 * self.keys.len() >= self.hslots.len() {
            self.rehash(2 * self.hslots.len());
        }
        (s as u32, true)
    }

    /// The values of `key`'s slot, created with zero values if there is none.
    #[inline]
    pub fn get_or_insert(&mut self, key: K, deg: u32) -> &mut [V] {
        let (s, _) = self.insert(key, deg);
        self.values_mut(s)
    }

    /// Rebuilds the index with `n_buckets` buckets.
    fn rehash(&mut self, n_buckets: usize) {
        self.hkeys = vec![K::default(); n_buckets];
        self.hslots = vec![EMPTY; n_buckets];
        for (s, &key) in self.keys.iter().enumerate() {
            let mut h = key.hash64() as usize & (n_buckets - 1);
            while self.hslots[h] != EMPTY {
                h = (h + 1) & (n_buckets - 1);
            }
            self.hkeys[h] = key;
            self.hslots[h] = s as u32;
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
        if 8 * self.keys.len() >= self.hslots.len() {
            self.hslots.fill(EMPTY);
        } else {
            // Find every bucket before emptying any, so that the probe sequences stay intact.
            let mut buckets = core::mem::take(&mut self.buckets);
            buckets.clear();
            buckets.extend(self.keys.iter().map(|&key| self.probe(key).0));
            for &h in &buckets {
                self.hslots[h] = EMPTY;
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

    /// The number of terms of degree at most `deg`, i.e. the end of their range.
    #[inline]
    pub fn len_upto(&self, deg: u32) -> usize {
        self.degs.partition_point(|&d| d <= deg)
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

/// The capacity of a [`LevelTable`] when it is first used.
const LEVEL_MIN_CAPACITY: usize = 256;

/// An open-addressing table of (key, value) entries, for the terms of one degree of an
/// exponential being built.
///
/// Unlike a [`SlotMap`] the entries live in the table itself, so that a probe touches the key
/// and the value together, and they keep no degree, since all of them have the same one.
/// Empty entries hold a sentinel key that no curve class packs to
/// ([`KeyLayout::sentinel`](super::key::KeyLayout::sentinel)). The table remembers which
/// entries are in use, in order of insertion, so that iterating over and clearing it cost the
/// number of entries rather than the capacity. Memory is only allocated on first use, and
/// the table grows when it is three quarters full.
#[derive(Debug, Clone)]
pub struct LevelTable<K, V> {
    entries: Vec<(K, V)>,
    /// The indices of the entries in use, in order of insertion.
    used: Vec<u32>,
    empty: K,
}

impl<K: Key, V: Copy + Default> LevelTable<K, V> {
    /// An empty table, with `empty` marking unused entries.
    pub fn new(empty: K) -> Self {
        Self {
            entries: Vec::new(),
            used: Vec::new(),
            empty,
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

    /// The entry holding `key`, or else the empty entry where it would go.
    #[inline]
    fn probe(&self, key: K) -> Result<usize, usize> {
        let mask = self.entries.len() - 1;
        let mut h = key.hash64() as usize & mask;
        loop {
            let k = self.entries[h].0;
            if k == key {
                return Ok(h);
            }
            if k == self.empty {
                return Err(h);
            }
            h = (h + 1) & mask;
        }
    }

    /// The value of `key`, if it has an entry.
    #[inline]
    pub fn get(&self, key: K) -> Option<&V> {
        if self.entries.is_empty() {
            return None;
        }
        self.probe(key).ok().map(|h| &self.entries[h].1)
    }

    /// The value of `key`, inserted as zero if it has no entry.
    #[inline]
    pub fn get_or_insert(&mut self, key: K) -> &mut V {
        debug_assert!(key != self.empty, "the sentinel key was inserted");
        if self.entries.is_empty() {
            self.entries = vec![(self.empty, V::default()); LEVEL_MIN_CAPACITY];
        }
        let h = match self.probe(key) {
            Ok(h) => h,
            Err(h) => {
                let h = if 4 * (self.used.len() + 1) > 3 * self.entries.len() {
                    self.grow();
                    self.probe(key).unwrap_err()
                } else {
                    h
                };
                self.entries[h] = (key, V::default());
                self.used.push(h as u32);
                h
            }
        };
        &mut self.entries[h].1
    }

    /// Doubles the capacity.
    fn grow(&mut self) {
        let capacity = 2 * self.entries.len();
        assert!(capacity <= 1 << 32, "too many entries for a level table");
        let old = core::mem::replace(
            &mut self.entries,
            vec![(self.empty, V::default()); capacity],
        );
        for i in self.used.iter_mut() {
            let entry = old[*i as usize];
            let mut h = entry.0.hash64() as usize & (capacity - 1);
            while self.entries[h].0 != self.empty {
                h = (h + 1) & (capacity - 1);
            }
            self.entries[h] = entry;
            *i = h as u32;
        }
    }

    /// The entries in use, in order of insertion.
    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = (K, V)> + '_ {
        self.used.iter().map(|&h| self.entries[h as usize])
    }

    /// Removes all entries, keeping the memory.
    pub fn clear(&mut self) {
        for &h in &self.used {
            self.entries[h as usize].0 = self.empty;
        }
        self.used.clear();
    }

    /// Removes all entries and frees the memory.
    pub fn release(&mut self) {
        self.entries = Vec::new();
        self.used = Vec::new();
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

    #[test]
    fn test_take() {
        let mut rng = TestRng::new(6);
        let layout = PackedLayout::new(&[20, 20, 20]).unwrap();
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
        poly.push(PackedKey(1), 1, &[1u64]);
        poly.push(PackedKey(0), 1, &[1u64]);
    }

    #[test]
    fn test_level_table() {
        let mut rng = TestRng::new(8);
        let layout = PackedLayout::new(&[100; 4]).unwrap();
        let mut table: LevelTable<PackedKey, u64> = LevelTable::new(layout.sentinel());
        assert_eq!(table.get(PackedKey(0)), None);
        assert_eq!(table.capacity(), 0);
        for n in [10, 3000, 200, 0, 5000] {
            let mut model: HashMap<PackedKey, u64> = HashMap::new();
            let mut order = Vec::new();
            for _ in 0..n {
                let c: Vec<i32> = (0..4).map(|_| rng.range(-12, 12) as i32).collect();
                let key = layout.pack(&c);
                let v = rng.next_u64();
                let entry = table.get_or_insert(key);
                *entry = entry.wrapping_add(v);
                let m = model.entry(key).or_insert_with(|| {
                    order.push(key);
                    0
                });
                *m = m.wrapping_add(v);
            }
            assert_eq!(table.len(), model.len());
            assert!(4 * table.len() <= 3 * table.capacity());
            let entries: Vec<(PackedKey, u64)> = table.iter().collect();
            let expected: Vec<(PackedKey, u64)> = order.iter().map(|k| (*k, model[k])).collect();
            assert_eq!(entries, expected);
            for (k, v) in &model {
                assert_eq!(table.get(*k), Some(v));
            }
            assert_eq!(table.get(layout.pack(&[50, 50, 50, 50])), None);
            table.clear();
            assert!(table.is_empty());
            assert_eq!(table.get(order.first().copied().unwrap_or_default()), None);
        }
        table.release();
        assert_eq!(table.capacity(), 0);
        *table.get_or_insert(PackedKey(0)) = 3;
        assert_eq!(table.iter().collect::<Vec<_>>(), [(PackedKey(0), 3)]);
    }
}
