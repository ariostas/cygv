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
//! are scattered into it. The level tables key the terms by linear hashes of their offsets
//! $m$ ([`level`](super::level)), so that the hash of $u + v$ is one addition, and only
//! the terms emitted are decoded into curve classes.
//!
//! The curves of each degree are applied in parallel, and the terms they emit are subtracted
//! from the series through shards with a lock each (see [`extract`]).
//!
//! With $\alpha$ along the rows of the charges, the exponentials of the rows of the most terms
//! are instead shared between the curves, along a tree of the curves by their intersections
//! with those rows (`tree`, `PLAN.md` phase 8e), and only the rest is built curve by curve.
//!
//! This stage only takes plain series and tables and returns a series, so that other
//! implementations of it, such as a GPU one, can take its place (`PLAN.md`, section 4.1).

use super::field::Ring;
use super::key::{Key, KeyError, KeyLayout};
use super::level::{
    level_box, BoxBounds, CompactTable, LevelBox, LevelLayout, LevelTable, LevelWord, OffsetKeys,
};
use super::lightcone::Lightcone;
use super::table::{shard_of, DegreeHeap, SlotMap, SparsePoly};
use super::Error;
use rayon::prelude::*;
use std::sync::Mutex;

mod tree;

/// Level tables larger than this are freed after each curve rather than kept for the next.
const LEVEL_KEEP: usize = 1 << 12;

/// How many terms ahead the scatter prefetches the entries it will write.
const PREFETCH: usize = 4;

/// The number of terms a worker buffers for a shard before it subtracts them.
const FLUSH: usize = 1024;

/// How many terms ahead the subtraction of a buffer prefetches the buckets of their keys.
const FOLD_PREFETCH: usize = 8;

/// What the curves pass on instead of their coefficients (see [`extract`]): given a curve
/// and its degree, it writes the values over them.
pub type Source<'a, E> = &'a dyn Fn(&[i32], u32, &mut [E]);

/// When the extraction applies a layer level-synchronously rather than curve by curve (see
/// [`extract`]). The defaults were tuned on a machine with four cores and hyperthreading.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Tuning {
    /// Layers with fewer curves than this per thread are applied level-synchronously when
    /// their curves are large.
    pub sync_max_curves: usize,
    /// Curves are large when $\alpha$ has at least this many terms up to the degree their
    /// exponentials reach.
    pub sync_min_terms: usize,
    /// Whether the level tables pack the offsets into the widest word the keys allow (128
    /// bits, or the width of keys wider than that) even when 64 bits fit, which only the tests
    /// need, since the geometries they can afford fit 64 bits.
    pub wide_levels: bool,
    /// The number of rows of the tree along which an extraction along rows shares the
    /// exponentials of the curves (see [`extract`]); zero builds every exponential curve by
    /// curve instead. `None` chooses the depth from the terms of the rows for a series with
    /// one value per term, and builds the exponentials curve by curve for one with several,
    /// which the tree extracts one component at a time, more slowly than the curves do.
    pub tree_depth: Option<usize>,
    /// The number of degrees whose terms the leaves of the tree take at a time, or `None` for
    /// all of them: the exponential of each curve is then built again for every window of
    /// that many degrees above it, which costs time and keeps the leaves to a fraction of their
    /// terms.
    pub leaf_window: Option<u32>,
    /// Whether the terms of the degrees of the tree that no power reaches beyond, which only
    /// sum up to the root, go to it directly, into tables split by hash that the threads fill
    /// through buffers, rather than into the nodes where they arise. That takes less memory
    /// and somewhat more time. Only for an exact ring, whose sums may be taken in any order.
    pub sinks: bool,
}

impl Tuning {
    /// Whether an extraction along rows of a series with `width` values per term goes along
    /// the tree.
    pub fn uses_tree(&self, width: usize) -> bool {
        match self.tree_depth {
            Some(depth) => depth > 0,
            None => width == 1,
        }
    }
}

impl Default for Tuning {
    fn default() -> Self {
        Self {
            sync_max_curves: 2,
            sync_min_terms: 20000,
            wide_levels: false,
            tree_depth: None,
            leaf_window: None,
            sinks: false,
        }
    }
}

/// The nonzero terms $\deg(u)\, L_u$ of $L = C \cdot \alpha$ up to some degree, in order of
/// degree, for the exponential of one curve, with the hashes of their offsets.
struct LTerms<W, E> {
    hashes: Vec<W>,
    vals: Vec<E>,
    /// The degrees present, each with the range of its terms.
    groups: Vec<(u32, usize, usize)>,
    /// The nonzero coordinates of the curve, with their images.
    coords: Vec<(usize, E)>,
}

impl<W: LevelWord, E: Copy + Default + PartialEq> LTerms<W, E> {
    fn new() -> Self {
        Self {
            hashes: Vec::new(),
            vals: Vec::new(),
            groups: Vec::new(),
            coords: Vec::new(),
        }
    }

    /// Computes the terms of degree at most `top` for `curve`, with `ints` the images of the
    /// integers up to `top`, `alpha` along `rows` if given (see [`extract`]), and the hashes
    /// of the offsets of its terms in `hashes`.
    #[allow(clippy::too_many_arguments)]
    fn build<R: Ring<Elem = E>, K: Key>(
        &mut self,
        ring: &R,
        curve: &[i32],
        alpha: &SparsePoly<K, E>,
        hashes: &[W],
        rows: Option<&[Vec<i32>]>,
        ints: &[E],
        top: u32,
    ) {
        self.coords.clear();
        match rows {
            None => self.coords.extend(
                curve
                    .iter()
                    .enumerate()
                    .filter(|&(_, &c)| c != 0)
                    .map(|(a, &c)| (a, ring.int(c as i64))),
            ),
            Some(rows) => self.coords.extend(
                rows.iter()
                    .map(|row| super::dot(curve, row))
                    .enumerate()
                    .filter(|&(_, x)| x != 0)
                    .map(|(r, x)| (r, ring.int(x))),
            ),
        }
        self.hashes.clear();
        self.vals.clear();
        self.groups.clear();
        for (u, &hash) in hashes.iter().enumerate().take(alpha.len_upto(top)) {
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
                    .push((deg, self.hashes.len(), self.hashes.len() + 1)),
            }
            self.hashes.push(hash);
            self.vals.push(ring.mul(v, ints[deg as usize]));
        }
    }

    /// Adds $\deg(u)\, L_u\, f_v$ to the entry of $u + v$ in `target` for every source
    /// $(v, f_v)$ and every term $u$ of group `q`, the offsets given by their hashes.
    fn scatter<R: Ring<Elem = E>>(
        &self,
        ring: &R,
        q: usize,
        sources: &[(W, E)],
        target: &mut LevelTable<W, E>,
    ) {
        let (_, j0, j1) = self.groups[q];
        let (hashes, vals) = (&self.hashes[j0..j1], &self.vals[j0..j1]);
        for &(hb, fb) in sources {
            for (j, (&hu, &lu)) in hashes.iter().zip(vals).enumerate() {
                if let Some(&hp) = hashes.get(j + PREFETCH) {
                    target.prefetch(hb.wadd(hp));
                }
                ring.fma(target.get_or_insert(hb.wadd(hu)), lu, fb);
            }
        }
    }
}

/// Collects the final terms $f_m$ of degree `e` of an exponential from `level`, where they
/// sit multiplied by `e` (unless `e` is zero), as the nonzero sources for the levels above,
/// and calls `emit(h, f_m)` for each one of positive degree, with `h` the hash of $m$.
///
/// Unless the ring is exact, the sources are sorted by hash, since the order in which a level
/// table holds its terms depends on its capacity, so that the sums the scatters take do not
/// depend on the history of the tables.
fn finalize<R: Ring, W: LevelWord>(
    ring: &R,
    level: &LevelTable<W, R::Elem>,
    e: u32,
    inv: &[R::Elem],
    sources: &mut Vec<(W, R::Elem)>,
    mut emit: impl FnMut(W, R::Elem),
) {
    sources.clear();
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
            emit(key, f);
        }
        sources.push((key, f));
    }
    if !R::EXACT {
        sources.sort_unstable_by_key(|&(h, _)| h);
    }
}

/// The working space for the exponentials, kept from one curve to the next.
struct ExpWorkspace<W, E> {
    /// The terms of each degree, from zero to the maximum degree.
    levels: Vec<LevelTable<W, E>>,
    l: LTerms<W, E>,
    /// The nonzero terms of the current degree.
    sources: Vec<(W, E)>,
}

impl<W: LevelWord, E: Copy + Default + PartialEq> ExpWorkspace<W, E> {
    fn new(marker: W, max_deg: u32) -> Self {
        Self {
            levels: (0..=max_deg).map(|_| LevelTable::new(marker)).collect(),
            l: LTerms::new(),
            sources: Vec::new(),
        }
    }

    /// Calls `emit(h, deg(m), f_m)` for every nonzero term of positive degree at most `top`
    /// of $f = \exp(C \cdot \alpha)$, in order of degree, with `h` the hash of $m$.
    ///
    /// `ints` and `inv` hold the images of the integers and their inverses up to `top`.
    #[allow(clippy::too_many_arguments)]
    fn exp<R: Ring<Elem = E>, K: Key>(
        &mut self,
        ring: &R,
        curve: &[i32],
        alpha: &SparsePoly<K, E>,
        hashes: &[W],
        rows: Option<&[Vec<i32>]>,
        ints: &[E],
        inv: &[E],
        top: u32,
        mut emit: impl FnMut(W, u32, E),
    ) {
        self.l.build(ring, curve, alpha, hashes, rows, ints, top);
        *self.levels[0].get_or_insert(W::default()) = ring.one();
        for e in 0..=top {
            let level = &mut self.levels[e as usize];
            if level.is_empty() {
                continue;
            }
            finalize(ring, level, e, inv, &mut self.sources, |m, f| emit(m, e, f));
            if level.capacity() > LEVEL_KEEP {
                level.release();
            } else {
                level.clear();
            }
            for (q, &(g, _, _)) in self.l.groups.iter().enumerate() {
                if g > top - e {
                    break;
                }
                let target = &mut self.levels[(e + g) as usize];
                self.l.scatter(ring, q, &self.sources, target);
            }
        }
    }
}

/// A shard of the instanton series as the extraction reduces it: the terms whose keys hash to
/// it, and a heap of those not taken yet.
struct Shard<K, E> {
    map: SlotMap<K, E>,
    heap: DegreeHeap,
    /// The final terms taken out of the map, nonzero, each part in order of degree and key.
    done: Vec<SparsePoly<K, E>>,
    /// Space for [`Shard::fold_ordered`]: the terms in order, by curve, buffer and index,
    /// and the hashes of their keys.
    order: Vec<(u32, u32, u32)>,
    hashes: Vec<u64>,
}

impl<K: Key, E: Copy + Default + PartialEq> Shard<K, E> {
    /// Takes the terms of degree below `below` out of the map, which are final once the
    /// curves below it are applied, when they are a ninth of it or more: the nonzero ones to
    /// [`Shard::done`], and the zero ones, most of them, away. Every term left is still in the
    /// heap.
    fn evict(&mut self, below: u32) {
        if 8 * self.map.len() < 9 * self.heap.len() + 4096 {
            return;
        }
        let width = self.map.width();
        let (mut fin, mut keep) = (Vec::new(), Vec::new());
        for s in 0..self.map.len() as u32 {
            if self.map.deg(s) >= below {
                keep.push(s);
            } else if self.map.values(s).iter().any(|v| *v != E::default()) {
                fin.push(s);
            }
        }
        fin.sort_unstable_by_key(|&s| (self.map.deg(s), self.map.key(s)));
        let mut part = SparsePoly::new(width);
        for &s in &fin {
            part.push(self.map.key(s), self.map.deg(s), self.map.values(s));
        }
        if !part.is_empty() {
            self.done.push(part);
        }
        let mut map = SlotMap::new(width, keep.len() + keep.len() / 4 + 16);
        let mut heap = DegreeHeap::new();
        for &s in &keep {
            let (t, _) = map.insert(self.map.key(s), self.map.deg(s));
            map.values_mut(t).copy_from_slice(self.map.values(s));
            heap.push(map.deg(t), t);
        }
        (self.map, self.heap) = (map, heap);
    }

    /// The nonzero terms, in parts.
    fn finish(mut self) -> Vec<SparsePoly<K, E>> {
        let rest = self.map.take();
        self.done.push(rest);
        self.done
    }

    /// Subtracts the terms of `out` and empties it.
    ///
    /// The keys are hashed first, so that the bucket of each one can be prefetched a few
    /// terms ahead of its insertion.
    fn fold<R: Ring<Elem = E>>(&mut self, ring: &R, out: &mut OutBuf<K, E>) {
        let width = self.map.width();
        out.hashes.clear();
        out.hashes.extend(out.keys.iter().map(|k| k.hash64()));
        for &h in out.hashes.iter().take(FOLD_PREFETCH) {
            self.map.prefetch(h);
        }
        for (i, (&key, &deg)) in out.keys.iter().zip(&out.degs).enumerate() {
            if let Some(&h) = out.hashes.get(i + FOLD_PREFETCH) {
                self.map.prefetch(h);
            }
            let (t, created) = self.map.insert_hashed(key, out.hashes[i], deg);
            if created {
                self.heap.push(deg, t);
            }
            let v = &out.vals[i * width..(i + 1) * width];
            for (x, &y) in self.map.values_mut(t).iter_mut().zip(v) {
                *x = ring.sub(*x, y);
            }
        }
        out.clear();
    }

    /// Subtracts the terms of `outs`, which record the curves that emitted them, in the order
    /// of the curves, and empties them. A curve emits at most one term per key, so every key
    /// takes its terms in the same order however they were spread over the buffers.
    fn fold_ordered<R: Ring<Elem = E>>(&mut self, ring: &R, outs: &mut [&mut OutBuf<K, E>]) {
        let width = self.map.width();
        self.order.clear();
        for (b, out) in outs.iter().enumerate() {
            let terms = out.srcs.iter().enumerate();
            self.order
                .extend(terms.map(|(i, &src)| (src, b as u32, i as u32)));
        }
        self.order.sort_unstable_by_key(|&(src, _, _)| src);
        self.hashes.clear();
        self.hashes.extend(
            self.order
                .iter()
                .map(|&(_, b, i)| outs[b as usize].keys[i as usize].hash64()),
        );
        for &h in self.hashes.iter().take(FOLD_PREFETCH) {
            self.map.prefetch(h);
        }
        for (j, &(_, b, i)) in self.order.iter().enumerate() {
            if let Some(&h) = self.hashes.get(j + FOLD_PREFETCH) {
                self.map.prefetch(h);
            }
            let (out, i) = (&outs[b as usize], i as usize);
            let deg = out.degs[i];
            let (t, created) = self.map.insert_hashed(out.keys[i], self.hashes[j], deg);
            if created {
                self.heap.push(deg, t);
            }
            let v = &out.vals[i * width..(i + 1) * width];
            for (x, &y) in self.map.values_mut(t).iter_mut().zip(v) {
                *x = ring.sub(*x, y);
            }
        }
        outs.iter_mut().for_each(|o| o.clear());
    }
}

/// The shards of the series `inst`, each with a heap of its terms.
fn new_shards<K: Key, E: Copy + Default + PartialEq>(
    inst: &SparsePoly<K, E>,
    n_shards: usize,
) -> Vec<Mutex<Shard<K, E>>> {
    let mut shards: Vec<Mutex<Shard<K, E>>> = (0..n_shards)
        .map(|_| {
            Mutex::new(Shard {
                map: SlotMap::new(inst.width(), (inst.len() + inst.len() / 8) / n_shards + 16),
                heap: DegreeHeap::new(),
                done: Vec::new(),
                order: Vec::new(),
                hashes: Vec::new(),
            })
        })
        .collect();
    for i in 0..inst.len() {
        let key = inst.keys()[i];
        let shard = shards[shard_of(key, n_shards)].get_mut().unwrap();
        let (s, _) = shard.map.insert(key, inst.degs()[i]);
        shard.map.values_mut(s).copy_from_slice(inst.values(i));
        shard.heap.push(inst.degs()[i], s);
    }
    shards
}

/// Terms to subtract from one shard, as a worker collects them.
struct OutBuf<K, E> {
    keys: Vec<K>,
    degs: Vec<u32>,
    vals: Vec<E>,
    /// The index in its group of the curve that emitted each term, when the sums must be
    /// taken in order (see [`extract`]).
    srcs: Vec<u32>,
    /// Space for the hashes of the keys, when they are subtracted.
    hashes: Vec<u64>,
}

impl<K, E> OutBuf<K, E> {
    fn new() -> Self {
        Self {
            keys: Vec::new(),
            degs: Vec::new(),
            vals: Vec::new(),
            srcs: Vec::new(),
            hashes: Vec::new(),
        }
    }

    fn clear(&mut self) {
        self.keys.clear();
        self.degs.clear();
        self.vals.clear();
        self.srcs.clear();
    }
}

/// Buffers $n_C f$ to be subtracted at `key`, of degree `deg`, in the buffer of its shard in
/// `out`, and subtracts the buffer if it is full. With `src`, the index of the curve, the
/// buffer records it and is only subtracted at the end of the group, in order.
#[inline]
#[allow(clippy::too_many_arguments)]
fn emit<R: Ring, K: Key>(
    ring: &R,
    shards: &[Mutex<Shard<K, R::Elem>>],
    out: &mut [OutBuf<K, R::Elem>],
    src: Option<u32>,
    key: K,
    deg: u32,
    n_c: &[R::Elem],
    f: R::Elem,
) {
    let s = shard_of(key, shards.len());
    let o = &mut out[s];
    o.keys.push(key);
    o.degs.push(deg);
    o.vals.extend(n_c.iter().map(|&n| ring.mul(n, f)));
    if let Some(src) = src {
        o.srcs.push(src);
    } else if o.keys.len() >= FLUSH {
        shards[s].lock().unwrap().fold(ring, o);
    }
}

/// The most terms that the curves of a group of a layer are expected to emit in all, as a
/// multiple of the terms of the instanton series, when they are kept until the group is done
/// (see [`extract`]).
const GROUP_TERMS: usize = 2;

/// The most curves in a group of a layer.
const GROUP_CURVES: usize = 4096;

/// Sorts the curves of a layer, with their values, `width` each, by key.
fn sort_layer<K: Key, E: Copy>(layer: &mut Vec<(K, u32)>, vals: &mut Vec<E>, width: usize) {
    let mut order: Vec<usize> = (0..layer.len()).collect();
    order.sort_unstable_by_key(|&i| layer[i].0);
    *vals = order
        .iter()
        .flat_map(|&i| vals[i * width..(i + 1) * width].iter().copied())
        .collect();
    *layer = order.iter().map(|&i| layer[i]).collect();
}

/// The end of the group of the curves of `layer` from `start` on: at most [`GROUP_CURVES`],
/// whose exponentials have about [`GROUP_TERMS`] times as many terms as `inst` in all, by the
/// terms of `inst` up to the degree each one reaches, but at least one curve.
fn group_end<K: Key>(layer: &[(K, u32)], start: usize, upto: &[usize], max_deg: u32) -> usize {
    let budget = GROUP_TERMS * upto[max_deg as usize].max(1 << 16);
    let mut end = start + 1;
    let mut terms = upto[(max_deg - layer[start].1) as usize];
    while end < layer.len() && end - start < GROUP_CURVES {
        terms += upto[(max_deg - layer[end].1) as usize];
        if terms > budget {
            break;
        }
        end += 1;
    }
    end
}

/// The result of an extraction.
type Extracted<K, E> = Result<SparsePoly<K, E>, Error>;

/// The workers of an extraction, one per thread, each behind a lock that only its thread
/// takes.
type Workers<W, K, E> = [Mutex<Worker<W, K, E>>];

/// What a worker thread keeps from one curve to the next.
struct Worker<W, K, E> {
    ws: ExpWorkspace<W, E>,
    curve: Vec<i32>,
    /// One buffer per shard.
    out: Vec<OutBuf<K, E>>,
}

/// A curve whose exponential is built together with those of the other curves of its layer,
/// one degree at a time, so that the scatters into the levels of every curve run in parallel.
struct SyncCurve<W, K, E> {
    key: K,
    deg: u32,
    /// The highest degree of its exponential.
    top: u32,
    l: LTerms<W, E>,
    /// The terms of each degree up to `top`, each locked by the one scatter that writes it.
    levels: Vec<Mutex<LevelTable<W, E>>>,
    sources: Vec<(W, E)>,
    /// One buffer per shard.
    out: Vec<OutBuf<K, E>>,
}

/// What the curves of a layer share while they are applied.
struct LayerCtx<'a, R: Ring, L: KeyLayout, W> {
    ring: &'a R,
    layout: &'a L,
    level: &'a LevelLayout<W>,
    /// The keys of the offsets, from their hashes.
    keys: &'a OffsetKeys<W, L::Key>,
    alpha: &'a SparsePoly<L::Key, R::Elem>,
    /// The hashes of the offsets of the terms of `alpha`.
    hashes: &'a [W],
    rows: Option<&'a [Vec<i32>]>,
    ints: &'a [R::Elem],
    inv: &'a [R::Elem],
    max_deg: u32,
    shards: &'a [Mutex<Shard<L::Key, R::Elem>>],
    /// Whether the terms record their curves, to be subtracted in order (see [`extract`]).
    ordered: bool,
    /// The curves, each with its coefficients, `width` values each.
    curves: &'a [(L::Key, u32)],
    vals: &'a [R::Elem],
    width: usize,
}

impl<R: Ring, L: KeyLayout, W: LevelWord> LayerCtx<'_, R, L, W> {
    fn n_c(&self, i: usize) -> &[R::Elem] {
        &self.vals[i * self.width..(i + 1) * self.width]
    }

    /// What the terms of curve `i` record.
    fn src(&self, i: usize) -> Option<u32> {
        self.ordered.then_some(i as u32)
    }

    /// Applies each curve on one thread, with the workers indexed by the threads.
    fn apply_by_curve(&self, workers: &Workers<W, L::Key, R::Elem>) {
        let apply = |i: usize| {
            let (kc, dc) = self.curves[i];
            let (n_c, src) = (self.n_c(i), self.src(i));
            let t = rayon::current_thread_index().unwrap_or(0);
            let mut worker = workers[t].lock().unwrap();
            let Worker { ws, curve, out } = &mut *worker;
            self.layout.unpack(kc, curve);
            ws.exp(
                self.ring,
                curve,
                self.alpha,
                self.hashes,
                self.rows,
                self.ints,
                self.inv,
                self.max_deg - dc,
                |h, e, f| {
                    let m = self.keys.key(h, e);
                    emit(self.ring, self.shards, out, src, kc + m, dc + e, n_c, f)
                },
            );
        };
        if workers.len() == 1 {
            (0..self.curves.len()).for_each(apply);
        } else {
            (0..self.curves.len()).into_par_iter().for_each(apply);
        }
    }

    /// Applies the curves level-synchronously: for each degree, every curve finalizes its
    /// terms of that degree, and then the scatters of every curve into each of its higher
    /// levels run as separate tasks, the largest first.
    fn apply_by_level(&self) -> Vec<SyncCurve<W, L::Key, R::Elem>> {
        let ring = self.ring;
        let marker = self.level.marker();
        let mut curves: Vec<SyncCurve<W, L::Key, R::Elem>> = self
            .curves
            .par_iter()
            .map(|&(key, deg)| {
                let top = self.max_deg - deg;
                let mut curve = vec![0; self.layout.dim()];
                self.layout.unpack(key, &mut curve);
                let mut l = LTerms::new();
                l.build(
                    ring,
                    &curve,
                    self.alpha,
                    self.hashes,
                    self.rows,
                    self.ints,
                    top,
                );
                let levels: Vec<_> = (0..=top)
                    .map(|_| Mutex::new(LevelTable::new(marker)))
                    .collect();
                *levels[0].lock().unwrap().get_or_insert(W::default()) = ring.one();
                SyncCurve {
                    key,
                    deg,
                    top,
                    l,
                    levels,
                    sources: Vec::new(),
                    out: (0..self.shards.len()).map(|_| OutBuf::new()).collect(),
                }
            })
            .collect();
        let max_top = curves.iter().map(|c| c.top).max().unwrap_or(0);
        let mut tasks: Vec<(usize, usize, usize)> = Vec::new();
        for e in 0..=max_top {
            curves.par_iter_mut().enumerate().for_each(|(i, c)| {
                c.sources.clear();
                if e > c.top {
                    return;
                }
                let level = c.levels[e as usize].get_mut().unwrap();
                if level.is_empty() {
                    return;
                }
                let (n_c, src) = (self.n_c(i), self.src(i));
                let (kc, dc) = (c.key, c.deg);
                let out = &mut c.out;
                finalize(ring, level, e, self.inv, &mut c.sources, |h, f| {
                    let m = self.keys.key(h, e);
                    emit(ring, self.shards, out, src, kc + m, dc + e, n_c, f)
                });
                // This level is done.
                *level = LevelTable::new(marker);
            });
            tasks.clear();
            for (i, c) in curves.iter().enumerate() {
                if c.sources.is_empty() {
                    continue;
                }
                for (q, &(g, j0, j1)) in c.l.groups.iter().enumerate() {
                    if g > c.top - e {
                        break;
                    }
                    tasks.push((c.sources.len() * (j1 - j0), i, q));
                }
            }
            // Largest first, and one task each, since rayon would otherwise run them in a few
            // runs in order, the first of which would take the largest one after the other.
            tasks.sort_unstable_by(|a, b| b.cmp(a));
            tasks.par_iter().with_max_len(1).for_each(|&(_, i, q)| {
                let c = &curves[i];
                let g = c.l.groups[q].0;
                let mut target = c.levels[(e + g) as usize].lock().unwrap();
                c.l.scatter(ring, q, &c.sources, &mut target);
            });
        }
        curves
    }
}

/// Extracts the coefficients $n_C$ of the instanton series `inst` in the basis of the
/// $q^C$, at every curve class of degree at most `max_deg` under `grading`, as a series with
/// as many values per term as `inst`.
///
/// `alpha` holds the $\alpha_a$, one value per coordinate of `layout`, and must have no
/// constant term. With `rows`, it instead holds $\beta$ along those rows $M_r$, one value per
/// row, with $\alpha_a = \sum_r M_{ra} \beta_r$, so that $C \cdot \alpha = \sum_r (M_r \cdot C)
/// \beta_r$: the same exactly, but for the majorant run each term $(M_r \cdot C) \beta_r$ is
/// bounded by itself, and vanishes with $M_r \cdot C$, rather than through
/// $\sum_a \lvert C_a\rvert \lvert\alpha_a\rvert$; `inv` holds the inverses of the integers up to `max_deg` at least. Only the
/// curves in `lightcone` are subtracted, so the coefficients are only right in there.
///
/// Along rows, the exponentials go along a tree (see the module `tree`) unless `tuning` says
/// otherwise ([`Tuning::uses_tree`]): $\exp(C \cdot \alpha) = \prod_r E_r^{M_r \cdot C}$ with
/// $E_r = \exp \beta_r$, and the powers of the rows of the most terms multiply the sums of the
/// curves that share their exponents, once for all of them. What follows describes the
/// exponentials built curve by curve, which the tree also builds along the other rows.
///
/// With `source`, each curve's $n_C$ is still what the result holds, but what is subtracted
/// is $q^C$ times the values `source` writes over it, given the curve and its degree. The
/// certifier's majorant run passes the values the candidate invariants imply, so that the
/// bounds do not compound from degree to degree.
///
/// It runs on the threads of the current rayon pool. Every term of $\exp(C \cdot \alpha)$
/// but the first has at least the lowest degree $\delta$ in the support of $\alpha$, so the
/// curves of degree in $[d, d + \delta)$ are all final once those of degree below $d$ are
/// subtracted, and they form a layer whose curves can be applied in parallel. Usually each
/// thread applies whole curves; a layer of few curves with large exponentials, such as the
/// first ones, instead builds them level-synchronously, with the scatters of every curve into
/// each of its levels as tasks. The series is split into shards by key, so that the terms the
/// threads emit can be subtracted under one lock per shard: each thread buffers them per
/// shard and subtracts a buffer when it is full, and the shards take the remaining buffers in
/// parallel at the end of the layer.
///
/// Unless the ring is exact ([`Ring::EXACT`]), the sums must not depend on the scheduling:
/// the curves of a layer are sorted by key and applied in groups whose exponentials have
/// about twice as many terms as `inst` in all (at most 4096 curves), the buffers recording
/// which curve emitted each term and kept until the group is done, when each shard takes its
/// terms in the order of the curves. Each curve emits at most one term per key, so every key
/// then takes its terms in the same order whatever the threads.
///
/// The level tables pack the offsets into 64 bits when the box of [`level_box`] fits, into
/// 128 bits otherwise, and for keys wider than 128 bits into a word as wide as theirs if that
/// is not enough either ([`Key::Level`]; [`Tuning::wide_levels`] forces the widest).
#[allow(clippy::too_many_arguments)]
pub fn extract<R: Ring, L: KeyLayout>(
    ring: &R,
    layout: &L,
    grading: &[i32],
    alpha: SparsePoly<L::Key, R::Elem>,
    rows: Option<&[Vec<i32>]>,
    inst: SparsePoly<L::Key, R::Elem>,
    inv: &[R::Elem],
    max_deg: u32,
    lightcone: &Lightcone<L::Key>,
    tuning: &Tuning,
    source: Option<Source<'_, R::Elem>>,
) -> Extracted<L::Key, R::Elem> {
    let h11 = layout.dim();
    let n_alpha = rows.map_or(h11, |r| r.len());
    if alpha.width() != n_alpha
        || grading.len() != h11
        || rows.is_some_and(|r| r.iter().any(|row| row.len() != h11))
        || inv.len() <= max_deg as usize
    {
        return Err(Error::DimensionMismatch);
    }
    assert!(
        alpha.degs().first().is_none_or(|&d| d > 0),
        "alpha has a constant term"
    );
    if let Some(rows) = rows.filter(|_| tuning.uses_tree(inst.width())) {
        return tree::extract_tree(
            ring, layout, grading, alpha, rows, inst, inv, max_deg, lightcone, tuning, source,
        );
    }
    // The narrowest level word that fits.
    let n_terms = alpha.len_upto(max_deg);
    let level_box = offset_box(layout, grading, &alpha, n_terms, max_deg).1;
    let wide = <<L::Key as Key>::Level as LevelWord>::BITS > u128::BITS;
    macro_rules! go {
        ($w:ty) => {
            extract_in::<R, L, $w>(
                ring, layout, grading, &alpha, rows, inst, inv, max_deg, lightcone, tuning, source,
            )
            .expect("the level layout fits")
        };
    }
    if !tuning.wide_levels && LevelLayout::<u64>::new(&level_box, grading).is_some() {
        return go!(u64);
    }
    if (!tuning.wide_levels || !wide) && LevelLayout::<u128>::new(&level_box, grading).is_some() {
        return go!(u128);
    }
    if wide && LevelLayout::<<L::Key as Key>::Level>::new(&level_box, grading).is_some() {
        return go!(<L::Key as Key>::Level);
    }
    Err(Error::Key(KeyError::TooWide {
        bits: level_box.bits(),
        max_bits: <<L::Key as Key>::Level as LevelWord>::BITS - 1,
    }))
}

/// The coordinates of the terms of `alpha` up to the `n_terms`-th, one after the other, and
/// the box of their sums up to `max_deg`.
fn offset_box<L: KeyLayout, E: Copy + Default + PartialEq>(
    layout: &L,
    grading: &[i32],
    alpha: &SparsePoly<L::Key, E>,
    n_terms: usize,
    max_deg: u32,
) -> (Vec<i32>, LevelBox) {
    let h11 = layout.dim();
    let mut terms = vec![0; n_terms * h11];
    for (u, t) in terms.chunks_exact_mut(h11.max(1)).enumerate() {
        layout.unpack(alpha.keys()[u], t);
    }
    let level_box = level_box(
        (0..n_terms).map(|u| (&terms[u * h11..(u + 1) * h11], alpha.degs()[u])),
        grading,
        max_deg,
    );
    (terms, level_box)
}

/// [`extract`] with the offsets packed into `W`, or `None` if they do not fit.
#[allow(clippy::too_many_arguments)]
fn extract_in<R: Ring, L: KeyLayout, W: LevelWord>(
    ring: &R,
    layout: &L,
    grading: &[i32],
    alpha: &SparsePoly<L::Key, R::Elem>,
    rows: Option<&[Vec<i32>]>,
    inst: SparsePoly<L::Key, R::Elem>,
    inv: &[R::Elem],
    max_deg: u32,
    lightcone: &Lightcone<L::Key>,
    tuning: &Tuning,
    source: Option<Source<'_, R::Elem>>,
) -> Option<Extracted<L::Key, R::Elem>> {
    let h11 = layout.dim();
    let n_terms = alpha.len_upto(max_deg);
    let (terms, level_box) = offset_box(layout, grading, alpha, n_terms, max_deg);
    let level = LevelLayout::<W>::new(&level_box, grading)?;
    let hashes: Vec<W> = terms
        .chunks_exact(h11.max(1))
        .take(n_terms)
        .map(|t| level.hash(t))
        .collect();
    // A coordinate that is always zero has no room for a one, and needs no unit.
    let units: Vec<L::Key> = (0..h11)
        .map(|t| {
            let mut e = vec![0; h11];
            e[t] = 1;
            layout.try_pack(&e).unwrap_or_default()
        })
        .collect();
    let keys = level.offset_keys(&units);
    Some(extract_with(
        ring, layout, &level, &keys, &hashes, alpha, rows, inst, inv, max_deg, lightcone, tuning,
        source,
    ))
}

/// [`extract`] with the level tables of `level`, and `hashes` the hashes of the offsets of
/// the terms of `alpha` up to `max_deg`.
#[allow(clippy::too_many_arguments)]
fn extract_with<R: Ring, L: KeyLayout, W: LevelWord>(
    ring: &R,
    layout: &L,
    level: &LevelLayout<W>,
    keys: &OffsetKeys<W, L::Key>,
    hashes: &[W],
    alpha: &SparsePoly<L::Key, R::Elem>,
    rows: Option<&[Vec<i32>]>,
    inst: SparsePoly<L::Key, R::Elem>,
    inv: &[R::Elem],
    max_deg: u32,
    lightcone: &Lightcone<L::Key>,
    tuning: &Tuning,
    source: Option<Source<'_, R::Elem>>,
) -> Extracted<L::Key, R::Elem> {
    let h11 = layout.dim();
    let width = inst.width();
    let ints: Vec<R::Elem> = (0..=max_deg).map(|n| ring.int(n as i64)).collect();
    let threads = rayon::current_num_threads();
    let n_shards = if threads == 1 { 1 } else { 4 * threads };

    let mut shards = new_shards(&inst, n_shards);
    // The shards hold the series now; the groups only need its number of terms by degree.
    let upto: Vec<usize> = (0..=max_deg).map(|d| inst.len_upto(d)).collect();
    drop(inst);
    let mut workers: Box<Workers<W, L::Key, R::Elem>> = (0..threads)
        .map(|_| {
            Mutex::new(Worker {
                ws: ExpWorkspace::new(level.marker(), max_deg),
                curve: vec![0; h11],
                out: (0..n_shards).map(|_| OutBuf::new()).collect(),
            })
        })
        .collect();

    let delta = alpha.degs().first().map_or(1, |&d| d.max(1));
    let mut curve = vec![0; h11];
    let mut layer: Vec<(L::Key, u32)> = Vec::new();
    let mut layer_vals: Vec<R::Elem> = Vec::new();
    while let Some(d) = shards
        .iter_mut()
        .filter_map(|s| s.get_mut().unwrap().heap.min_deg())
        .min()
    {
        let top = d.saturating_add(delta);
        layer.clear();
        layer_vals.clear();
        for shard in shards.iter_mut() {
            let shard = shard.get_mut().unwrap();
            while shard.heap.min_deg().is_some_and(|e| e < top) {
                let slot = shard.heap.pop().unwrap();
                let (kc, dc) = (shard.map.key(slot), shard.map.deg(slot));
                let n_c = shard.map.values(slot);
                // A curve of the maximum degree has nothing left to change.
                if dc == 0 || dc >= max_deg || n_c.iter().all(|&x| ring.is_zero(x)) {
                    continue;
                }
                layout.unpack(kc, &mut curve);
                // A curve outside the lightcones only changes curves outside them.
                if !lightcone.contains(kc, &curve) {
                    continue;
                }
                let start = layer_vals.len();
                layer_vals.extend_from_slice(n_c);
                if let Some(source) = source {
                    let vals = &mut layer_vals[start..];
                    source(&curve, dc, vals);
                    if vals.iter().all(|&x| ring.is_zero(x)) {
                        layer_vals.truncate(start);
                        continue;
                    }
                }
                layer.push((kc, dc));
            }
        }
        if layer.is_empty() {
            continue;
        }
        if !R::EXACT {
            sort_layer(&mut layer, &mut layer_vals, width);
        }
        let mut start = 0;
        while start < layer.len() {
            let end = if R::EXACT {
                layer.len()
            } else {
                group_end(&layer, start, &upto, max_deg)
            };
            let ctx = LayerCtx {
                ring,
                layout,
                level,
                keys,
                alpha,
                hashes,
                rows,
                ints: &ints,
                inv,
                max_deg,
                shards: &shards,
                ordered: !R::EXACT,
                curves: &layer[start..end],
                vals: &layer_vals[start * width..end * width],
                width,
            };
            let by_level = threads > 1
                && end - start < tuning.sync_max_curves.saturating_mul(threads)
                && alpha.len_upto(max_deg.saturating_sub(d)) >= tuning.sync_min_terms;
            let mut synced = if by_level {
                ctx.apply_by_level()
            } else {
                ctx.apply_by_curve(&workers);
                Vec::new()
            };

            let mut outs: Vec<Vec<&mut OutBuf<L::Key, R::Elem>>> =
                (0..n_shards).map(|_| Vec::new()).collect();
            let worker_outs = workers.iter_mut().map(|w| &mut w.get_mut().unwrap().out);
            for out in worker_outs.chain(synced.iter_mut().map(|c| &mut c.out)) {
                for (s, o) in out.iter_mut().enumerate() {
                    outs[s].push(o);
                }
            }
            shards
                .par_iter_mut()
                .zip(outs)
                .for_each(|(shard, mut outs)| {
                    let shard = shard.get_mut().unwrap();
                    if R::EXACT {
                        for o in outs {
                            shard.fold(ring, o);
                        }
                    } else {
                        shard.fold_ordered(ring, &mut outs);
                    }
                });
            start = end;
        }
        shards
            .par_iter_mut()
            .for_each(|s| s.get_mut().unwrap().evict(top));
    }
    let parts = shards
        .into_iter()
        .flat_map(|s| s.into_inner().unwrap().finish())
        .collect();
    Ok(SparsePoly::merge(width, parts))
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
        let layout = PackedLayout::<u128>::new(&[30, 30]).unwrap();
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
            let got = extract(
                &lf,
                &layout,
                &[1, 1],
                alpha.clone(),
                None,
                inst.clone(),
                &inv,
                max_deg,
                &Lightcone::All,
                &Tuning::default(),
                None,
            )
            .unwrap();
            let got: HashMap<PackedKey, Vec<E>> = (0..got.len())
                .map(|i| (got.keys()[i], got.values(i).to_vec()))
                .collect();
            let expected: HashMap<PackedKey, Vec<E>> =
                n.iter().map(|(c, v)| (layout.pack(c), v.clone())).collect();
            assert_eq!(got, expected, "width {width}");

            // The same with the offsets in 128 bits.
            let wide = extract_in::<_, _, u128>(
                &lf,
                &layout,
                &[1, 1],
                &alpha,
                None,
                inst.clone(),
                &inv,
                max_deg,
                &Lightcone::All,
                &Tuning::default(),
                None,
            )
            .unwrap()
            .unwrap();
            let wide: HashMap<PackedKey, Vec<E>> = (0..wide.len())
                .map(|i| (wide.keys()[i], wide.values(i).to_vec()))
                .collect();
            assert_eq!(wide, expected, "width {width}, 128 bits");

            // A source that gives each curve its own coefficient changes nothing, and one that
            // gives zero leaves the series as it was.
            let with = |source: Source<'_, E>| {
                let got = extract(
                    &lf,
                    &layout,
                    &[1, 1],
                    alpha.clone(),
                    None,
                    inst.clone(),
                    &inv,
                    max_deg,
                    &Lightcone::All,
                    &Tuning::default(),
                    Some(source),
                )
                .unwrap();
                (0..got.len())
                    .map(|i| (got.keys()[i], got.values(i).to_vec()))
                    .collect::<HashMap<PackedKey, Vec<E>>>()
            };
            let own = |c: &[i32], _: u32, v: &mut [E]| v.copy_from_slice(&n[&[c[0], c[1]]]);
            assert_eq!(with(&own), expected, "width {width}");
            let unchanged: HashMap<PackedKey, Vec<E>> = (0..inst.len())
                .map(|i| (inst.keys()[i], inst.values(i).to_vec()))
                .collect();
            let zero = |_: &[i32], _: u32, v: &mut [E]| v.fill(E::default());
            assert_eq!(with(&zero), unchanged, "width {width}");

            // The same along rows: beta on three rows, and alpha_a = sum_r M_ra beta_r.
            let rows = vec![vec![1, 0], vec![2, -1], vec![0, 3]];
            let mut beta_map = SlotMap::new(3, 16);
            let mut alpha_map = SlotMap::new(2, 16);
            for u in 0..alpha.len() {
                let (k, d) = (alpha.keys()[u], alpha.degs()[u]);
                let b: Vec<E> = (0..3).map(|_| lf.int(rng.range(-5, 5))).collect();
                beta_map.get_or_insert(k, d).copy_from_slice(&b);
                for (a, v) in alpha_map.get_or_insert(k, d).iter_mut().enumerate() {
                    for (row, &x) in rows.iter().zip(&b) {
                        lf.fma(v, lf.int(row[a] as i64), x);
                    }
                }
            }
            let (beta, alpha_m) = (beta_map.take(), alpha_map.take());
            let run = |a: &SparsePoly<PackedKey, E>, r: Option<&[Vec<i32>]>, tuning: &Tuning| {
                let got = extract(
                    &lf,
                    &layout,
                    &[1, 1],
                    a.clone(),
                    r,
                    inst.clone(),
                    &inv,
                    max_deg,
                    &Lightcone::All,
                    tuning,
                    None,
                )
                .unwrap();
                (0..got.len())
                    .map(|i| (got.keys()[i], got.values(i).to_vec()))
                    .collect::<HashMap<PackedKey, Vec<E>>>()
            };
            let expected_rows = run(&alpha_m, None, &Tuning::default());
            // Curve by curve, and along trees of every depth, in 64 and 128 bits, with the
            // leaves taking every degree or windows of them, and with sinks or without.
            for tree_depth in [None, Some(0), Some(1), Some(2), Some(3)] {
                for (wide_levels, leaf_window, sinks) in [
                    (false, None, false),
                    (true, None, false),
                    (false, Some(1), false),
                    (false, None, true),
                    (false, Some(2), true),
                    (true, Some(1), true),
                ] {
                    let tuning = Tuning {
                        tree_depth,
                        wide_levels,
                        leaf_window,
                        sinks,
                        ..Tuning::default()
                    };
                    assert_eq!(
                        run(&beta, Some(&rows), &tuning),
                        expected_rows,
                        "width {width}, {tuning:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_errors() {
        let lf = ring();
        let layout = PackedLayout::<u128>::new(&[30, 30]).unwrap();
        let inv = lf.inverses(5).unwrap();
        let empty = SparsePoly::new(2);
        let inst = SparsePoly::new(1);
        assert_eq!(
            extract(
                &lf,
                &layout,
                &[1, 1],
                empty.clone(),
                None,
                inst.clone(),
                &inv,
                5,
                &Lightcone::All,
                &Tuning::default(),
                None,
            ),
            Err(Error::DimensionMismatch)
        );
        assert_eq!(
            extract(
                &lf,
                &layout,
                &[1, 1],
                SparsePoly::new(3),
                None,
                inst.clone(),
                &inv,
                4,
                &Lightcone::All,
                &Tuning::default(),
                None,
            ),
            Err(Error::DimensionMismatch)
        );
        assert!(extract(
            &lf,
            &layout,
            &[1, 1],
            empty.clone(),
            None,
            inst.clone(),
            &inv,
            4,
            &Lightcone::All,
            &Tuning::default(),
            None,
        )
        .unwrap()
        .is_empty());
    }
}
