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
//! The curves of each degree are applied in parallel, and the terms they emit are subtracted
//! from the series through shards with a lock each (see [`extract`]).
//!
//! This stage only takes plain series and tables and returns a series, so that other
//! implementations of it, such as a GPU one, can take its place (`PLAN.md`, section 4.1).

use super::field::Ring;
use super::key::{Key, KeyLayout};
use super::lightcone::Lightcone;
use super::table::{shard_of, DegreeHeap, LevelTable, SlotMap, SparsePoly};
use super::Error;
use rayon::prelude::*;
use std::sync::Mutex;

/// Level tables larger than this are freed after each curve rather than kept for the next.
const LEVEL_KEEP: usize = 1 << 15;

/// How many terms ahead the scatter prefetches the entries it will write.
const PREFETCH: usize = 4;

/// The number of terms a worker buffers for a shard before it subtracts them.
const FLUSH: usize = 1024;

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
}

impl Default for Tuning {
    fn default() -> Self {
        Self {
            sync_max_curves: 2,
            sync_min_terms: 20000,
        }
    }
}

/// The nonzero terms $\deg(u)\, L_u$ of $L = C \cdot \alpha$ up to some degree, in order of
/// degree, for the exponential of one curve.
struct LTerms<K, E> {
    keys: Vec<K>,
    vals: Vec<E>,
    /// The degrees present, each with the range of its terms.
    groups: Vec<(u32, usize, usize)>,
    /// The nonzero coordinates of the curve, with their images.
    coords: Vec<(usize, E)>,
}

impl<K: Key, E: Copy + Default + PartialEq> LTerms<K, E> {
    fn new() -> Self {
        Self {
            keys: Vec::new(),
            vals: Vec::new(),
            groups: Vec::new(),
            coords: Vec::new(),
        }
    }

    /// Computes the terms of degree at most `top` for `curve`, with `ints` the images of the
    /// integers up to `top`.
    fn build<R: Ring<Elem = E>>(
        &mut self,
        ring: &R,
        curve: &[i32],
        alpha: &SparsePoly<K, E>,
        ints: &[E],
        top: u32,
    ) {
        self.coords.clear();
        self.coords.extend(
            curve
                .iter()
                .enumerate()
                .filter(|&(_, &c)| c != 0)
                .map(|(a, &c)| (a, ring.int(c as i64))),
        );
        self.keys.clear();
        self.vals.clear();
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
                    .push((deg, self.keys.len(), self.keys.len() + 1)),
            }
            self.keys.push(alpha.keys()[u]);
            self.vals.push(ring.mul(v, ints[deg as usize]));
        }
    }

    /// Adds $\deg(u)\, L_u\, f_v$ to the entry of $u + v$ in `target` for every source
    /// $(v, f_v)$ and every term $u$ of group `q`.
    fn scatter<R: Ring<Elem = E>>(
        &self,
        ring: &R,
        q: usize,
        sources: &[(K, E)],
        target: &mut LevelTable<K, E>,
    ) {
        let (_, j0, j1) = self.groups[q];
        let (keys, vals) = (&self.keys[j0..j1], &self.vals[j0..j1]);
        for &(kb, fb) in sources {
            for (j, (&ku, &lu)) in keys.iter().zip(vals).enumerate() {
                if let Some(&kp) = keys.get(j + PREFETCH) {
                    target.prefetch(kb + kp);
                }
                ring.fma(target.get_or_insert(kb + ku), lu, fb);
            }
        }
    }
}

/// Collects the final terms $f_m$ of degree `e` of an exponential from `level`, where they
/// sit multiplied by `e` (unless `e` is zero), as the nonzero sources for the levels above,
/// and calls `emit(m, f_m)` for each one of positive degree.
fn finalize<R: Ring, K: Key>(
    ring: &R,
    level: &LevelTable<K, R::Elem>,
    e: u32,
    inv: &[R::Elem],
    sources: &mut Vec<(K, R::Elem)>,
    mut emit: impl FnMut(K, R::Elem),
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
}

/// The working space for the exponentials, kept from one curve to the next.
struct ExpWorkspace<K, E> {
    /// The terms of each degree, from zero to the maximum degree.
    levels: Vec<LevelTable<K, E>>,
    l: LTerms<K, E>,
    /// The nonzero terms of the current degree.
    sources: Vec<(K, E)>,
}

impl<K: Key, E: Copy + Default + PartialEq> ExpWorkspace<K, E> {
    fn new(sentinel: K, max_deg: u32) -> Self {
        Self {
            levels: (0..=max_deg).map(|_| LevelTable::new(sentinel)).collect(),
            l: LTerms::new(),
            sources: Vec::new(),
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
        self.l.build(ring, curve, alpha, ints, top);
        *self.levels[0].get_or_insert(K::default()) = ring.one();
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
}

impl<K: Key, E: Copy + Default + PartialEq> Shard<K, E> {
    /// Subtracts the terms of `out` and empties it.
    fn fold<R: Ring<Elem = E>>(&mut self, ring: &R, out: &mut OutBuf<K, E>) {
        let width = self.map.width();
        for (i, (&key, &deg)) in out.keys.iter().zip(&out.degs).enumerate() {
            let (t, created) = self.map.insert(key, deg);
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
}

/// Terms to subtract from one shard, as a worker collects them.
struct OutBuf<K, E> {
    keys: Vec<K>,
    degs: Vec<u32>,
    vals: Vec<E>,
}

impl<K, E> OutBuf<K, E> {
    fn new() -> Self {
        Self {
            keys: Vec::new(),
            degs: Vec::new(),
            vals: Vec::new(),
        }
    }

    fn clear(&mut self) {
        self.keys.clear();
        self.degs.clear();
        self.vals.clear();
    }
}

/// Buffers $n_C f$ to be subtracted at `key`, of degree `deg`, in the buffer of its shard in
/// `out`, and subtracts the buffer if it is full.
#[inline]
fn emit<R: Ring, K: Key>(
    ring: &R,
    shards: &[Mutex<Shard<K, R::Elem>>],
    out: &mut [OutBuf<K, R::Elem>],
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
    if o.keys.len() >= FLUSH {
        shards[s].lock().unwrap().fold(ring, o);
    }
}

/// What a worker thread keeps from one curve to the next.
struct Worker<K, E> {
    ws: ExpWorkspace<K, E>,
    curve: Vec<i32>,
    /// One buffer per shard.
    out: Vec<OutBuf<K, E>>,
}

/// A curve whose exponential is built together with those of the other curves of its layer,
/// one degree at a time, so that the scatters into the levels of every curve run in parallel.
struct SyncCurve<K, E> {
    key: K,
    deg: u32,
    /// The highest degree of its exponential.
    top: u32,
    l: LTerms<K, E>,
    /// The terms of each degree up to `top`, each locked by the one scatter that writes it.
    levels: Vec<Mutex<LevelTable<K, E>>>,
    sources: Vec<(K, E)>,
    /// One buffer per shard.
    out: Vec<OutBuf<K, E>>,
}

/// What the curves of a layer share while they are applied.
struct LayerCtx<'a, R: Ring, L: KeyLayout> {
    ring: &'a R,
    layout: &'a L,
    alpha: &'a SparsePoly<L::Key, R::Elem>,
    ints: &'a [R::Elem],
    inv: &'a [R::Elem],
    max_deg: u32,
    shards: &'a [Mutex<Shard<L::Key, R::Elem>>],
    /// The curves, each with its coefficients, `width` values each.
    curves: &'a [(L::Key, u32)],
    vals: &'a [R::Elem],
    width: usize,
}

impl<R: Ring, L: KeyLayout> LayerCtx<'_, R, L> {
    fn n_c(&self, i: usize) -> &[R::Elem] {
        &self.vals[i * self.width..(i + 1) * self.width]
    }

    /// Applies each curve on one thread, with the workers indexed by the threads.
    fn apply_by_curve(&self, workers: &[Mutex<Worker<L::Key, R::Elem>>]) {
        let apply = |i: usize| {
            let (kc, dc) = self.curves[i];
            let n_c = self.n_c(i);
            let t = rayon::current_thread_index().unwrap_or(0);
            let mut worker = workers[t].lock().unwrap();
            let Worker { ws, curve, out } = &mut *worker;
            self.layout.unpack(kc, curve);
            ws.exp(
                self.ring,
                curve,
                self.alpha,
                self.ints,
                self.inv,
                self.max_deg - dc,
                |m, e, f| emit(self.ring, self.shards, out, kc + m, dc + e, n_c, f),
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
    fn apply_by_level(&self) -> Vec<SyncCurve<L::Key, R::Elem>> {
        let ring = self.ring;
        let sentinel = self.layout.sentinel();
        let mut curves: Vec<SyncCurve<L::Key, R::Elem>> = self
            .curves
            .par_iter()
            .map(|&(key, deg)| {
                let top = self.max_deg - deg;
                let mut curve = vec![0; self.layout.dim()];
                self.layout.unpack(key, &mut curve);
                let mut l = LTerms::new();
                l.build(ring, &curve, self.alpha, self.ints, top);
                let levels: Vec<_> = (0..=top)
                    .map(|_| Mutex::new(LevelTable::new(sentinel)))
                    .collect();
                *levels[0].lock().unwrap().get_or_insert(L::Key::default()) = ring.one();
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
                let n_c = self.n_c(i);
                let (kc, dc) = (c.key, c.deg);
                let out = &mut c.out;
                finalize(ring, level, e, self.inv, &mut c.sources, |m, f| {
                    emit(ring, self.shards, out, kc + m, dc + e, n_c, f)
                });
                // This level is done.
                *level = LevelTable::new(sentinel);
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
            tasks.sort_unstable_by(|a, b| b.cmp(a));
            tasks.par_iter().for_each(|&(_, i, q)| {
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
/// $q^C$, at every curve class of degree at most `max_deg`, as a series with as many values
/// per term as `inst`.
///
/// `alpha` holds the $\alpha_a$, one value per coordinate of `layout`, and must have no
/// constant term; `inv` holds the inverses of the integers up to `max_deg` at least. Only the
/// curves in `lightcone` are subtracted, so the coefficients are only right in there.
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
#[allow(clippy::too_many_arguments)]
pub fn extract<R: Ring, L: KeyLayout>(
    ring: &R,
    layout: &L,
    alpha: &SparsePoly<L::Key, R::Elem>,
    inst: &SparsePoly<L::Key, R::Elem>,
    inv: &[R::Elem],
    max_deg: u32,
    lightcone: &Lightcone<L::Key>,
    tuning: &Tuning,
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
    let threads = rayon::current_num_threads();
    let n_shards = if threads == 1 { 1 } else { 4 * threads };

    let mut shards: Vec<Mutex<Shard<L::Key, R::Elem>>> = (0..n_shards)
        .map(|_| {
            Mutex::new(Shard {
                map: SlotMap::new(width, 2 * inst.len() / n_shards),
                heap: DegreeHeap::new(),
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
    let mut workers: Vec<Mutex<Worker<L::Key, R::Elem>>> = (0..threads)
        .map(|_| {
            Mutex::new(Worker {
                ws: ExpWorkspace::new(layout.sentinel(), max_deg),
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
                layer.push((kc, dc));
                layer_vals.extend_from_slice(n_c);
            }
        }
        if layer.is_empty() {
            continue;
        }

        let ctx = LayerCtx {
            ring,
            layout,
            alpha,
            ints: &ints,
            inv,
            max_deg,
            shards: &shards,
            curves: &layer,
            vals: &layer_vals,
            width,
        };
        let by_level = threads > 1
            && layer.len() < tuning.sync_max_curves.saturating_mul(threads)
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
        shards.par_iter_mut().zip(outs).for_each(|(shard, outs)| {
            let shard = shard.get_mut().unwrap();
            for o in outs {
                shard.fold(ring, o);
            }
        });
    }
    let parts = shards
        .into_iter()
        .map(|s| s.into_inner().unwrap().map.take())
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
            let got = extract(
                &lf,
                &layout,
                &alpha,
                &inst,
                &inv,
                max_deg,
                &Lightcone::All,
                &Tuning::default(),
            )
            .unwrap();
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
            extract(
                &lf,
                &layout,
                &empty,
                &inst,
                &inv,
                5,
                &Lightcone::All,
                &Tuning::default()
            ),
            Err(Error::DimensionMismatch)
        );
        assert_eq!(
            extract(
                &lf,
                &layout,
                &SparsePoly::new(3),
                &inst,
                &inv,
                4,
                &Lightcone::All,
                &Tuning::default()
            ),
            Err(Error::DimensionMismatch)
        );
        assert!(extract(
            &lf,
            &layout,
            &empty,
            &inst,
            &inv,
            4,
            &Lightcone::All,
            &Tuning::default()
        )
        .unwrap()
        .is_empty());
    }
}
