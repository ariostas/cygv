//! The extraction along a tree of the rows of the charges, which shares the exponentials of
//! the rows of the most terms between the curves (`PLAN.md` phase 8e, `ALGORITHMS.md` idea 1).
//!
//! Along the rows $M_r$, $C \cdot \alpha = \sum_r m_r(C)\, \beta_r$ with $m_r(C) = M_r \cdot C$,
//! so with $E_r = \exp \beta_r$ the exponential of a curve factors as
//! $\exp(C \cdot \alpha) = \prod_r E_r^{m_r(C)}$. Each $\beta_r$ has the terms of the
//! classes with no negative intersection and those whose one negative intersection is with
//! row $r$, so a few rows hold most of the terms. Take the $k$ rows of the most terms,
//! $r_1, \dots, r_k$, and let $L'_C$ be $C \cdot \alpha$ along the others. The series to
//! subtract, $\sum_C n_C\, z^C \exp(C \cdot \alpha)$, is then evaluated as a tree over the
//! prefixes $(m_{r_1}(C), \dots, m_{r_j}(C))$ of the curves: a node of depth $j$ holds
//! $N = \sum_C n_C\, z^C \prod_{i > j} E_{r_i}^{m_{r_i}(C)} \exp(L'_C)$ over its curves, which
//! is $\sum_m E_{r_{j+1}}^m N_m$ over its children $N_m$, and the root holds the series. Each
//! curve builds $\exp(L'_C)$ as [`extract`] builds its exponential, but adds $n_C\, z^C$ times
//! it to its leaf, and every node with $m \ne 0$ is multiplied by the power $E^m$ of its row
//! once, for all of its curves, instead of once per curve.
//!
//! The products are relaxed. Every term of a power but the first has positive degree, so the
//! terms of a node of degree $e$ are final once the curves of degree at most $e$ are applied,
//! and those of degree below $e$ scattered against the powers. The curves are applied a layer
//! of degrees $[d, d + \delta)$ at a time, as by [`extract`]. Before a layer, the terms of its
//! degrees are summed up the tree (the first term, 1, of each power) and subtracted from the
//! series at the root; the curves of the layer are then read, and their terms added to their
//! leaves, and $n_C\, z^C$ to every node above them but the root (the series keeps it); and
//! then the terms of the degrees of the layer, now final, are scattered against the other terms
//! of the power of each node into its parent. The nodes keep their terms in level tables of
//! absolute classes, keyed by linear hashes over a box that holds every class the extraction
//! touches.
//!
//! Each kind of work runs in parallel: the curves of a layer, as by [`extract`]; the sums by
//! parent; the subtractions by table; and the products by target table, each writing one table
//! of one node. Unless the ring is exact, every sum is taken in an order that does not depend
//! on the threads: the curves of a layer are sorted, their terms kept until they are all built
//! and then added in that order, and the subtractions buffered until the layer's sums are done.
//!
//! A series with several values per term (an n-fold) is extracted one component at a time,
//! since the components do not mix.

use super::*;
use std::collections::HashMap;

/// When [`Tuning::tree_depth`] leaves the depth open, the rows of the tree are those with at
/// least this fraction of the terms of the row of the most terms, at most [`AUTO_MAX_DEPTH`] of
/// them.
const AUTO_SHARE: f64 = 0.1;

/// The most rows of the tree when [`Tuning::tree_depth`] leaves the depth open. The nodes cost
/// more with each row, while the rows below the first two hold few terms: on the benchmark
/// models a third row saves pairs but no time, except on one thread at $h^{1,1} = 10$.
const AUTO_MAX_DEPTH: usize = 2;

/// The terms of positive degree of a power $E_r^m$, in order of degree, as hashes and values,
/// with the range of each degree.
struct Power<W, E> {
    hashes: Vec<W>,
    vals: Vec<E>,
    groups: Vec<(u32, usize, usize)>,
}

/// A power, by the depth of the node whose row it is of, less one, and its exponent.
type PowerKey = (usize, i64);

/// The series of a row of the tree, with the hashes of its terms.
type RowSeries<K, W, E> = (SparsePoly<K, E>, Vec<W>);

/// A worker thread, behind a lock that only its thread takes, but when the terms it kept are
/// added to the leaves.
type Worker<K, W, E> = Mutex<TreeWorker<K, W, E>>;

/// The final terms of a node of one degree: the node, the degree and the terms.
type NodeTerms<W, E> = (usize, u32, Vec<(W, E)>);

/// A product into a table: the target node and degree, the index of the terms of the source
/// among the [`NodeTerms`], and the range of the terms of the power.
type Task = ((usize, u32), usize, usize, usize);

/// A node of the tree: the curves with a prefix $(m_1, \dots, m_j)$ of the exponents of the
/// rows of the tree, with the terms known so far in a table per degree.
struct Node<W, E> {
    parent: usize,
    depth: usize,
    /// The last entry of the prefix, the exponent of the power of its row.
    m: i64,
    tabs: Vec<Mutex<LevelTable<W, E>>>,
}

/// A curve of a layer.
#[derive(Clone, Copy)]
struct TreeCurve<K, W, E> {
    key: K,
    deg: u32,
    /// What is subtracted at the curve: $n_C$, or what the source gives.
    val: E,
    leaf: usize,
    /// The hash of the class.
    hash: W,
}

/// What a worker thread keeps from one curve to the next.
struct TreeWorker<K, W, E> {
    ws: ExpWorkspace<W, E>,
    curve: Vec<i32>,
    /// The terms of the curves built so far, as hashes, degrees and values, when they are
    /// kept until the layer is done.
    terms: Vec<(W, u32, E)>,
    /// The index of each curve built so far, with the range of its terms.
    curves: Vec<(usize, usize, usize)>,
    /// The tables of the leaves that the terms filled, as node and degree.
    filled: Vec<(usize, u32)>,
    /// One buffer per shard.
    out: Vec<OutBuf<K, E>>,
}

/// The node of `prefix`, created with its ancestors if needed.
fn node_of<W: LevelWord, E: Copy + Default>(
    prefix: &[i64],
    nodes: &mut Vec<Node<W, E>>,
    index: &mut HashMap<Vec<i64>, usize>,
    marker: W,
    max_deg: u32,
) -> usize {
    if let Some(&n) = index.get(prefix) {
        return n;
    }
    let parent = node_of(&prefix[..prefix.len() - 1], nodes, index, marker, max_deg);
    nodes.push(Node {
        parent,
        depth: prefix.len(),
        m: prefix[prefix.len() - 1],
        tabs: new_tabs(marker, max_deg),
    });
    index.insert(prefix.to_vec(), nodes.len() - 1);
    nodes.len() - 1
}

/// The tables of a node, one per degree up to `max_deg`, empty.
fn new_tabs<W: LevelWord, E: Copy + Default>(
    marker: W,
    max_deg: u32,
) -> Vec<Mutex<LevelTable<W, E>>> {
    (0..=max_deg)
        .map(|_| Mutex::new(LevelTable::new(marker)))
        .collect()
}

/// Adds the terms of `sources` to `target`.
fn add_all<R: Ring, W: LevelWord>(
    ring: &R,
    sources: impl IntoIterator<Item = (W, R::Elem)>,
    target: &mut LevelTable<W, R::Elem>,
) {
    for (h, v) in sources {
        let t = target.get_or_insert(h);
        *t = ring.add(*t, v);
    }
}

/// Adds $f_v l_j$ at $v + u_j$ to `target` for every source $(v, f_v)$ and every term
/// $(u_j, l_j)$ of `hs` and `vs`, the classes given by their hashes.
fn product<R: Ring, W: LevelWord>(
    ring: &R,
    sources: &[(W, R::Elem)],
    hs: &[W],
    vs: &[R::Elem],
    target: &mut LevelTable<W, R::Elem>,
) {
    for &(hb, fb) in sources {
        for (j, (&hu, &lu)) in hs.iter().zip(vs).enumerate() {
            if let Some(&hp) = hs.get(j + PREFETCH) {
                target.prefetch(hb.wadd(hp));
            }
            ring.fma(target.get_or_insert(hb.wadd(hu)), lu, fb);
        }
    }
}

/// [`extract`] along the tree of the rows `rows`, with `beta` along them.
#[allow(clippy::too_many_arguments)]
pub(super) fn extract_tree<R: Ring, L: KeyLayout>(
    ring: &R,
    layout: &L,
    grading: &[i32],
    beta: &SparsePoly<L::Key, R::Elem>,
    rows: &[Vec<i32>],
    inst: &SparsePoly<L::Key, R::Elem>,
    inv: &[R::Elem],
    max_deg: u32,
    lightcone: &Lightcone<L::Key>,
    tuning: &Tuning,
    source: Option<Source<'_, R::Elem>>,
) -> Extracted<L::Key, R::Elem> {
    let h11 = layout.dim();
    // Every class the extraction touches is a sum of terms of beta and of the series.
    let mut coords = Vec::new();
    let mut degs = Vec::new();
    let mut c = vec![0; h11];
    for p in [beta, inst] {
        for i in 0..p.len_upto(max_deg) {
            if p.degs()[i] > 0 {
                layout.unpack(p.keys()[i], &mut c);
                coords.extend_from_slice(&c);
                degs.push(p.degs()[i]);
            }
        }
    }
    let terms = degs
        .iter()
        .enumerate()
        .map(|(i, &d)| (&coords[i * h11..(i + 1) * h11], d));
    let level_box = level_box(terms, grading, max_deg);
    drop((coords, degs));
    let run = |level: &dyn Fn(usize, &SparsePoly<_, _>) -> Extracted<_, _>| {
        let width = inst.width();
        if width == 1 {
            return level(0, inst);
        }
        // One component at a time.
        let mut map = SlotMap::new(width, inst.len());
        for k in 0..width {
            let mut part = SparsePoly::new(1);
            for i in 0..inst.len() {
                let v = inst.values(i)[k];
                if !ring.is_zero(v) {
                    part.push(inst.keys()[i], inst.degs()[i], &[v]);
                }
            }
            let got = level(k, &part)?;
            for i in 0..got.len() {
                map.get_or_insert(got.keys()[i], got.degs()[i])[k] = got.values(i)[0];
            }
        }
        Ok(map.take())
    };
    // The narrowest level word that fits, as in the extraction curve by curve.
    macro_rules! with {
        ($level:expr) => {
            run(&|k, part| {
                tree_in(
                    ring,
                    layout,
                    &$level,
                    beta,
                    rows,
                    part,
                    inv,
                    max_deg,
                    lightcone,
                    tuning,
                    source,
                    k,
                    inst.width(),
                )
            })
        };
    }
    type Wide<L> = <<L as KeyLayout>::Key as Key>::Level;
    let wide = <Wide<L> as LevelWord>::BITS > u128::BITS;
    if !tuning.wide_levels {
        if let Some(level) = LevelLayout::<u64>::new(&level_box, grading) {
            return with!(level);
        }
    }
    if !tuning.wide_levels || !wide {
        if let Some(level) = LevelLayout::<u128>::new(&level_box, grading) {
            return with!(level);
        }
    }
    if wide {
        if let Some(level) = LevelLayout::<Wide<L>>::new(&level_box, grading) {
            return with!(level);
        }
    }
    Err(Error::Key(KeyError::TooWide {
        bits: level_box.bits(),
        max_bits: <Wide<L> as LevelWord>::BITS - 1,
    }))
}

/// The rows of the tree for `beta` along `rows`, those of the most terms up to `max_deg` first:
/// `depth` of them, or as [`AUTO_SHARE`] says.
fn tree_rows<R: Ring, K: Key>(
    ring: &R,
    beta: &SparsePoly<K, R::Elem>,
    n_rows: usize,
    max_deg: u32,
    depth: Option<usize>,
) -> Vec<usize> {
    let mut counts = vec![0usize; n_rows];
    for u in 0..beta.len_upto(max_deg) {
        for (n, &v) in counts.iter_mut().zip(beta.values(u)) {
            if !ring.is_zero(v) {
                *n += 1;
            }
        }
    }
    let mut order: Vec<usize> = (0..n_rows).filter(|&r| counts[r] > 0).collect();
    order.sort_by_key(|&r| core::cmp::Reverse(counts[r]));
    let depth = depth.unwrap_or_else(|| {
        let top = order.first().map_or(0, |&r| counts[r]) as f64;
        order
            .iter()
            .take_while(|&&r| counts[r] as f64 >= AUTO_SHARE * top)
            .count()
            .min(AUTO_MAX_DEPTH)
    });
    order.truncate(depth);
    order
}

/// [`extract_tree`] for component `k` of `width`, given in `inst`, with the classes hashed by
/// `level`.
#[allow(clippy::too_many_arguments)]
fn tree_in<R: Ring, L: KeyLayout, W: LevelWord>(
    ring: &R,
    layout: &L,
    level: &LevelLayout<W>,
    beta: &SparsePoly<L::Key, R::Elem>,
    rows: &[Vec<i32>],
    inst: &SparsePoly<L::Key, R::Elem>,
    inv: &[R::Elem],
    max_deg: u32,
    lightcone: &Lightcone<L::Key>,
    tuning: &Tuning,
    source: Option<Source<'_, R::Elem>>,
    k: usize,
    width: usize,
) -> Extracted<L::Key, R::Elem> {
    let h11 = layout.dim();
    let marker = level.marker();
    let units: Vec<L::Key> = (0..h11)
        .map(|t| {
            let mut e = vec![0; h11];
            e[t] = 1;
            layout.try_pack(&e).unwrap_or_default()
        })
        .collect();
    let keys = level.offset_keys(&units);
    let ints: Vec<R::Elem> = (0..=max_deg).map(|n| ring.int(n as i64)).collect();
    let mut c = vec![0; h11];
    let n_beta = beta.len_upto(max_deg);
    let hashes: Vec<W> = (0..n_beta)
        .map(|u| {
            layout.unpack(beta.keys()[u], &mut c);
            level.hash(&c)
        })
        .collect();

    // The rows of the tree, the series of each one, and beta along the other rows, without the
    // terms where they all vanish.
    let tree = tree_rows(ring, beta, rows.len(), max_deg, tuning.tree_depth);
    let depth = tree.len();
    let row_series: Vec<RowSeries<L::Key, W, R::Elem>> = tree
        .iter()
        .map(|&r| {
            let mut b = SparsePoly::new(1);
            let mut h = Vec::new();
            for (u, &hash) in hashes.iter().enumerate() {
                let v = beta.values(u)[r];
                if !ring.is_zero(v) {
                    b.push(beta.keys()[u], beta.degs()[u], &[v]);
                    h.push(hash);
                }
            }
            (b, h)
        })
        .collect();
    let rest: Vec<usize> = (0..rows.len()).filter(|r| !tree.contains(r)).collect();
    let rest_rows: Vec<Vec<i32>> = rest.iter().map(|&r| rows[r].clone()).collect();
    let mut rest_beta = SparsePoly::new(rest.len());
    let mut rest_hashes = Vec::new();
    let mut v = vec![R::Elem::default(); rest.len()];
    for (u, &hash) in hashes.iter().enumerate() {
        for (x, &r) in v.iter_mut().zip(&rest) {
            *x = beta.values(u)[r];
        }
        if v.iter().any(|&x| !ring.is_zero(x)) {
            rest_beta.push(beta.keys()[u], beta.degs()[u], &v);
            rest_hashes.push(hash);
        }
    }

    // The power $E_j^m$ of the row of depth j + 1, up to degree `top`.
    let power = |j: usize, m: i64, top: u32| {
        let mut p = Power {
            hashes: Vec::new(),
            vals: Vec::new(),
            groups: Vec::new(),
        };
        let (b, h) = &row_series[j];
        let mut ws = ExpWorkspace::new(marker, top);
        ws.exp(ring, &[m as i32], b, h, None, &ints, inv, top, |h, e, f| {
            match p.groups.last_mut() {
                Some(g) if g.0 == e => g.2 += 1,
                _ => p.groups.push((e, p.hashes.len(), p.hashes.len() + 1)),
            }
            p.hashes.push(h);
            p.vals.push(f);
        });
        p
    };
    let mut powers: HashMap<PowerKey, Power<W, R::Elem>> = HashMap::new();

    let threads = rayon::current_num_threads();
    let n_shards = if threads == 1 { 1 } else { 4 * threads };
    let mut shards: Vec<Mutex<Shard<L::Key, R::Elem>>> = (0..n_shards)
        .map(|_| {
            Mutex::new(Shard {
                map: SlotMap::new(1, 2 * inst.len() / n_shards),
                heap: DegreeHeap::new(),
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
    let mut workers: Vec<Worker<L::Key, W, R::Elem>> = (0..threads)
        .map(|_| {
            Mutex::new(TreeWorker {
                ws: ExpWorkspace::new(marker, max_deg),
                curve: vec![0; h11],
                terms: Vec::new(),
                curves: Vec::new(),
                filled: Vec::new(),
                out: (0..n_shards).map(|_| OutBuf::new()).collect(),
            })
        })
        .collect();
    let mut nodes: Vec<Node<W, R::Elem>> = vec![Node {
        parent: usize::MAX,
        depth: 0,
        m: 0,
        tabs: new_tabs(marker, max_deg),
    }];
    let mut index: HashMap<Vec<i64>, usize> = HashMap::from([(Vec::new(), 0)]);
    // The nodes whose tables of each degree have terms, possibly with repetitions.
    let mut active: Vec<Vec<usize>> = vec![Vec::new(); max_deg as usize + 1];
    let one = ring.one();
    let in_lightcone = |key: L::Key, c: &mut [i32]| {
        layout.unpack(key, c);
        lightcone.contains(key, c)
    };
    let mut values = vec![R::Elem::default(); width];

    let delta = beta.degs().first().map_or(1, |&d| d.max(1));
    let mut layer: Vec<TreeCurve<L::Key, W, R::Elem>> = Vec::new();
    let mut d = 1;
    while d <= max_deg {
        let top = (d + delta).min(max_deg + 1);

        // Sum the terms of the degrees of the layer up the tree, and subtract the root and its
        // children from the series.
        for x in d..top {
            let xi = x as usize;
            for j in (2..=depth).rev() {
                let mut kids: Vec<(usize, usize)> = active[xi]
                    .iter()
                    .filter(|&&n| nodes[n].depth == j)
                    .map(|&n| (nodes[n].parent, n))
                    .collect();
                kids.sort_unstable();
                kids.dedup();
                let nodes = &nodes;
                let filled: Vec<usize> = kids
                    .chunk_by(|a, b| a.0 == b.0)
                    .collect::<Vec<_>>()
                    .into_par_iter()
                    .filter_map(|g| {
                        let p = g[0].0;
                        let mut t = nodes[p].tabs[xi].lock().unwrap();
                        let was_empty = t.is_empty();
                        for &(_, n) in g {
                            add_all(ring, nodes[n].tabs[xi].lock().unwrap().iter(), &mut t);
                        }
                        (was_empty && !t.is_empty()).then_some(p)
                    })
                    .collect();
                active[xi].extend(filled);
            }
            let mut emitted: Vec<usize> = active[xi]
                .iter()
                .copied()
                .filter(|&n| nodes[n].depth == 1)
                .chain([0])
                .collect();
            emitted.sort_unstable();
            emitted.dedup();
            let (nodes, shards, workers, keys) = (&nodes, &shards, &workers, &keys);
            emitted.par_iter().enumerate().for_each(|(i, &n)| {
                let t = nodes[n].tabs[xi].lock().unwrap();
                let w = rayon::current_thread_index().unwrap_or(0);
                let mut worker = workers[w].lock().unwrap();
                let TreeWorker { curve, out, .. } = &mut *worker;
                let src = (!R::EXACT).then_some(i as u32);
                for (h, f) in t.iter() {
                    if ring.is_zero(f) {
                        continue;
                    }
                    let key = keys.key(h, x);
                    if in_lightcone(key, curve) {
                        emit(ring, shards, out, src, key, x, &[one], f);
                    }
                }
            });
        }
        for x in d..top {
            nodes[0].tabs[x as usize].get_mut().unwrap().release();
        }
        let mut outs: Vec<Vec<&mut OutBuf<L::Key, R::Elem>>> =
            (0..n_shards).map(|_| Vec::new()).collect();
        for worker in workers.iter_mut() {
            for (s, o) in worker.get_mut().unwrap().out.iter_mut().enumerate() {
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

        // The curves of the layer, with their leaves, sorted by leaf and then key.
        layer.clear();
        for shard in shards.iter_mut() {
            let shard = shard.get_mut().unwrap();
            while shard.heap.min_deg().is_some_and(|e| e < top) {
                let slot = shard.heap.pop().unwrap();
                let (key, deg) = (shard.map.key(slot), shard.map.deg(slot));
                let mut val = shard.map.values(slot)[0];
                // A curve of the maximum degree has nothing left to change, and one outside
                // the lightcones only changes curves outside them.
                if deg == 0 || deg >= max_deg || ring.is_zero(val) || !in_lightcone(key, &mut c) {
                    continue;
                }
                if let Some(source) = source {
                    source(&c, deg, &mut values);
                    val = values[k];
                    if ring.is_zero(val) {
                        continue;
                    }
                }
                layer.push(TreeCurve {
                    key,
                    deg,
                    val,
                    leaf: 0,
                    hash: level.hash(&c),
                });
            }
        }
        // The nodes are numbered as they appear, in an order that does not depend on the
        // threads, since the sums go in their order.
        layer.sort_unstable_by_key(|t| t.key);
        for t in layer.iter_mut() {
            layout.unpack(t.key, &mut c);
            let prefix: Vec<i64> = tree
                .iter()
                .map(|&r| super::super::dot(&c, &rows[r]))
                .collect();
            t.leaf = node_of(&prefix, &mut nodes, &mut index, marker, max_deg);
        }

        // Build the exponentials of the curves along the other rows. An exact ring adds their
        // terms to the leaves right away; otherwise they are kept, and added in order.
        let apply = |i: usize| {
            let t = &layer[i];
            let w = rayon::current_thread_index().unwrap_or(0);
            let mut worker = workers[w].lock().unwrap();
            let TreeWorker {
                ws,
                curve,
                terms,
                curves,
                filled,
                ..
            } = &mut *worker;
            layout.unpack(t.key, curve);
            let start = terms.len();
            ws.exp(
                ring,
                curve,
                &rest_beta,
                &rest_hashes,
                Some(&rest_rows),
                &ints,
                inv,
                max_deg - t.deg,
                |h, e, f| terms.push((t.hash.wadd(h), t.deg + e, ring.mul(t.val, f))),
            );
            if R::EXACT {
                add_terms(ring, &nodes[t.leaf], t.leaf, &terms[start..], filled);
                terms.truncate(start);
            } else {
                curves.push((i, start, terms.len()));
            }
        };
        if threads == 1 {
            (0..layer.len()).for_each(apply);
        } else {
            (0..layer.len()).into_par_iter().for_each(apply);
        }
        if !R::EXACT {
            // Each leaf takes the terms of its curves in their order.
            let mut kept: Vec<(usize, usize, usize, usize)> = Vec::new();
            for (w, worker) in workers.iter_mut().enumerate() {
                let worker = worker.get_mut().unwrap();
                kept.extend(worker.curves.drain(..).map(|(i, s, e)| (i, w, s, e)));
            }
            kept.sort_unstable_by_key(|&(i, ..)| (layer[i].leaf, i));
            let (nodes, workers, layer) = (&nodes, &workers, &layer);
            let filled: Vec<(usize, u32)> = kept
                .chunk_by(|a, b| layer[a.0].leaf == layer[b.0].leaf)
                .collect::<Vec<_>>()
                .into_par_iter()
                .flat_map_iter(|g| {
                    let leaf = layer[g[0].0].leaf;
                    let mut filled = Vec::new();
                    for &(_, w, s, e) in g {
                        // Only this task reads the worker's terms now.
                        let worker = workers[w].lock().unwrap();
                        add_terms(ring, &nodes[leaf], leaf, &worker.terms[s..e], &mut filled);
                    }
                    filled
                })
                .collect();
            for (n, x) in filled {
                active[x as usize].push(n);
            }
        }
        for worker in workers.iter_mut() {
            let worker = worker.get_mut().unwrap();
            worker.terms.clear();
            for (n, x) in worker.filled.drain(..) {
                active[x as usize].push(n);
            }
        }
        // The curves themselves, in every node above them but the root.
        for t in &layer {
            let mut n = t.leaf;
            while n != 0 {
                let tab = nodes[n].tabs[t.deg as usize].get_mut().unwrap();
                if tab.is_empty() {
                    active[t.deg as usize].push(n);
                }
                add_all(ring, [(t.hash, t.val)], tab);
                n = nodes[n].parent;
            }
        }

        // The terms of the degrees of the layer are final: take them out of the nodes, but
        // those outside the lightcones, which only change classes outside them.
        let mut srcs: Vec<NodeTerms<W, R::Elem>> = Vec::new();
        for x in d..top {
            let mut list = core::mem::take(&mut active[x as usize]);
            list.retain(|&n| n != 0);
            list.sort_unstable();
            list.dedup();
            let (nodes, keys) = (&nodes, &keys);
            srcs.par_extend(list.into_par_iter().filter_map(|n| {
                let mut t = nodes[n].tabs[x as usize].lock().unwrap();
                let mut c = vec![0; h11];
                let s: Vec<(W, R::Elem)> = t
                    .iter()
                    .filter(|&(h, f)| {
                        !ring.is_zero(f)
                            && (lightcone.is_all() || in_lightcone(keys.key(h, x), &mut c))
                    })
                    .collect();
                t.release();
                (nodes[n].m != 0 && !s.is_empty()).then_some((n, x, s))
            }));
        }
        // The powers they need, as far as the lowest degree of their terms needs them.
        let mut need: Vec<(usize, i64, u32)> = srcs
            .iter()
            .map(|&(n, x, _)| (nodes[n].depth - 1, nodes[n].m, max_deg - x))
            .filter(|&(j, m, _)| !powers.contains_key(&(j, m)))
            .collect();
        need.sort_unstable_by_key(|&(j, m, top)| (j, m, core::cmp::Reverse(top)));
        need.dedup_by_key(|&mut (j, m, _)| (j, m));
        let computed: Vec<(PowerKey, Power<W, R::Elem>)> = need
            .into_par_iter()
            .map(|(j, m, top)| ((j, m), power(j, m, top)))
            .collect();
        powers.extend(computed);
        // The products into the parents, as tasks by target table, the largest first.
        let mut tasks: Vec<Task> = Vec::new();
        for (si, (n, x, _)) in srcs.iter().enumerate() {
            let pw = &powers[&(nodes[*n].depth - 1, nodes[*n].m)];
            for &(g, j0, j1) in &pw.groups {
                if g > max_deg - x {
                    break;
                }
                tasks.push(((nodes[*n].parent, x + g), si, j0, j1));
            }
        }
        tasks.sort_by_key(|t| t.0);
        let mut by_target: Vec<(usize, &[Task])> = tasks
            .chunk_by(|a, b| a.0 == b.0)
            .map(|g| {
                let cost = g
                    .iter()
                    .map(|&(_, si, j0, j1)| srcs[si].2.len() * (j1 - j0))
                    .sum();
                (cost, g)
            })
            .collect();
        by_target.sort_by_key(|&(cost, _)| core::cmp::Reverse(cost));
        let filled: Vec<(usize, u32)> = {
            let (nodes, powers, srcs) = (&nodes, &powers, &srcs);
            by_target
                .into_par_iter()
                .filter_map(|(_, g)| {
                    let (p, y) = g[0].0;
                    let mut t = nodes[p].tabs[y as usize].lock().unwrap();
                    let was_empty = t.is_empty();
                    for &(_, si, j0, j1) in g {
                        let (n, _, s) = &srcs[si];
                        let pw = &powers[&(nodes[*n].depth - 1, nodes[*n].m)];
                        product(ring, s, &pw.hashes[j0..j1], &pw.vals[j0..j1], &mut t);
                    }
                    (was_empty && !t.is_empty()).then_some((p, y))
                })
                .collect()
        };
        for (n, y) in filled {
            active[y as usize].push(n);
        }
        d = top;
    }
    let parts = shards
        .into_iter()
        .map(|s| s.into_inner().unwrap().map.take())
        .collect();
    Ok(SparsePoly::merge(1, parts))
}

/// Adds `terms`, as hashes, degrees and values in order of degree, to the tables of `node`,
/// whose index is `n`, and records in `filled` the tables that had no terms before.
fn add_terms<R: Ring, W: LevelWord>(
    ring: &R,
    node: &Node<W, R::Elem>,
    n: usize,
    terms: &[(W, u32, R::Elem)],
    filled: &mut Vec<(usize, u32)>,
) {
    for run in terms.chunk_by(|a, b| a.1 == b.1) {
        let x = run[0].1;
        let mut t = node.tabs[x as usize].lock().unwrap();
        if t.is_empty() {
            filled.push((n, x));
        }
        add_all(ring, run.iter().map(|&(h, _, f)| (h, f)), &mut t);
    }
}
