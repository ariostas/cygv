# Big algorithmic ideas for `cygv` — research and design report

Written on branch `gv-modular-rust` on 2026-10-13 (after `cbc3878`), by an agent asked for ideas
on the scale of the move from big rationals to word-sized primes and the CRT; `PLAN.md` section
5.3 tracks the ideas taken up. The measurements and prototypes it cites were scratch material and
are not in the repository:

- `measure/`: a Rust program, depending on the crate by path, that ran the modular stages on a
  bench geometry and printed sizes, cost models and histograms (outputs `h11_8_d20_*.txt` and
  `h11_11_d18_*.txt`).
- `tree_extraction_poc.py`: an exact-rational proof of concept of the top idea (shared-exponential
  extraction), checked against the per-curve extraction and against planted invariants.
- `lagrange_good_check.py`: an independent check of the multivariate Lagrange–Good formula.
- A `perf record` of the one-thread degree-20 run of `h11_8`.
- Three literature searches, summarized in section 6.

Labels: **[measured]** = measured on this machine with the code in `measure/`; **[verified]** =
checked by an exact computation; **[literature]** = established, with a source; **[estimate]** = cost
model or inference; **[speculation]** = plausible, unchecked; **[implemented]** = in the Rust code
now (see the outcome below).

**Outcome.** Ideas 1 and 2 were implemented in phase 8e of `PLAN.md` (`modular/extraction/tree.rs`),
with outputs identical to the earlier extraction: 3–6× end to end on the bench models at degree
20–22 and 13× on `h11_10` at 20, on the 4-core development machine. Two predictions below were off:
computing the powers E_r^m in full cost more than everything else until each one was computed
only up to the degree its node needs, and a product costs 35–55 ns a pair (tables near the root out
of cache) against the 14 ns of the scatter, so trees deeper than two rows save pairs but not time.
`PLAN.md` 8e has the measurements, and 5.3 the follow-ups. Idea 3 was then prototyped in the pass: the
formula holds, modulo the primes, at all 11 targets tried, but it needs an exponential for every
class it returns, so it cannot replace the extraction for full runs or for the target-point mode
as it is (every invariant in the lightcones). For the targets alone it was 1–2× faster than the
extraction on their lightcones and about 2.5× smaller, most of its time going to Δ and IΔ, which
the targets share; not adopted (`PLAN.md` 5.3).

---

## 0. Summary

1. **There is no known formulation that avoids the mirror-map composition.** GKZ/Picard–Fuchs,
   Givental's I/J functions and Birkhoff factorization, the Yukawa-coupling route, toric residues,
   quasimaps: every route ends with "re-expand a z-series in q = z e^α", which is exactly the
   extraction. No fast multivariate degree-truncated composition exists (fast composition is
   univariate only: Brent–Kung, Kinoshita–Li 2024). Dense/FFT methods lose to the current sparse
   pairwise arithmetic at h¹¹ ≈ 8–11, D ≈ 20 by 1–2 orders of magnitude because the cone truncation
   already makes a sparse product cost ≈ M²/C(2r, r) instead of M². Intersection numbers are not an
   avoidable input; the pipeline uses only the single bilinear form K_w = Σ_t w_t κ_t·· (threefolds),
   which `phase.rs` already derives from the fan. **[literature + estimate]**

2. **The extraction is 92% of the run and its structure is: one exponential per nonzero curve, each
   exponential a product of per-divisor "local mirror maps".** Because the one-negative period
   coefficients are c₁ = ±c₀-factor · Q_s, one has α = Qᵀβ with β_r scalar series supported on the
   small cones K_r, so exp(C·α) = Π_r E_r^{Q_r·C} with E_r = exp β_r fixed series. The current Euler
   recurrence already exploits this implicitly (L = C·α has only the terms of the divisors C meets),
   which is why its pair count is 30× below a naive model. **[measured]**

3. **The big lever is sharing those products across curves.** Curves sharing the exponents of the
   dominant divisors share the expensive factors; the subtraction Σ_C n_C z^C Π_r E_r^{m_r(C)} can be
   evaluated as a Horner tree over divisor exponents with *relaxed* (degree-online) products, which
   respects the degree-order constraint of the extraction. Exact pair counts on the real data
   **[measured]**: sharing only the two dominant divisors cuts the scatter work **14.4×** on `h11_8`
   at degree 20 (5.03e8 → 3.50e7 pairs) and **6.7×** on `h11_11` at degree 18 (7.1e7 → 1.06e7), and
   the terms folded into the series **121× / 29×**. Half of the remaining cost is per-curve
   exponentials over the other divisors, which deeper sharing removes; projected **10–30× on the
   extraction** at these sizes, growing with degree. Correctness of the online scheme is
   **[verified]** in `tree_extraction_poc.py`. This is idea 1 below, now implemented.

4. **For target-point / lightcone mode** two cheaper wins: prune the exponentials to the lightcone
   (**[measured]** 3.9× pairs, 50× emitted terms on `h11_8` with 3 targets; idea 2, now
   implemented), and the Lagrange–Good per-target formula (**[verified]**), which makes targets
   independent (parallel, memory-local).

5. Ideas with orders of magnitude *only in other regimes*: NTT products via randomized cyclic
   projections at D ≳ 30 (memory makes that regime rare), lightcone-chunked extraction as a
   time–memory tradeoff for the memory-bound degree-26–28 regime, symmetry orbits.

6. Rejected: Kronecker/box FFT, Taylor/θ expansion, triangular/Horner composition, Newton inversion
   of the mirror map, evaluation–interpolation over F_p, p-adic lifting, GKZ/PF/Yukawa/Birkhoff,
   WDVV, "track only candidate points" (nonzero GV occur at 169k non-candidate classes), topological
   vertex/localization, modular bootstraps (special fibrations only).

---

## 1. Where the time goes today **[measured]**

`h11_8` (h¹¹ = 8, 12 divisors), `max_deg` 20, one thread, release CLI with the bench cone data:
23.4 s (3 lanes; includes the lane probe at degree 15). `perf` (self time):

| function | share |
|---|---|
| `LTerms::scatter` (Euler recurrence: hash-add target, 3-lane fma) | 48.8% |
| `apply_by_curve` closure = `finalize` + `emit` inlined (decode hash → key, buffer) | 17.0% |
| `Shard::fold` + `SlotMap::insert_hashed` (subtract buffers into the series) | 13.9% |
| `emit` | 6.0% |
| `LTerms::build` | 4.2% |
| `LevelTable::insert` | 1.9% |
| probe pass (`Lanes<2>` symbols) | ≈ 1.5% |
| period, instanton, multicover, CRT, YAML | < 3% |

So the extraction is ≈ 92%, split ≈ 55% exponential construction and ≈ 37% emitting/folding the
exponentials' terms into the series.

Sizes at degree 20 (`measure`):

| quantity | `h11_8` @20 | `h11_11` @18 |
|---|---|---|
| lattice points of the Mori cone up to `max_deg` | 14.4 M | 94.9 M |
| … by number of negative intersections 0..7 | 25 / 21 k / 620 k / 3.46 M / 5.99 M / 3.70 M / 589 k / 15 k | 2 / 11 k / 507 k / 5.6 M / 22.8 M / 36.1 M / 23.3 M / 6.5 M |
| candidates (≤ 2 negatives) | 55 767 | 83 108 |
| K_∅ points (c₀ support) | 25 | 2 |
| α terms (≤ 1 negative) | 21 289 | 10 986 |
| semigroup generated by supp α (support of any exponential) | 2.48 M | 3.89 M |
| terms of I | 55 283 | 81 710 |
| curves with nonzero n_C (exponentials computed) | 225 207 | 113 661 |
| … of which **not** candidates (composite classes with nonzero GV) | 169 440 | 30 553 |
| union of all emitted supports (points ever touched, layers ≤ 14/13) | 366 k | 323 k |
| Σ_C \|exp(C·α)\| (terms emitted) | 4.44e7 | 9.42e6 |
| Euler scatter pairs Σ_C \|S_C\| ⊗ \|L_C\| | 5.6e8 (5.03e8 in layers 1–14) | 7.1e7 (PLAN 8d: 7.23e7) |
| extraction time (one thread, 2 lanes) | 13.4–16 s | 2.7–3.0 s |

Per pair the scatter costs ≈ 14 ns ≈ 37 cycles with 3 lanes: instruction-bound (IPC 1.7, PLAN 8d),
about 20 instructions per lane-fma including the table probe. The hash tables of a curve's
exponential fit in L2, so this is not memory latency.

**Cost by layer** (`h11_8_d20_layers_full.txt`, exact over all curves): layers 6–14 carry 83% of the
pairs, layers 1–5 16%, layers ≥ 15 < 2%. The expensive curves are the *hundreds to thousands* of
mid-degree curves, not the ten degree-1 curves (whose exponentials are the largest but few).

**Why the naive model overcounts 30×.** The period formulas give, at a point u with exactly one
negative intersection (with divisor s), c₁(u) = ±c₀-factor × Q_s, and α = c₀⁻¹c₁ keeps that
direction because c₀⁻¹ is supported on K_∅ (sums stay one-negative at s). Hence α(u) ∝ Q_s and
L_u = C·α(u) ∝ Q_s·C vanishes unless C meets divisor s. So:

- supp L_C = {u ∈ K_∅} ∪ ⋃_{s : Q_s·C ≠ 0} {one-negative points at s}; e.g. a degree-1 curve of
  `h11_8` with 4–5 nonzero intersections has |L_C| = 7457 of α's 21 289 terms;
- supp exp(C·α) ⊂ semigroup generated by K_∅ ∪ {K_s : Q_s·C ≠ 0} ⊂ K_{T(C)} (points negative
  only on the divisors C meets): 25.6 k points at top degree 19 instead of 2.48 M.

The α terms per divisor are very unequal: `h11_8` [0, 0, 42, 102, 45, 102, 286, 153, 644, 4928,
974, **13989**] + 24 zero-negative; `h11_11` [0, 0, 4, 54, 48, 8, 4, 31, 667, 291, 34, 2749, 529,
1176, **5390**, 1]. One divisor carries 50–66% of α, the top two 80%+. **[measured]**

This is the structure every idea below builds on.

---

## 2. Ranked ideas

### Idea 1 — Shared-exponential extraction: Horner tree over divisor exponents with relaxed products **[implemented]**

**Status: implemented in phase 8e of `PLAN.md` (`modular/extraction/tree.rs`), the default for
threefolds, with at most two rows shared (see the outcome above). Before that: strongest candidate;
cost model measured on real data; algorithm verified exactly on synthetic data.**

**Mathematics.** Write the extraction identity with α along the rows of the charge matrix (the
anticanonical rows j and the divisor rows r; the code already has this as `c1_rows`/`beta`/`rows`
for the certifier): α_a = Σ_r Q_{ra} β_r, so C·α = Σ_r m_r(C) β_r with m_r(C) = Q_r·C ∈ ℤ and

    exp(C·α) = Π_r E_r^{m_r(C)},     E_r := exp β_r.

Each β_r is a *scalar* series supported on K_r (one-negative at r, plus K_∅) — it is the local
mirror map of the surface D_r ∩ X in GKZ homogeneous coordinates — and E_r lives on the same small
cone (13 994 points for the biggest divisor of `h11_8` at degree 20; the next 4 929; most under
300). The series the extraction must subtract is

    T(G) = Σ_C n_C z^C Π_r E_r^{m_r(C)}.

Order the rows r₁, r₂, …, r_k (largest |K_r| first). A node of depth j is a prefix
(m_{r₁}, …, m_{r_j}) and holds

    N_prefix = Σ_{C with that prefix} n_C z^C Π_{i>j} E_{r_i}^{m_{r_i}(C)},
    N_prefix = Σ_m E_{r_{j+1}}^m · N_{prefix, m}       (distributivity),

the root is T(G) and the leaves (full exponent vector = one curve, since Q has rank h¹¹) are
n_C z^C. The expensive factors (the dominant divisors) are at the top, where there are only
~#distinct values of Q_{r₁}·C ≈ 40–60 nodes, and everything below a node is shared by all curves
with that prefix.

**The degree-order constraint.** n_C is only known after every curve of lower degree has been
subtracted, and T(G) is needed only degree by degree. Since every E^m = 1 + (positive degree), the
products can be *relaxed* (online in degree): processing degree e,

1. the final degree-e terms of a node = its pending degree-e terms (scatters received from the
   children's lower-degree finals) + its children's degree-e finals; the degree-e leaves are
   still unknown, and they would pass up the chain unchanged (times the constant 1 of each E^m),
   so the root's degree-e finals *without them* are exactly the correction to subtract from I at
   degree e;
2. read n_C = I[C] − correction[C] for every class C of degree e in the union of both supports;
3. insert the new leaves (n_C z^C added to the degree-e finals of their k ancestors — k additions
   per curve);
4. scatter every node's final degree-e terms against the positive-degree terms of its E^m into the
   parent's pending tables of degrees e + g (the current `LTerms::scatter` loop, with the node's
   terms as sources and E^m's terms as the "L" group), then free degree e.

This is the current level-table machinery (hashes of offsets add; `OffsetKeys` decodes only what is
emitted) with two changes: the sources of a scatter are a node's terms instead of one curve's
exponential, and the L-terms are the fixed series E^m instead of C·α. The per-curve Euler recurrence
disappears entirely if the tree goes down to the smallest rows; a hybrid (tree over the top j rows,
per-curve exponentials over the rest attached as leaves) is simpler and already captures most of
the gain.

**Measured cost model** (exact pair counts over *all* curves with nonzero n_C in the cost-carrying
layers; `h11_8_d20_layers_full.txt`, `h11_11_d18_layers_full.txt`; pairs = Σ_e |source level e| ×
|E-terms of degree ≤ top − e|, the quantity the scatter loop executes):

| | `h11_8` @20, layers 1–14 | `h11_11` @18, layers 1–13 |
|---|---|---|
| Euler (today) | 5.03e8 | 7.10e7 |
| hybrid, share divisor r₁ only, products per layer | 1.46e8 (3.4×) | 3.19e7 (2.2×) |
| hybrid, share r₁, r₂, products per layer | 7.46e7 (6.7×) | 1.81e7 (3.9×) |
| hybrid, share r₁, r₂, **relaxed** (cross-layer) products | **3.50e7 (14.4×)** | **1.06e7 (6.7×)** |
| … of which per-curve exponentials over the other rows | 1.76e7 | 5.6e6 |
| … products at level 2 (353–808 nodes) | 6.8e6 | 2.0e6 |
| … products at level 1 (39–53 nodes) | 1.06e7 | 3.0e6 |
| terms emitted/folded into I: today Σ_C \|exp\| → union | 4.44e7 → 3.66e5 (121×) | 9.42e6 → 3.23e5 (29×) |

A plain product tree per curve (balanced split of the divisors, no sharing) already gives 3.8× /
2.8× **[measured]**; the sharing is what gets the rest: among the 6 314 curves of degree ≤ 12 of
`h11_8` there are 187 691 (curve, divisor-pair) incidences but only 7 948 distinct
(s, t, Q_s·C, Q_t·C) and 320 distinct (s, Q_s·C) **[measured]**.

Projection **[estimate]**: with the tree carried down to all rows, the remaining per-curve part
(half of the 2-level total) becomes deep nodes with tiny supports and tiny E's; total scatter work
≈ 2e7 for `h11_8`@20, i.e. ≈ 25× below today, and the fold work ≈ 100× below. Since the extraction
is 92% of the run, that is roughly a 10–20× end-to-end speed-up at these sizes if the per-pair
constant holds (same loop body). The ratio grows with degree: the number of top-level nodes grows
like D while the number of curves grows like D^{h¹¹−1}.

**Why it beats the alternatives.** The composition literature (section 6) reorganizes T⁻¹ as dense
operations on the whole cone and loses the small supports; this scheme keeps every product on the
cones K_T the curves actually live on, and shares exactly the factors that are shared. It subsumes
the C's curve-class replay (PLAN 5.3, ~10%).

**Requirements.** β along rows (exists: `period(…, rows=true)`, `instanton::beta`); power tables
E_r^m for the m that occur (≤ ~120 values per row; each one is a one-divisor exponential, cheap:
< 3e5 pairs each **[measured]** via the leaf costs); a node structure with per-degree level tables
keyed by the offset hash (the tables of a node hold absolute points, so the hash box must cover the
touched set rather than one exponential — the `u128` level words already exist for that); a
bottom-up finalize/merge per degree. The grading, keys, CRT, multicover, lightcones and n-fold
widths are untouched.

**What it covers / breaks.** Threefolds and n-folds (leaves carry M values), GV and GW, complete
intersections (anticanonical rows are just more rows), vex phases (negative m_0 is a negative power;
stratum curves are leaves like any other), target points (prune nodes to the lightcone, idea 2).
Certifier: the sums are the same terms reassociated, so bounds stay valid with `BoundRing`; the
"source" override applies at the leaves; reproducibility needs a fixed merge order of children
(natural). Exactness: unchanged (same ring operations).

**Risks.** (i) Per-pair constant: node tables are larger than a curve's level tables (level-1 nodes
hold 743 k points in all vs 366 k for the whole series today), so more cache misses; the current fold
into the shards has the same issue and costs ≈ 50 ns per term, so even a 3× worse constant leaves
a large gain. (ii) Memory: node tables ≈ 2–3× today's I table at `h11_8`@20; per-thread level tables
go away. (iii) Parallelism: nodes of a level are independent (39–808 of them); the root's children
scatter into the root through shards as today; at 64 threads the few level-1 products need splitting
by source range. (iv) Choice of row order and tree depth per geometry — a cost model from the
candidate counts per K_r (available before the extraction) can choose it.

**First experiment (1–2 days).** In a scratch crate against the current modules: implement the
1-level hybrid (share only r₁ = argmax |K_r|): per layer, group curves by m = Q_{r₁}·C, build each
curve's exponential of C·α − m β_{r₁} with the existing `ExpWorkspace`, sum them into one series
A_m per group (absolute keys), and scatter A_m against E_{r₁}^m into the shards — with cross-layer
accumulation (keep the A_m across layers and scatter only their newly finalized degree). Check the
GVs against `compute_gv` on `h11_8`@20 and `h11_11`@18 and compare the scatter counters with the
table above (expect ≈ 3× pairs, ≈ 10× emitted terms). Then the 2-level version.

### Idea 2 — Lightcone pruning inside the exponentials (target-point mode) **[implemented]**

**Status: implemented in phase 8e of `PLAN.md`, in the tree of idea 1 (its nodes' sources and the
emitted terms; the curve-by-curve extraction does not prune): 3.6× fewer product pairs on `h11_8`
with three degree-16 targets. Before that: measured; trivial to implement; large gain only with
targets.**

With targets P, a level-table entry at offset v of curve C is useful only if C + v is in the
lightcone (if it is not, no C + v + u is: the lightcone is a lower set). Today the pruning is at
three points outside the exponentials (PLAN 6.3 left this "unmeasured"). Test membership per
*final* level entry (hyperplane test, `Lightcone::contains`), not per scatter pair, and drop the
non-members from the sources and the emits.

**Measured** (`h11_8`, 3 targets of degree 16, 1 422 nonzero curves in their lightcones): scatter
pairs 3.2e7 → 8.2e6 (**3.9×**), emitted terms 2.65e6 → 5.25e4 (**50×**). With emit/fold ≈ 40% of the
time, expect ≈ 5–8× on the extraction in this mode; more for fewer or lower targets. Also the enabler
of idea 5. No effect on full runs. Works for every variant (the lightcone is a property of the keys).

### Idea 3 — Lagrange–Good per-target extraction **[measured, not adopted]**

**Status: prototyped in the Rust pass and measured (see the outcome above and `PLAN.md` 5.3): it
matches the extraction at every target tried, but only pays for a mode that returns the targets
alone, by 1–2× in time and about 2.5× in memory in an untuned prototype. Before that: formula
established [literature: Good 1960; Gessel 1987; Bender–Richmond 1998]; verified here exactly; cost
estimated, not measured.**

For q_a = z_a e^{α_a(z)} and any series F,

    [q^C] F(z(q)) = [z^C] F(z) · e^{−C·α(z)} · Δ(z),   Δ = det(δ_ab + θ_b α_a),  θ_b = z_b ∂_b.

With F = I this gives n_C = [z^C] (I Δ) e^{−C·α} for each target *independently*: one exponential
per target over its lower set (lightcone), no sequential subtraction, no dependence on other
curves' invariants. `lagrange_good_check.py` confirms it exactly (n = 3, D = 7, 119 coefficients);
the subagents confirmed it on the quintic and P(1,1,2,2,2)[8] with real period data **[verified]**.

Cost per target ≈ |S_P ∩ lightcone| ⊗ |L_P| — about the same total as idea 2 for the curves below P
**[estimate]** — plus Δ once: a k × k determinant of series, k = rank of the face the targets span
(Sylvester reduces it to the face), O(k³) series products on the lightcone; cheap for 2-faces, heavy
for faces of rank 32–40. Gains: perfect parallelism and distribution (one target per machine),
memory ∝ lightcone, and a natural fit for "potent rays" and 2-faces at large h¹¹. Not a win for full
truncations: every class would need an exponential, including the ~90% with zero invariant
(DKMMR: 360 M classes examined for 2 500 nonzero at h¹¹ = 101). First experiment: Python on
`h11_9` with the C's lightcone test cases (PLAN 6b), compare with idea 2.

### Idea 4 — Quasi-linear products by NTT with randomized cyclic projections (high-degree regime)

**Status: literature + estimate; crossover too high for today's typical runs.**

Map z^m ↦ x^{(e·m) mod N} τ^{deg m} with N ≈ 1.2–2× the slice size and 2–3 independent random
linear forms e; multiply by cyclic NTT in x and univariate truncation in τ; peel collisions using the
known output support (van der Hoeven–Lecerf 2013 "sparse polynomial and series multiplication";
Arnold–Roche 2015; Giorgi–Grenet–Perret du Cray 2020). Cost ≈ 2 000 operations per output coefficient
at n = 9, D = 20 **[estimate, subagent]**, against ≈ 50 today for the biggest exponentials and ≈ 100–
700 per output term for the level-1 products of idea 1. Crossover D ≈ 30 for h¹¹ ≈ 8–9; 10–100× at
D ≥ 40 — but memory (6 GB at degree 28 for h¹¹ = 10 in the C) makes that regime rare. Real use:
small h¹¹ (≤ 4–5) at very high degree (DKMMR's h¹¹ = 2 at degree 200 took 6 h in their code;
`h11_9_plike`-style gradings), and later the root-level products of idea 1, which are the largest
objects. Kronecker substitution into a box is hopeless (box/cone ratio ≈ 10⁶ at h¹¹ = 8, D = 20,
coordinate bounds up to 60 **[measured]**).

### Idea 5 — Lightcone-chunked extraction as a time–memory tradeoff

**Status: estimate.** Peak memory ∝ the touched set (366 k points at `h11_8`@20 but growing like
D^{h¹¹}: the C reports 3–6 GB at degree 26–28 for h¹¹ = 10–11). Partition the top-degree classes into
groups, run the extraction restricted (idea 2) to the union of each group's lightcones, and merge.
Memory ∝ the largest group's lightcone; time ∝ Σ_groups (restricted cost), i.e. the overlap of the
lightcones (all share the low-degree region) is the overhead — plausibly 2–5× time for 5–20× memory
**[speculation]**. Only worth it in the memory-bound regime; idea 1's structure (node tables)
complicates it.

### Idea 6 — Symmetries of the polytope

Automorphisms of the polytope that preserve the fan, q (up to row permutation), κ and the grading
permute curve classes and leave GV invariant; exponentials of an orbit are permutations of one
another. Gain ≤ |G| on the extraction's exponential work; typical Kreuzer–Skarke polytopes at h¹¹ ≈
10 have trivial or ℤ₂ groups, special ones (the 491 polytope) more. Idea 1's sharing already captures
most of what symmetry would give. Low priority **[estimate]**.

### Constant factors worth knowing (not big ideas)

- Lanes: a run needing 5–7 primes does two passes because `NL ≤ 4`; `NL` up to 8 (compile time
  ≈ +8 s per lane count) makes it one pass, ≈ 1.6× on such runs **[estimate from PLAN 6.5]**.
- SIMD modular arithmetic: 62-bit Montgomery is scalar on x86; 50-bit primes with the double-precision
  FMA mulmod (or AVX-512 IFMA 52-bit) vectorize the 3–4 lanes of a slot, ≈ 1.5–2.5× on the arithmetic
  half of the scatter, at the price of 20% more primes **[estimate]**.
- The lane probe at 3/4 of the degree costs ≈ (3/4)^{4.5} ≈ 27% of a run **[estimate from the
  measured degree scaling 18.1 → 27.5 s for 20 → 22]**; running the first pass at 4 lanes instead
  costs about the same, so nothing to gain here.
- 64-bit keys where they fit, structure-of-arrays level tables (PLAN 5.3).

---

## 3. Ideas considered and rejected

| idea | why rejected |
|---|---|
| Kronecker substitution / box or TFT NTT for the whole series | box/cone blow-up ≈ 10⁶ at h¹¹ = 8, D = 20 **[measured bounds]**; TFT swells 2^d; vdH–L: "hopeless for n ≥ 6" **[literature]** |
| Dense multivariate exp/log/composition by Newton | gains nothing over the Euler recurrence with schoolbook products; multivariate doubling multiplies size by 2^n so only the last step counts **[literature, subagent]** |
| Taylor expansion G = Σ_μ (−β)^μ/μ! θ^μ I | C(D/δ + r, r) ≈ 3e6 full products at r = 8, D = 20 **[estimate]** |
| Horner / triangular single-coordinate composition (q = q^(n)∘…∘q^(1)) | ≈ n(n+1)/2 · D full products ≈ 1e12 pairs at `h11_8`@20 vs 5e8 today; wins only at D ≳ 50 **[estimate, subagent]** |
| Newton inversion of the mirror map then one composition | composition has the same per-term exponential structure (z = q e^{−β}); ≥ 1× today's cost **[estimate]** |
| Evaluation–interpolation over F_p (Prony / Ben-Or–Tiwari) | truncation does not commute with evaluation; in the transform model the composition is a multivariate multipoint evaluation at τ-adic points (Kedlaya–Umans territory) **[literature]** |
| p-adic / Hensel lifting instead of CRT | 62-bit primes invert every denominator; the pipeline has no reusable Jacobian; mod p^k words cost what k primes cost; lanes already amortize the table work **[literature: Lairez–Vaccon; estimate]** |
| GKZ / Picard–Fuchs ideal, Yukawa couplings Y(z) rational, Birkhoff factorization, quasimaps | all need z(q) and a composition; the PF ideal is intractable at h¹¹ ≳ 4 (HKTY, HLY, DKMMR); Yukawa needs h¹¹³/6 rational functions with the A-discriminant in the denominators **[literature]** |
| Avoiding intersection numbers | they come from the fan in ≈ 0.3 s (`phase.rs`); only K_w = Σ_t w_t κ_t·· enters (threefolds); reconstructing them from consistency (INSTANTON does this) costs more |
| WDVV / associativity recursion | automatic for CY3 (only H²×H²→H⁴), gives no recursion **[literature]** |
| Track I only on candidate points / prune by reachability of candidates | invalid: 169 440 of the 225 207 nonzero n_C of `h11_8`@20 are at non-candidate (composite) classes **[measured]** |
| Topological vertex / localization / blowup equations for the compact X | non-compact toric only, or exponential in degree; local pieces give the face results cygv already computes by restriction **[literature]** |
| Modular / Noether–Lefschetz / Jacobi-form bootstraps | special fibrations, fixed base degree; not a general algorithm **[literature]** |
| Reducing α's support by a change of coordinates (rows) | the Euler recurrence already drops the zero L-terms, so `rows` changes nothing for the exact run **[measured: \|L_C\| ≪ \|α\|]** |
| Hosono's scalar prepotential F = −½∫ log W(x(q), J) | needs x(q) and an h¹¹²-product cubic contraction; strictly more than the gradient route **[literature, subagent verified the formula]** |

---

## 4. Answers to the specific questions

**Purely toric data.** Yes in the sense that fan → (rays, SR ideal) → κ, Mori cone, strata, grading is
exactly what `phase.rs` does, cheaply even at h¹¹ = 491 (0.3 s for the wall classes, 0.26 s for κ;
only the Mori cone's double description is expensive there and is skipped when not needed). No
GKZ/PF/I-function route removes κ: the genus-0 invariants depend on the cohomology ring of X, which
enters through the ħ⁻² term of the I-function / the Frobenius second derivatives. A notable
simplification is that the modular pipeline uses only the single symmetric matrix
K_w = Σ_t w_t κ_t·· (and S₂ ∝ Q_s K_w Q_u at two-negative points), not κ itself.

**Different extraction bases / triangular systems / Lagrange inversion / Newton.** The extraction is
the composition T⁻¹; the only structurally different exact formulations are (a) Lagrange–Good
(idea 3, per target) and (b) the factorization exp(C·α) = Π_r E_r^{Q_r·C} with shared Horner
products (idea 1). Logarithmic-derivative formulations (Euler operator identities) relate T of
different series but never avoid T.

**Fast truncated multivariate arithmetic.** At h¹¹ ≈ 8–11, D ≈ 20 the cone truncation makes
schoolbook products cost ≈ 50–700 operations per output coefficient; quasi-linear methods cost
≈ 2 000 (vdH–L's measured constants, subagent estimate). Crossover D ≈ 25–30; see idea 4.

**Number theory.** CRT with word-size primes, early termination and the certifier is the right
design; p-adic lifting has no advantage here. More lanes per pass is the one lever.

**Structure.** The exploitable structure is the per-divisor factorization of the mirror map (idea 1),
lightcones (ideas 2, 3, 5), and — already in the code — faces, the ≤ 2-negative candidate cones and
the vanishing-form skip. Closed forms (toric curves, contractible curves, Weyl curves: DKMMR §6,
Bryan–Katz–Leung) cover single classes and short rays, not the bulk.

---

## 5. Prototypes and what they show

1. `measure/` **[measured]**: the size table and cost models of sections 1–2 (Euler pairs reproduce
   PLAN's counted 7.23e7 scatter steps for `h11_11`@18 to within 3%, which validates the model), the
   per-curve product-tree costs, the sharing counts, the lightcone-pruning counts, the layer-level
   and cross-layer (relaxed) factorization costs. Run: `MEASURE_SKIP_MORI=1 MEASURE_LAYER_SAMPLE=1000000
   measure/target/release/measure benches/data/h11_8.yaml benches/data/cones/h11_8.yaml 20 1`.
2. `tree_extraction_poc.py` **[verified]**: exact-rational implementation of idea 1's relaxed Horner
   tree (rows ordered, nodes = exponent prefixes, degree-online finalize/read/insert/push) on random
   data (3 coordinates, 5 rows, degrees 8–11, several seeds): recovers the planted n_C exactly and
   agrees with the per-curve extraction. (On these unstructured toy inputs the tree does *not* save
   pairs — the saving comes from the real geometry's shared exponents and small cones, which is why
   the cost model was run on the real data.)
3. `lagrange_good_check.py` **[verified]**: idea 3's formula, including that the determinant is
   necessary.

---

## 6. Literature (from three subagent searches; key items, all with sources)

Established, read in the sources by the subagents unless marked otherwise:

- HKTY hep-th/9308122, hep-th/9406055 (Yukawa route via PF: "exceedingly tedious" beyond few
  moduli); HLY alg-geom/9511001 (Gröbner bases of the toric ideal; extra operators needed).
- DKMMR "Computational Mirror Symmetry" arXiv:2303.00757: Algorithm 1 is degree-by-degree
  subtraction; PF reduction "intractable even for modest numbers of moduli"; faces for h¹¹ ≳ 100;
  Table 4 (h¹¹ = 491, degree 25: 1.7 M classes, 627 nonzero, 9 min); §6 closed forms for toric curves.
- Givental alg-geom/9701016; Coates–Givental math/0110142; Coates–Corti–Iritani–Tseng math/0702234
  (J = I/F, τ = G/F; "we do not know how to invert the mirror map"); Guest math/0206212 and
  Guest–Sakai 0810.4236 (Birkhoff factorization needs D-module Gröbner bases); Iritani math/0411111;
  J. Wang arXiv:1910.14440 (non-convex phases, used for vex strata).
- Lagrange–Good: Good, Proc. Cambridge Philos. Soc. 56 (1960); Gessel, JCTA 45 (1987) and
  arXiv:1609.05988; Bender–Richmond, EJC 5 R33 (1998).
- Fast series arithmetic: Lecerf–Schost 2003 (quasi-linear total-degree products, char 0); van der
  Hoeven "Relax, but don't be too lazy" JSC 34 (2002) (weighted slicing, relaxed products); van der
  Hoeven–Lecerf, JSC 2013 "On the bit-complexity of sparse polynomial and series multiplication"
  (finite-field version, measured 10× *slower* than naive at n = 5, δ = 20); van der Hoeven, Notes
  on the TFT (multivariate TFT swells 2^d); Brent–Kung JACM 1978 and the 1977 multivariate paper;
  Kinoshita–Li arXiv:2404.05177 (near-linear composition, univariate only); Neiger–Salvy–Schost–
  Villard arXiv:2110.08354; Arnold–Roche arXiv:1501.05296; Giorgi–Grenet–Perret du Cray
  arXiv:2001.11959; Nakos arXiv:1901.09355.
- p-adic: Lairez–Vaccon arXiv:1602.00244; Berthomieu–van der Hoeven–Lecerf (relaxed p-adics);
  Kontsevich–Schwarz–Vologodsky hep-th/0603106 (integrality, not an algorithm).
- Structure: CKYZ hep-th/9903053 (local pieces inside compact CYs); Klemm–Maulik–Pandharipande–
  Scheidegger arXiv:0807.2477 (K3 fibrations, Noether–Lefschetz); Huang–Katz–Klemm arXiv:1501.04891
  (elliptic fibrations, Jacobi forms); Bryan–Katz–Leung math/9911056 (contractible curves);
  arXiv:2107.09064 and arXiv:2406.13751 (large-h¹¹ production use: potent rays, faces);
  MacFadden–Sheridan arXiv:2512.14817 (vex phases); Lukas–Ruehle arXiv:2210.09369 and
  Candelas–de la Ossa–Kuusela–McGovern arXiv:2111.02440 (symmetries).
- Software: INSTANTON (Klemm/Kreuzer, Mathematica; explicit inversion of the mirror map, hep-th/0410018);
  Gathmann's GROWI (hypersurfaces in Pⁿ); CYTools (`compute_gvs` wraps cygv).

What the literature does *not* contain: a fast multivariate degree-truncated composition, a
many-moduli PF/Yukawa algorithm, or any genus-0 recursion in q for CY3.
