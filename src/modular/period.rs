//! The fundamental period and its derivatives at each candidate.
//!
//! At a curve class $C$ with $d_r = C \cdot Q_r$ and $d_{0,j} = C \cdot Q_{0,j}$, the
//! coefficient of the fundamental period is
//! $c_0 = \prod_j d_{0,j}! / \prod_r d_r!$, which vanishes as soon as some $d_r < 0$, and the
//! derivatives with respect to the $h^{1,1}$ deformation parameters give $c_{1,a}$ and
//! $c_{2,ab}$. The formulas are those of the reference pipeline's `compute_c_*neg`
//! (`reference::fundamental_period`), with three differences, following the C program:
//!
//! - factorials and harmonic numbers come from tables ([`Tables`]) instead of being computed
//!   per curve;
//! - the second derivatives are contracted with the bilinear forms right away,
//!   $S_2^{(k)} = \sum_{a \le b} K^{(k)}_{ab}\, c_{2,ab}$ with half weight on the diagonal,
//!   instead of being kept per pair of indices;
//! - curves with negative anticanonical degree are rejected when they have fewer than two
//!   negative intersections, instead of wrapping around, and contribute nothing otherwise,
//!   except in a vex phase. There a curve of anticanonical degree $m < 0$ whose negative
//!   intersections $-k_i$ are with the divisors of a stratum $S$ adds $U_0\, (w \cdot d_S)$ to
//!   $S_2$, with
//!   $U_0 = (-1)^{e} \prod_{i \in S} (k_i - 1)! \,/\, \big(\prod_{i \notin S} b_i!\,
//!   (\lvert m \rvert - 1)!\big)$, $e = \lvert m \rvert - 1 + \sum_{i \in S} (k_i - 1)$ and $b_i$
//!   the other intersections, as in the C program and the reference's `compute_c_vex`. (The
//!   reference adds $c_0^{-1} U_0\, d_{S,t}$ to each instanton correction $t$, which is the
//!   same once contracted with the grading vector.)

use super::field::Ring;
use super::key::KeyLayout;
use super::table::SparsePoly;
use super::{dot, Error, Geometry};
use rayon::prelude::*;

/// Factorials, their inverses, inverses and harmonic numbers of the integers below a bound.
#[derive(Debug, Clone)]
pub struct Tables<E> {
    /// $n!$
    pub fact: Vec<E>,
    /// $1 / n!$
    pub inv_fact: Vec<E>,
    /// $1 / n$, and zero for $n = 0$.
    pub inv: Vec<E>,
    /// $H_n = \sum_{i=1}^n 1/i$
    pub h1: Vec<E>,
    /// $H^{(2)}_n = \sum_{i=1}^n 1/i^2$
    pub h2: Vec<E>,
}

impl<E: Copy + Default> Tables<E> {
    /// The tables for $0 \le n <$ `len`, with `len` at least 3, so that $1/2$ is in them.
    pub fn new<R: Ring<Elem = E>>(ring: &R, len: usize) -> Result<Self, Error> {
        let len = len.max(3);
        let inv = ring.inverses(len).ok_or(Error::NotInvertible)?;
        let mut fact = vec![ring.one(); len];
        let mut inv_fact = vec![ring.one(); len];
        let mut h1 = vec![E::default(); len];
        let mut h2 = vec![E::default(); len];
        for n in 1..len {
            fact[n] = ring.mul(fact[n - 1], ring.int(n as i64));
            inv_fact[n] = ring.mul(inv_fact[n - 1], inv[n]);
            h1[n] = ring.add(h1[n - 1], inv[n]);
            h2[n] = ring.add(h2[n - 1], ring.mul(inv[n], inv[n]));
        }
        Ok(Self {
            fact,
            inv_fact,
            inv,
            h1,
            h2,
        })
    }

    /// The number of entries of each table.
    #[inline]
    pub fn len(&self) -> usize {
        self.fact.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.fact.is_empty()
    }
}

/// The fundamental period data on the candidates, as sparse series.
#[derive(Debug, Clone)]
pub struct PeriodData<K, E> {
    /// $c_0$, one value per term.
    pub c0: SparsePoly<K, E>,
    /// $c_{1,a}$, $h^{1,1}$ values per term.
    pub c1: SparsePoly<K, E>,
    /// $S_2^{(k)}$, one value per bilinear form.
    pub s2: SparsePoly<K, E>,
    /// If asked for, $c_1$ along the rows of the charges, the anticanonical ones and then those
    /// of the divisors, so that $c_{1,a} = \sum_j Q^0_{ja} c_{1,j} + \sum_r Q_{ra} c_{1,r}$
    /// (see [`period`]).
    pub c1_rows: Option<SparsePoly<K, E>>,
    /// The tables used, which the later stages share. They cover every degree up to the
    /// maximum degree.
    pub tables: Tables<E>,
}

/// The ring images of the geometry, and the tables, for evaluating the period at a curve.
struct Evaluator<'a, R: Ring> {
    ring: &'a R,
    geometry: &'a Geometry,
    tables: &'a Tables<R::Elem>,
    /// Per coordinate $a$, the divisors $r$ with $Q_{ra} \ne 0$ and the image of $Q_{ra}$.
    q_cols: Vec<Vec<(usize, R::Elem)>>,
    /// The images of the anticanonical charges, row-major.
    q0: Vec<R::Elem>,
    /// The images of the charges, row-major.
    q: Vec<R::Elem>,
    /// The pairs $a \le b$ on which some form does not vanish.
    pairs: Vec<(usize, usize)>,
    /// Per pair, the images of $K^{(k)}_{ab}$ for every $k$, halved when $a = b$.
    k_pairs: Vec<R::Elem>,
    /// Per stratum, the image of $w \cdot d_S$.
    strata: Vec<R::Elem>,
}

/// Scratch space for [`Evaluator::point`].
struct Scratch<E> {
    dq: Vec<i64>,
    dq0: Vec<i64>,
    a: Vec<E>,
}

impl<'a, R: Ring> Evaluator<'a, R> {
    fn new(
        ring: &'a R,
        geometry: &'a Geometry,
        tables: &'a Tables<R::Elem>,
    ) -> Result<Self, Error> {
        let h11 = geometry.h11();
        let ndiv = geometry.n_divisors();
        let q_cols = (0..h11)
            .map(|a| {
                (0..ndiv)
                    .filter(|&r| geometry.divisor_q(r)[a] != 0)
                    .map(|r| (r, ring.int(geometry.divisor_q(r)[a] as i64)))
                    .collect()
            })
            .collect();
        let q0 = (0..geometry.n_anticanonical())
            .flat_map(|j| geometry.divisor_q0(j).iter().map(|&x| ring.int(x as i64)))
            .collect();
        let q = (0..ndiv)
            .flat_map(|r| geometry.divisor_q(r).iter().map(|&x| ring.int(x as i64)))
            .collect();
        let half = tables.inv[2];
        let mut pairs = Vec::new();
        let mut k_pairs = Vec::new();
        for a in 0..h11 {
            for b in a..h11 {
                let ks: Vec<i64> = (0..geometry.n_forms())
                    .map(|k| geometry.form(k)[a * h11 + b])
                    .collect();
                if ks.iter().all(|&x| x == 0) {
                    continue;
                }
                pairs.push((a, b));
                for x in ks {
                    let x = ring.int(x);
                    k_pairs.push(if a == b { ring.mul(x, half) } else { x });
                }
            }
        }
        let mut strata = Vec::new();
        for stratum in geometry.strata() {
            let mut v = R::Elem::default();
            for &(num, den) in &stratum.weighted_class {
                let inv = ring.inv(ring.int(den)).ok_or(Error::NotInvertible)?;
                ring.fma(&mut v, ring.int(num), inv);
            }
            strata.push(v);
        }
        Ok(Self {
            ring,
            geometry,
            tables,
            q_cols,
            q0,
            q,
            pairs,
            k_pairs,
            strata,
        })
    }

    fn scratch(&self) -> Scratch<R::Elem> {
        Scratch {
            dq: vec![0; self.geometry.n_divisors()],
            dq0: vec![0; self.geometry.n_anticanonical()],
            a: vec![R::Elem::default(); self.geometry.h11()],
        }
    }

    /// Adds $\sum_{(a, b)} K^{(k)}_{ab}\, f(a, b)$ to `s2[k]` for every form $k$, with the
    /// diagonal halved.
    #[inline]
    fn contract(&self, s2: &mut [R::Elem], mut f: impl FnMut(usize, usize) -> R::Elem) {
        let m = s2.len();
        for (i, &(a, b)) in self.pairs.iter().enumerate() {
            let v = f(a, b);
            for (s, &k) in s2.iter_mut().zip(&self.k_pairs[i * m..(i + 1) * m]) {
                self.ring.fma(s, k, v);
            }
        }
    }

    /// Computes $c_0$, $c_{1,a}$ and $S_2^{(k)}$ at the curve class `c` into `c0`, `c1` and
    /// `s2`, which must have one, $h^{1,1}$ and $M$ entries, and $c_1$ along the rows into
    /// `c1_rows` unless it is empty.
    fn point(
        &self,
        c: &[i32],
        sc: &mut Scratch<R::Elem>,
        c0: &mut R::Elem,
        c1: &mut [R::Elem],
        s2: &mut [R::Elem],
        c1_rows: &mut [R::Elem],
    ) -> Result<(), Error> {
        let ring = self.ring;
        let t = self.tables;
        let g = self.geometry;
        let h11 = g.h11();
        let zero = R::Elem::default();
        *c0 = zero;
        c1.fill(zero);
        s2.fill(zero);
        c1_rows.fill(zero);

        let mut neg = [0usize; 3];
        let mut nneg = 0;
        for (r, dq) in sc.dq.iter_mut().enumerate() {
            *dq = dot(c, g.q(r));
            if *dq < 0 {
                if nneg < 3 {
                    neg[nneg] = r;
                }
                nneg += 1;
            }
        }
        let mut anticanonical_negative = false;
        for (j, dq0) in sc.dq0.iter_mut().enumerate() {
            *dq0 = dot(c, g.q0(j));
            anticanonical_negative |= *dq0 < 0;
        }
        if anticanonical_negative {
            if nneg < 2 {
                return Err(Error::NegativeAnticanonicalDegree);
            }
            // Strata only exist on hypersurfaces, which have a single anticanonical row.
            let stratum = if nneg == 3 { g.stratum(&neg) } else { None };
            if let Some(i) = stratum {
                let m = -sc.dq0[0];
                let mut u0 = t.inv_fact[(m - 1) as usize];
                let mut e = m - 1;
                for &d in &sc.dq {
                    let f = if d < 0 {
                        e += -d - 1;
                        t.fact[(-d - 1) as usize]
                    } else {
                        t.inv_fact[d as usize]
                    };
                    u0 = ring.mul(u0, f);
                }
                if e % 2 != 0 {
                    u0 = ring.neg(u0);
                }
                s2[0] = ring.mul(u0, self.strata[i]);
            }
            return Ok(());
        }
        if nneg > 2 {
            return Ok(());
        }

        // prod_j d0_j! / prod_r d_r!, where 1/d! for d < 0 is read as (-d - 1)!; the sign
        // that goes with it is applied below.
        let mut fact = ring.one();
        for &d in &sc.dq0 {
            fact = ring.mul(fact, t.fact[d as usize]);
        }
        for &d in &sc.dq {
            let f = if d >= 0 {
                t.inv_fact[d as usize]
            } else {
                t.fact[(-d - 1) as usize]
            };
            fact = ring.mul(fact, f);
        }

        // A_a = sum_j H_{d0_j} Q0_ja - sum_r H_{n_r} Q_ra, n_r = d_r, or -d_r - 1 if negative.
        if nneg < 2 {
            for (a, aa) in sc.a.iter_mut().enumerate() {
                let mut v = zero;
                for (j, &d) in sc.dq0.iter().enumerate() {
                    ring.fma(&mut v, t.h1[d as usize], self.q0[j * h11 + a]);
                }
                for &(r, qra) in &self.q_cols[a] {
                    let d = sc.dq[r];
                    let n = if d < 0 { -d - 1 } else { d };
                    v = ring.sub(v, ring.mul(t.h1[n as usize], qra));
                }
                *aa = v;
            }
        }
        let a = &sc.a;
        let q = |r: usize, a: usize| self.q[r * h11 + a];

        let n0 = sc.dq0.len();
        match nneg {
            0 => {
                *c0 = fact;
                for (c1a, &aa) in c1.iter_mut().zip(a) {
                    *c1a = ring.mul(fact, aa);
                }
                if !c1_rows.is_empty() {
                    for (v, &d) in c1_rows.iter_mut().zip(&sc.dq0) {
                        *v = ring.mul(fact, t.h1[d as usize]);
                    }
                    for (v, &d) in c1_rows[n0..].iter_mut().zip(&sc.dq) {
                        let n = if d < 0 { -d - 1 } else { d };
                        *v = ring.neg(ring.mul(fact, t.h1[n as usize]));
                    }
                }
                self.contract(s2, |x, y| {
                    let mut v = ring.mul(a[x], a[y]);
                    for (j, &d) in sc.dq0.iter().enumerate() {
                        let qq = ring.mul(self.q0[j * h11 + x], self.q0[j * h11 + y]);
                        v = ring.sub(v, ring.mul(t.h2[d as usize], qq));
                    }
                    for &(r, qrx) in &self.q_cols[x] {
                        let qry = q(r, y);
                        if !ring.is_zero(qry) {
                            ring.fma(&mut v, t.h2[sc.dq[r] as usize], ring.mul(qrx, qry));
                        }
                    }
                    v
                });
                for s in s2.iter_mut() {
                    *s = ring.mul(fact, *s);
                }
            }
            1 => {
                let s = neg[0];
                let fs = if sc.dq[s] % 2 == 0 {
                    ring.neg(fact)
                } else {
                    fact
                };
                for (x, c1a) in c1.iter_mut().enumerate() {
                    *c1a = ring.mul(fs, q(s, x));
                }
                if !c1_rows.is_empty() {
                    c1_rows[n0 + s] = fs;
                }
                self.contract(s2, |x, y| {
                    ring.add(ring.mul(a[x], q(s, y)), ring.mul(a[y], q(s, x)))
                });
                for v in s2.iter_mut() {
                    *v = ring.mul(fs, *v);
                }
            }
            _ => {
                let (s, u) = (neg[0], neg[1]);
                let fs = if (sc.dq[s] + sc.dq[u]) % 2 == 0 {
                    fact
                } else {
                    ring.neg(fact)
                };
                self.contract(s2, |x, y| {
                    ring.add(ring.mul(q(s, x), q(u, y)), ring.mul(q(s, y), q(u, x)))
                });
                for v in s2.iter_mut() {
                    *v = ring.mul(fs, *v);
                }
            }
        }
        Ok(())
    }
}

/// Computes the fundamental period data at the candidates, which come with their degrees in
/// increasing order of degree and then key, as [`candidates`](super::candidates) gives them.
///
/// With `rows`, $c_1$ is also given along the rows of the charges: in $c_{1,a} = c_0 A_a$,
/// $A_a = \sum_j H_{Q^0_j \cdot C}\, Q^0_{ja} - \sum_r H_{n_r} Q_{ra}$ is split into its terms, and
/// with one negative intersection, $Q_s \cdot C < 0$, only row $s$ is nonzero. Exactly, that
/// is the same series; for the majorant run it is what keeps the exponents of the extraction
/// tight (see [`extract`](super::extraction::extract)).
pub fn period<R: Ring, L: KeyLayout>(
    ring: &R,
    layout: &L,
    geometry: &Geometry,
    candidates: &[(L::Key, u32)],
    max_deg: u32,
    rows: bool,
) -> Result<PeriodData<L::Key, R::Elem>, Error> {
    let (h11, dim) = (geometry.h11(), geometry.dim());
    if layout.dim() != dim {
        return Err(Error::DimensionMismatch);
    }
    let chunk = candidates
        .len()
        .div_ceil(4 * rayon::current_num_threads())
        .max(256);

    // The tables must reach every |C . Q_r| and |C . Q0_j| that is looked up, and every degree.
    let max_arg = candidates
        .par_chunks(chunk)
        .map(|cands| {
            let mut c = vec![0; dim];
            let mut max_arg = 0;
            for &(key, _) in cands {
                layout.unpack(key, &mut c);
                let rows = (0..geometry.n_divisors())
                    .map(|r| geometry.q(r))
                    .chain((0..geometry.n_anticanonical()).map(|j| geometry.q0(j)));
                for row in rows {
                    max_arg = max_arg.max(dot(&c, row).unsigned_abs());
                }
            }
            max_arg
        })
        .max()
        .unwrap_or(0)
        .max(max_deg as u64);
    let tables = Tables::new(ring, max_arg as usize + 2)?;

    let eval = Evaluator::new(ring, geometry, &tables)?;
    let m = geometry.n_forms();
    let na = if rows {
        geometry.n_anticanonical() + geometry.n_divisors()
    } else {
        0
    };
    // Each chunk of candidates gives the terms of a range of keys, in order.
    let parts = candidates
        .par_chunks(chunk)
        .map(|cands| {
            crate::interrupt::check(Error::Interrupted)?;
            let mut c = vec![0; dim];
            let mut sc = eval.scratch();
            let mut c0 = SparsePoly::new(1);
            let mut c1 = SparsePoly::new(h11);
            let mut s2 = SparsePoly::new(m);
            let mut c1r = SparsePoly::new(na);
            let mut v0 = R::Elem::default();
            let mut v1 = vec![R::Elem::default(); h11];
            let mut v2 = vec![R::Elem::default(); m];
            let mut vr = vec![R::Elem::default(); na];
            for &(key, deg) in cands {
                layout.unpack(key, &mut c);
                eval.point(&c, &mut sc, &mut v0, &mut v1, &mut v2, &mut vr)?;
                if !ring.is_zero(v0) {
                    c0.push(key, deg, &[v0]);
                }
                if v1.iter().any(|&x| !ring.is_zero(x)) {
                    c1.push(key, deg, &v1);
                }
                if v2.iter().any(|&x| !ring.is_zero(x)) {
                    s2.push(key, deg, &v2);
                }
                if vr.iter().any(|&x| !ring.is_zero(x)) {
                    c1r.push(key, deg, &vr);
                }
            }
            Ok((c0, c1, s2, c1r))
        })
        .collect::<Result<Vec<_>, Error>>()?;
    let (mut p0, mut p1, mut p2, mut pr) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for (a, b, c, d) in parts {
        p0.push(a);
        p1.push(b);
        p2.push(c);
        pr.push(d);
    }
    let c0 = SparsePoly::concat(1, p0);
    let c1 = SparsePoly::concat(h11, p1);
    let s2 = SparsePoly::concat(m, p2);
    let c1_rows = SparsePoly::concat(na, pr);
    Ok(PeriodData {
        c0,
        c1,
        s2,
        c1_rows: rows.then_some(c1_rows),
        tables,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modular::field::{primes, LaneField};
    use crate::modular::residue;
    use crate::CYKind;
    use malachite::{Integer, Rational};

    #[test]
    fn test_tables() {
        let ps = [primes().next().unwrap(), 1_000_003];
        let lf = LaneField::new(ps);
        let tables = Tables::new(&lf, 40).unwrap();
        assert_eq!(tables.len(), 40);
        let mut fact = Integer::from(1);
        let (mut h1, mut h2) = (Rational::default(), Rational::default());
        for n in 0..40u32 {
            if n > 0 {
                fact *= Integer::from(n);
                h1 += Rational::from_unsigneds(1, n);
                h2 += Rational::from_unsigneds(1, n * n);
            }
            for (l, f) in lf.fields().iter().enumerate() {
                let p = f.p();
                let at = |v: &Vec<_>| lf.residues(v[n as usize])[l];
                assert_eq!(at(&tables.fact), residue(&Rational::from(fact.clone()), p));
                assert_eq!(
                    at(&tables.inv_fact),
                    residue(&Rational::from_integers(Integer::from(1), fact.clone()), p)
                );
                let inv = if n == 0 {
                    0
                } else {
                    residue(&Rational::from_unsigneds(1, n), p)
                };
                assert_eq!(at(&tables.inv), inv);
                assert_eq!(at(&tables.h1), residue(&h1, p));
                assert_eq!(at(&tables.h2), residue(&h2, p));
            }
        }
        // At least up to 2, whatever is asked for.
        assert_eq!(Tables::new(&lf, 0).unwrap().len(), 3);
        // 7! is zero modulo 7, so the table cannot reach it.
        assert!(Tables::new(&LaneField::new([7]), 8).is_err());
    }

    #[test]
    fn test_negative_anticanonical_degree() {
        use crate::modular::key::{KeyLayout, PackedLayout};
        let lf = LaneField::new([primes().next().unwrap()]);
        let layout = PackedLayout::<u128>::new(&[2]).unwrap();
        let cands = [
            (layout.pack(&[0]), 0),
            (layout.pack(&[1]), 1),
            (layout.pack(&[2]), 2),
        ];
        // C . Q0 = -1 with no negative intersection.
        let g = Geometry::new(CYKind::Threefold, 1, vec![1], vec![-1], vec![1]).unwrap();
        assert_eq!(
            period(&lf, &layout, &g, &cands, 2, false).unwrap_err(),
            Error::NegativeAnticanonicalDegree
        );
        // With two negative intersections the curve contributes nothing.
        let g = Geometry::new(CYKind::Threefold, 1, vec![-1, -1], vec![-2], vec![1]).unwrap();
        let pd = period(&lf, &layout, &g, &cands, 2, false).unwrap();
        assert_eq!(pd.c0.len(), 1);
        assert!(pd.c1.is_empty() && pd.s2.is_empty());
    }

    /// $c_1$ along the rows of the charges gives $c_1$ back, on the toy threefold.
    #[test]
    fn test_c1_rows() {
        use crate::modular::candidates::candidates_by_enumeration;
        use crate::modular::key::PackedLayout;
        use nalgebra::DMatrix;
        use std::collections::HashMap;
        let q = DMatrix::from_column_slice(6, 2, &[1, 1, 1, 0, 1, 2, 0, 0, -1, 1, 1, -1]);
        let intnums = HashMap::from([
            ((0, 0, 0), 2),
            ((1, 0, 0), 1),
            ((0, 1, 1), -1),
            ((1, 1, 1), 5),
        ]);
        let w = [3, -1];
        let g = Geometry::from_toric(&q, &[], &intnums, CYKind::Threefold, &w).unwrap();
        let layout = PackedLayout::<u128>::new(&[20, 20]).unwrap();
        let gens: [&[i32]; 2] = [&[0, -1], &[1, 2]];
        let cands = candidates_by_enumeration(&layout, gens, &w, 8, &g).unwrap();
        let lf = LaneField::new([primes().next().unwrap()]);
        assert!(period(&lf, &layout, &g, &cands, 8, false)
            .unwrap()
            .c1_rows
            .is_none());
        let pd = period(&lf, &layout, &g, &cands, 8, true).unwrap();
        let rows = pd.c1_rows.as_ref().unwrap();
        assert_eq!(rows.width(), g.n_anticanonical() + g.n_divisors());
        let terms = |p: &SparsePoly<_, _>| -> HashMap<_, Vec<_>> {
            (0..p.len())
                .map(|i| (p.keys()[i], p.values(i).to_vec()))
                .collect()
        };
        let (c1, rows) = (terms(&pd.c1), terms(rows));
        assert!(!c1.is_empty());
        let m: Vec<&[i32]> = (0..g.n_anticanonical())
            .map(|j| g.q0(j))
            .chain((0..g.n_divisors()).map(|r| g.q(r)))
            .collect();
        for (key, row_vals) in &rows {
            let expected: Vec<_> = (0..2)
                .map(|a| {
                    let mut v = Default::default();
                    for (row, &x) in m.iter().zip(row_vals) {
                        lf.fma(&mut v, lf.int(row[a] as i64), x);
                    }
                    v
                })
                .collect();
            let got = c1.get(key).cloned().unwrap_or(vec![Default::default(); 2]);
            assert_eq!(got, expected);
        }
        assert!(c1.keys().all(|k| rows.contains_key(k)));
    }
}
