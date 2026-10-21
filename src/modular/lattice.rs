//! Curve classes in a reduced basis of the lattice they span.
//!
//! The curve classes of a run are sums of the generators of its candidates, so they lie in
//! the lattice $\Lambda$ those span. On a face of the Mori cone that lattice has a lower rank
//! than $h^{1,1}$, and the coordinates in the input basis are many and large where a basis
//! of $\Lambda$ needs few and small ones: a face of rank 24 at $h^{1,1} = 491$ needs about
//! 300 bits in the input coordinates at degree 20, and 120 in an LLL-reduced basis of
//! $\Lambda$. So when the rank is lower, the run keys the curves by their coordinates in such
//! a basis, as the C program of PR #102 does (`setup_lattice`): the curve side of the
//! charges and the grading vector become their values on the basis, the lightcones'
//! inequalities too, and the output goes back to the input coordinates. Unlike a
//! projection this is exact: a class and its coordinates determine each other, and since a
//! multiple $n C'$ of a class of $\Lambda$ has $n$ dividing its coordinates exactly when
//! $C' \in \Lambda$, the multicover's divisors are the same (a $C/n$ outside $\Lambda$ is no
//! sum of generators, and has no invariant to contribute).
//!
//! [`Lattice::spanned`] takes the vectors one at a time, shortest first. A vector whose
//! coordinates in the basis so far are found, in floating point against the inverse of the
//! Gram matrix and then checked exactly, is in the lattice already. One independent of the
//! basis modulo a prime of 61 bits, and so over the rationals, joins it. One in its span
//! with fractional coordinates $a / D$ (exactly, by rationals) joins it too: with a
//! unimodular $U$ taking $a$ to $(g, 0, \dots, 0)$, the lattice and the vector $v$ span the
//! lattice of basis $v / g$ and the rows $2, \dots, k$ of $U^{-1} B$. After each, the basis is
//! reduced again by the integral LLL algorithm (Cohen, *A Course in Computational Algebraic
//! Number Theory*, algorithm 2.6.7) with $\delta = 0.99$, the C's, exactly, so the basis does
//! not depend on the platform. Where the C takes an echelon form of all the vectors, in `i64`,
//! this only ever handles the basis and one more vector, since echelon forms can grow far
//! beyond their input.

use super::field::Ring;
use super::table::SparsePoly;
use super::Error;
use malachite::base::num::arithmetic::traits::{DivMod, Gcd, Lcm};
use malachite::{Integer, Natural, Rational};
use std::collections::HashSet;
use std::sync::OnceLock;

/// A basis of the lattice spanned by some vectors of $\mathbb{Z}^{h^{1,1}}$, and the
/// coordinates in it.
#[derive(Debug, Clone)]
pub struct Lattice {
    h11: usize,
    /// The reduced basis, one row of `h11` entries per basis vector.
    basis: Vec<i64>,
    /// The inverse of its Gram matrix $G = B B^T$, row-major.
    gram_inverse: Vec<f64>,
    /// The same exactly, as $N / d$, computed when first needed.
    exact_inverse: OnceLock<(Vec<Integer>, Integer)>,
}

/// An overflow in computing a basis.
fn overflow() -> Error {
    Error::LatticeTooLarge
}

/// A row of big integers.
type Row = Vec<Integer>;

fn row(v: &[i32]) -> Row {
    v.iter().map(|&x| Integer::from(x)).collect()
}

/// Where a vector lies with respect to a lattice.
enum Position {
    /// In the lattice, with these coordinates.
    In(Vec<Integer>),
    /// In its span, with coordinates $a / D$ for some $D > 1$ with $\gcd(a, D) = 1$; this is
    /// $a$.
    Span(Vec<Integer>),
    /// Outside its span.
    Outside,
}

impl Lattice {
    /// The lattice spanned by `vectors`, each of `h11` entries, if it has a rank below `h11`
    /// (zero included), and `None` if it has full rank.
    pub fn spanned<'a>(
        vectors: impl IntoIterator<Item = &'a [i32]>,
        h11: usize,
    ) -> Result<Option<Self>, Error> {
        let mut seen = HashSet::new();
        let mut vectors: Vec<&[i32]> = vectors.into_iter().filter(|&v| seen.insert(v)).collect();
        if vectors.iter().any(|v| v.len() != h11) {
            return Err(Error::DimensionMismatch);
        }
        // The shortest first, which makes a basis more likely to span the rest.
        vectors.sort_by_key(|v| v.iter().map(|&x| x as i64 * x as i64).sum::<i64>());
        let mut lattice = Self::new(Vec::new(), h11)?;
        let mut span = Span::new(h11);
        for v in vectors {
            if lattice.fast_coords(v).is_some() {
                continue;
            }
            let rows: Vec<Row> = if span.add(v) {
                if span.rank() == h11 {
                    return Ok(None);
                }
                let mut rows = lattice.rows();
                rows.push(row(v));
                rows
            } else {
                match lattice.position(v) {
                    Position::In(_) => continue,
                    Position::Span(a) => lattice.joined(v, a),
                    // The rank modulo the prime was too low, which a prime of 61 bits makes
                    // vanishingly unlikely.
                    Position::Outside => return Err(Error::LatticeTooLarge),
                }
            };
            lattice = Self::new(lll(rows), h11)?;
        }
        Ok(Some(lattice))
    }

    /// The lattice with the independent `rows` as its basis.
    fn new(rows: Vec<Row>, h11: usize) -> Result<Self, Error> {
        let k = rows.len();
        let basis: Vec<i64> = rows
            .iter()
            .flatten()
            .map(|x| i64::try_from(x).map_err(|_| overflow()))
            .collect::<Result<_, _>>()?;
        let mut gram = vec![0f64; k * k];
        for i in 0..k {
            for j in 0..k {
                let (a, b) = (
                    &basis[i * h11..(i + 1) * h11],
                    &basis[j * h11..(j + 1) * h11],
                );
                gram[i * k + j] = a.iter().zip(b).map(|(&x, &y)| (x * y) as f64).sum();
            }
        }
        Ok(Self {
            h11,
            basis,
            gram_inverse: invert(gram, k),
            exact_inverse: OnceLock::new(),
        })
    }

    /// The basis as rows of big integers.
    fn rows(&self) -> Vec<Row> {
        (0..self.rank())
            .map(|j| {
                self.basis_vector(j)
                    .iter()
                    .map(|&x| Integer::from(x))
                    .collect()
            })
            .collect()
    }

    /// The rank of the lattice, i.e. the number of coordinates of a class.
    pub fn rank(&self) -> usize {
        self.basis.len().checked_div(self.h11).unwrap_or(0)
    }

    /// The dimension $h^{1,1}$ of the space it lies in.
    pub fn h11(&self) -> usize {
        self.h11
    }

    /// Basis vector `j`, in the input coordinates.
    pub fn basis_vector(&self, j: usize) -> &[i64] {
        &self.basis[j * self.h11..(j + 1) * self.h11]
    }

    /// The coordinates of `x` in the reduced basis, or `None` if it is not in the lattice or
    /// they do not fit an `i32`.
    pub fn coords(&self, x: &[i32]) -> Option<Vec<i32>> {
        if x.len() != self.h11 {
            return None;
        }
        if let Some(c) = self.fast_coords(x) {
            return c.into_iter().map(|c| i32::try_from(c).ok()).collect();
        }
        match self.position(x) {
            Position::In(c) => c.iter().map(|c| i32::try_from(c).ok()).collect(),
            _ => None,
        }
    }

    /// The coordinates of `x`, from those in floating point rounded, if they give `x` exactly.
    fn fast_coords(&self, x: &[i32]) -> Option<Vec<i64>> {
        let (n, k) = (self.h11, self.rank());
        // c G = x B^T, in floating point, and then exactly c B = x.
        let nonzero: Vec<(usize, i64)> = x
            .iter()
            .enumerate()
            .filter(|&(_, &v)| v != 0)
            .map(|(t, &v)| (t, v as i64))
            .collect();
        let y: Vec<f64> = (0..k)
            .map(|j| {
                let b = self.basis_vector(j);
                nonzero.iter().map(|&(t, v)| (v * b[t]) as f64).sum()
            })
            .collect();
        let mut c = vec![0i64; k];
        for (j, cj) in c.iter_mut().enumerate() {
            let v: f64 = (0..k).map(|i| y[i] * self.gram_inverse[i * k + j]).sum();
            if !v.is_finite() || v.abs() > (1u64 << 52) as f64 {
                return None;
            }
            *cj = v.round() as i64;
        }
        let mut s: Vec<i128> = x.iter().map(|&v| v as i128).collect();
        for (j, &cj) in c.iter().enumerate().filter(|&(_, &cj)| cj != 0) {
            for (st, &b) in s.iter_mut().zip(&self.basis[j * n..(j + 1) * n]) {
                *st -= cj as i128 * b as i128;
            }
        }
        s.iter().all(|&r| r == 0).then_some(c)
    }

    /// Where `x` lies, exactly: its coordinates are $c = x B^T G^{-1}$ if it is in the span.
    fn position(&self, x: &[i32]) -> Position {
        let (n, k) = (self.h11, self.rank());
        let (inv, d) = self
            .exact_inverse
            .get_or_init(|| exact_inverse(&self.rows()));
        let y: Vec<Integer> = (0..k)
            .map(|j| {
                let b = self.basis_vector(j);
                Integer::from(
                    x.iter()
                        .zip(b)
                        .map(|(&v, &b)| v as i128 * b as i128)
                        .sum::<i128>(),
                )
            })
            .collect();
        // c = a' / d with a' = y N.
        let a: Vec<Integer> = (0..k)
            .map(|j| (0..k).map(|i| &y[i] * &inv[i * k + j]).sum())
            .collect();
        // Back in the input coordinates, d x = a' B.
        let in_span = (0..n).all(|t| {
            let s: Integer = a
                .iter()
                .enumerate()
                .map(|(j, aj)| aj * Integer::from(self.basis[j * n + t]))
                .sum();
            s == d * Integer::from(x[t])
        });
        if !in_span {
            return Position::Outside;
        }
        let g = a
            .iter()
            .fold(Natural::from(0u32), |g, aj| g.gcd(aj.unsigned_abs_ref()))
            .gcd(d.unsigned_abs_ref());
        let g = Integer::from(g);
        let a: Vec<Integer> = a.into_iter().map(|aj| aj / &g).collect();
        let d = d / &g;
        if d == 1 {
            Position::In(a)
        } else {
            Position::Span(a)
        }
    }

    /// A basis of the lattice spanned by this one and `v`, of coordinates $a / D$ in it with
    /// $D > 1$ and $\gcd(a, D) = 1$.
    ///
    /// With $U$ unimodular and $a U = g e_1$, $g = \gcd(a)$, which is prime to $D$, the
    /// coordinates of the joined lattice, $\mathbb{Z}^k + \mathbb{Z} a / D$, times $D U$ are
    /// $\mathbb{Z} g e_1 + D \mathbb{Z}^k = \mathbb{Z} e_1 \oplus D \mathbb{Z}^{k-1}$. So it has
    /// the basis $e_1 U^{-1} / D = a / (g D)$, i.e. the vector $v / g$, and the rows
    /// $2, \dots, k$ of $U^{-1}$.
    fn joined(&self, v: &[i32], mut a: Vec<Integer>) -> Vec<Row> {
        let k = a.len();
        // U^{-1}, by the inverses of the column operations that take a to g e_1.
        let mut u_inv: Vec<Row> = (0..k)
            .map(|i| (0..k).map(|j| Integer::from((i == j) as i32)).collect())
            .collect();
        for j in 1..k {
            if a[j] == 0 {
                continue;
            }
            // (a_0, a_j) M = (g, 0) with M = [[x, -a_j/g], [y, a_0/g]], and rows 0 and j of
            // U^{-1} become M^{-1} = [[a_0/g, a_j/g], [-y, x]] times them.
            let (g, x, y) = extended_gcd(&a[0], &a[j]);
            let (s, t) = (&a[0] / &g, &a[j] / &g);
            let (r0, rj) = (u_inv[0].clone(), u_inv[j].clone());
            u_inv[0] = r0.iter().zip(&rj).map(|(p, q)| &s * p + &t * q).collect();
            u_inv[j] = r0.iter().zip(&rj).map(|(p, q)| &x * q - &y * p).collect();
            a[0] = g;
            a[j] = Integer::from(0);
        }
        let g = if a[0] < 0 { -&a[0] } else { a[0].clone() };
        let b = self.rows();
        let mut rows: Vec<Row> = vec![row(v).into_iter().map(|x| x / &g).collect()];
        for r in &u_inv[1..] {
            rows.push(
                (0..self.h11)
                    .map(|t| r.iter().zip(&b).map(|(c, bj)| c * &bj[t]).sum())
                    .collect(),
            );
        }
        rows
    }

    /// The class with coordinates `c` in the reduced basis, in the input coordinates, or
    /// `None` if they do not fit an `i32`.
    pub fn to_input(&self, c: &[i32], out: &mut [i32]) -> Option<()> {
        debug_assert_eq!((c.len(), out.len()), (self.rank(), self.h11));
        for (t, o) in out.iter_mut().enumerate() {
            let mut s = 0i64;
            for (j, &cj) in c.iter().enumerate() {
                if cj != 0 {
                    s = s.checked_add((cj as i64).checked_mul(self.basis[j * self.h11 + t])?)?;
                }
            }
            *o = i32::try_from(s).ok()?;
        }
        Some(())
    }

    /// The values of a linear form, given by its `h11` coefficients, on the reduced basis:
    /// its coefficients in the coordinates of the lattice.
    pub fn restrict(&self, form: &[i32]) -> Result<Vec<i32>, Error> {
        if form.len() != self.h11 {
            return Err(Error::DimensionMismatch);
        }
        (0..self.rank())
            .map(|j| {
                let s: i128 = self
                    .basis_vector(j)
                    .iter()
                    .zip(form)
                    .map(|(&b, &w)| b as i128 * w as i128)
                    .sum();
                i32::try_from(s).map_err(|_| overflow())
            })
            .collect()
    }
}

/// The restriction of each term of a series of $h^{1,1}$ components, in `ring`, to the lattice
/// of basis `basis` (one row of $h^{1,1}$ entries per basis vector): the coefficients of
/// $C \cdot \alpha = \sum_j c_j \alpha'_j$ for a class $C$ of coordinates $c$, with
/// $\alpha'_j = \sum_t b_{jt} \alpha_t$.
pub fn restrict_series<R: Ring, K: super::key::Key>(
    ring: &R,
    basis: &[i64],
    alpha: &SparsePoly<K, R::Elem>,
) -> Result<SparsePoly<K, R::Elem>, Error> {
    let h11 = alpha.width();
    if h11 == 0 || !basis.len().is_multiple_of(h11) {
        return Err(Error::DimensionMismatch);
    }
    let k = basis.len() / h11;
    // Per basis vector, its nonzero entries.
    let entries: Vec<Vec<(usize, R::Elem)>> = basis
        .chunks_exact(h11)
        .map(|b| {
            b.iter()
                .enumerate()
                .filter(|&(_, &x)| x != 0)
                .map(|(t, &x)| (t, ring.int(x)))
                .collect()
        })
        .collect();
    let mut out = SparsePoly::new(k);
    let mut v = vec![R::Elem::default(); k];
    for i in 0..alpha.len() {
        let a = alpha.values(i);
        for (vj, e) in v.iter_mut().zip(&entries) {
            *vj = R::Elem::default();
            for &(t, b) in e {
                ring.fma(vj, b, a[t]);
            }
        }
        if v.iter().any(|&x| !ring.is_zero(x)) {
            out.push(alpha.keys()[i], alpha.degs()[i], &v);
        }
    }
    Ok(out)
}

/// A prime of 61 bits, for the ranks.
const P: u64 = (1 << 61) - 1;

/// The span of some vectors modulo [`P`], as an echelon form with a one at each pivot.
struct Span {
    h11: usize,
    rows: Vec<(usize, Vec<u64>)>,
}

impl Span {
    fn new(h11: usize) -> Self {
        Self {
            h11,
            rows: Vec::new(),
        }
    }

    fn rank(&self) -> usize {
        self.rows.len()
    }

    /// Adds `v` if it is independent of the span modulo the prime, and so over the
    /// rationals, and tells whether it was.
    fn add(&mut self, v: &[i32]) -> bool {
        let mul = |a: u64, b: u64| ((a as u128 * b as u128) % P as u128) as u64;
        let mut v: Vec<u64> = v
            .iter()
            .map(|&x| (x as i64).rem_euclid(P as i64) as u64)
            .collect();
        debug_assert_eq!(v.len(), self.h11);
        for (p, r) in &self.rows {
            let f = v[*p];
            if f != 0 {
                for (x, &y) in v.iter_mut().zip(r) {
                    *x = (*x + P - mul(f, y)) % P;
                }
            }
        }
        let Some(p) = v.iter().position(|&x| x != 0) else {
            return false;
        };
        let inv = pow(v[p], P - 2);
        v.iter_mut().for_each(|x| *x = mul(*x, inv));
        self.rows.push((p, v));
        true
    }
}

/// $a^e \bmod P$.
fn pow(mut a: u64, mut e: u64) -> u64 {
    let mut r = 1u64;
    while e > 0 {
        if e & 1 == 1 {
            r = ((r as u128 * a as u128) % P as u128) as u64;
        }
        a = ((a as u128 * a as u128) % P as u128) as u64;
        e >>= 1;
    }
    r
}

/// The inverse of the symmetric positive definite `k` × `k` matrix `m`, row-major, by
/// Gauss–Jordan elimination with partial pivoting.
fn invert(mut m: Vec<f64>, k: usize) -> Vec<f64> {
    let mut inv = vec![0f64; k * k];
    for i in 0..k {
        inv[i * k + i] = 1.0;
    }
    for col in 0..k {
        let p = (col..k)
            .max_by(|&a, &b| m[a * k + col].abs().total_cmp(&m[b * k + col].abs()))
            .unwrap();
        for j in 0..k {
            m.swap(col * k + j, p * k + j);
            inv.swap(col * k + j, p * k + j);
        }
        let d = m[col * k + col];
        for j in 0..k {
            m[col * k + j] /= d;
            inv[col * k + j] /= d;
        }
        for r in (0..k).filter(|&r| r != col) {
            let f = m[r * k + col];
            if f != 0.0 {
                for j in 0..k {
                    m[r * k + j] -= f * m[col * k + j];
                    inv[r * k + j] -= f * inv[col * k + j];
                }
            }
        }
    }
    inv
}

/// The inverse of the Gram matrix of the independent `rows` exactly, as $N / d$ with $N$
/// integer, row-major, and $d > 0$.
fn exact_inverse(rows: &[Row]) -> (Vec<Integer>, Integer) {
    let k = rows.len();
    let dot =
        |x: &[Integer], y: &[Integer]| -> Integer { x.iter().zip(y).map(|(a, c)| a * c).sum() };
    let mut m: Vec<Rational> = Vec::with_capacity(k * k);
    for a in rows {
        for b in rows {
            m.push(Rational::from(dot(a, b)));
        }
    }
    let mut inv: Vec<Rational> = (0..k * k)
        .map(|i| Rational::from((i / k == i % k) as u32))
        .collect();
    // Gauss–Jordan; a Gram matrix of independent rows is positive definite, so the diagonal
    // stays nonzero without pivoting.
    for col in 0..k {
        let d = m[col * k + col].clone();
        for j in 0..k {
            m[col * k + j] /= &d;
            inv[col * k + j] /= &d;
        }
        for r in (0..k).filter(|&r| r != col) {
            let f = m[r * k + col].clone();
            if f != 0u32 {
                for j in 0..k {
                    let x = &f * &m[col * k + j];
                    m[r * k + j] -= x;
                    let x = &f * &inv[col * k + j];
                    inv[r * k + j] -= x;
                }
            }
        }
    }
    let d = inv
        .iter()
        .fold(Natural::from(1u32), |l, x| l.lcm(x.denominator_ref()));
    let n = inv
        .iter()
        .map(|x| {
            let (num, den) = x.numerator_and_denominator_ref();
            let num = Integer::from(num.clone()) * Integer::from(&d / den);
            if *x < 0u32 {
                -num
            } else {
                num
            }
        })
        .collect();
    (n, Integer::from(d))
}

/// $(g, x, y)$ with $g = \gcd(a, b) = x a + y b > 0$, for $(a, b) \ne (0, 0)$.
fn extended_gcd(a: &Integer, b: &Integer) -> (Integer, Integer, Integer) {
    let (mut r0, mut r1) = (a.clone(), b.clone());
    let (mut x0, mut x1, mut y0, mut y1) = (
        Integer::from(1),
        Integer::from(0),
        Integer::from(0),
        Integer::from(1),
    );
    while r1 != 0 {
        let q = (&r0).div_mod(&r1).0;
        let r = &r0 - &q * &r1;
        (r0, r1) = (r1, r);
        let x = &x0 - &q * &x1;
        (x0, x1) = (x1, x);
        let y = &y0 - &q * &y1;
        (y0, y1) = (y1, y);
    }
    if r0 < 0 {
        (-r0, -x0, -y0)
    } else {
        (r0, x0, y0)
    }
}

/// $\lfloor a / b + 1/2 \rfloor$ for $b > 0$.
fn round_div(a: &Integer, b: &Integer) -> Integer {
    (Integer::from(2) * a + b).div_mod(Integer::from(2) * b).0
}

/// The integral LLL reduction of the independent `rows`, with $\delta = 99/100$.
///
/// Cohen's algorithm 2.6.7, in which $d_i$ is the Gram determinant of the first $i$ rows and
/// $\lambda_{ij} = d_j \mu_{ij}$, all integers, so it runs in exact arithmetic.
fn lll(mut b: Vec<Row>) -> Vec<Row> {
    let k = b.len();
    if k == 0 {
        return b;
    }
    let dot =
        |x: &[Integer], y: &[Integer]| -> Integer { x.iter().zip(y).map(|(a, c)| a * c).sum() };
    // d[0] = 1 and d[i + 1] for row i; lambda[i][j] for j < i.
    let mut d = vec![Integer::from(1); k + 1];
    let mut lambda = vec![vec![Integer::from(0); k]; k];
    let (p, q) = (Integer::from(99), Integer::from(100));

    // The Gram–Schmidt data of row i, from those of the rows before it.
    let gram = |i: usize, b: &[Row], d: &mut [Integer], lambda: &mut [Vec<Integer>]| {
        for j in 0..=i {
            let mut u = dot(&b[i], &b[j]);
            for l in 0..j {
                u = (&d[l + 1] * u - &lambda[i][l] * &lambda[j][l]) / &d[l];
            }
            if j < i {
                lambda[i][j] = u;
            } else {
                d[i + 1] = u;
            }
        }
    };
    // b_i -= q b_l with q the rounded mu_il, for l < i, if that is not zero.
    let reduce = |i: usize, l: usize, b: &mut [Row], d: &[Integer], lambda: &mut [Vec<Integer>]| {
        let twice = Integer::from(2) * &lambda[i][l];
        let abs = if twice < 0 { -twice } else { twice };
        if abs <= d[l + 1] {
            return;
        }
        let qi = round_div(&lambda[i][l], &d[l + 1]);
        let (lo, hi) = b.split_at_mut(i);
        for (x, y) in hi[0].iter_mut().zip(&lo[l]) {
            *x -= &qi * y;
        }
        lambda[i][l] -= &qi * &d[l + 1];
        let (lo, hi) = lambda.split_at_mut(i);
        for (x, y) in hi[0][..l].iter_mut().zip(&lo[l][..l]) {
            *x -= &qi * y;
        }
    };

    gram(0, &b, &mut d, &mut lambda);
    let mut kmax = 0;
    let mut i = 1;
    while i < k {
        if i > kmax {
            kmax = i;
            gram(i, &b, &mut d, &mut lambda);
        }
        reduce(i, i - 1, &mut b, &d, &mut lambda);
        // Lovász: q d_i d_{i-2} < p d_{i-1}^2 - q lambda^2, in the indexing of d from 1.
        let lhs = &q * &d[i + 1] * &d[i - 1];
        let rhs = &p * &d[i] * &d[i] - &q * &lambda[i][i - 1] * &lambda[i][i - 1];
        if lhs < rhs {
            // Swap rows i - 1 and i.
            b.swap(i, i - 1);
            let (lo, hi) = lambda.split_at_mut(i);
            lo[i - 1][..i - 1].swap_with_slice(&mut hi[0][..i - 1]);
            let lam = lambda[i][i - 1].clone();
            let bb = (&d[i - 1] * &d[i + 1] + &lam * &lam) / &d[i];
            for row in &mut lambda[i + 1..=kmax] {
                let tv = row[i].clone();
                row[i] = (&d[i + 1] * &row[i - 1] - &lam * &tv) / &d[i];
                row[i - 1] = (&bb * &tv + &lam * &row[i]) / &d[i + 1];
            }
            d[i] = bb;
            i = (i - 1).max(1);
        } else {
            for l in (0..i - 1).rev() {
                reduce(i, l, &mut b, &d, &mut lambda);
            }
            i += 1;
        }
    }
    b
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modular::TestRng;

    /// Every class of the lattice, as a combination of the vectors, has coordinates that give
    /// it back, and the coordinates add.
    fn check(vectors: &[Vec<i32>], h11: usize, rank: usize, rng: &mut TestRng) -> Lattice {
        let lattice = Lattice::spanned(vectors.iter().map(|v| &v[..]), h11)
            .unwrap()
            .unwrap();
        assert_eq!(lattice.rank(), rank);
        let mut out = vec![0; h11];
        for _ in 0..200 {
            let lambda: Vec<i32> = vectors.iter().map(|_| rng.range(-3, 3) as i32).collect();
            let x: Vec<i32> = (0..h11)
                .map(|t| vectors.iter().zip(&lambda).map(|(v, &l)| v[t] * l).sum())
                .collect();
            let c = lattice.coords(&x).unwrap();
            lattice.to_input(&c, &mut out).unwrap();
            assert_eq!(out, x);
            // A form takes the same value on the class and on its coordinates.
            let w: Vec<i32> = (0..h11).map(|_| rng.range(-5, 5) as i32).collect();
            let wl = lattice.restrict(&w).unwrap();
            let dot = |a: &[i32], b: &[i32]| -> i64 {
                a.iter().zip(b).map(|(&x, &y)| x as i64 * y as i64).sum()
            };
            assert_eq!(dot(&c, &wl), dot(&x, &w));
        }
        lattice
    }

    #[test]
    fn test_full_rank() {
        let v = [vec![1, 0], vec![1, 1], vec![3, 5]];
        assert!(Lattice::spanned(v.iter().map(|v| &v[..]), 2)
            .unwrap()
            .is_none());
        // Full rank but not all of Z^2 is still full rank.
        let v = [vec![2, 0], vec![0, 2]];
        assert!(Lattice::spanned(v.iter().map(|v| &v[..]), 2)
            .unwrap()
            .is_none());
    }

    #[test]
    fn test_sublattices() {
        let mut rng = TestRng::new(3);
        // A line, with a multiple of the primitive vector.
        let l = check(&[vec![2, 4, -6], vec![3, 6, -9]], 3, 1, &mut rng);
        assert_eq!(
            l.basis_vector(0)
                .iter()
                .map(|x| x.abs())
                .collect::<Vec<_>>(),
            [1, 2, 3]
        );
        assert_eq!(l.coords(&[1, 2, 4]), None);
        assert_eq!(l.coords(&[1, 2, -3]).map(|c| c[0].abs()), Some(1));
        // A plane whose lattice is not saturated: (1, 1, 0) is in its span but not in it.
        let l = check(
            &[vec![2, 0, 0], vec![0, 2, 0], vec![1, 1, 0]],
            3,
            2,
            &mut rng,
        );
        assert!(l.coords(&[1, 0, 0]).is_none());
        let l = check(&[vec![2, 0, 0], vec![0, 2, 0]], 3, 2, &mut rng);
        assert!(l.coords(&[1, 1, 0]).is_none());
        // Nothing.
        let l = Lattice::spanned(std::iter::empty(), 4).unwrap().unwrap();
        assert_eq!(l.rank(), 0);
        assert_eq!(l.coords(&[0; 4]), Some(vec![]));
        assert_eq!(l.coords(&[0, 1, 0, 0]), None);
    }

    /// Vectors whose first independent ones span a sublattice of finite index, so that the
    /// others join it, give the lattice they span, no more.
    #[test]
    fn test_joins() {
        let mut rng = TestRng::new(5);
        let l = check(
            &[vec![2, 0, 0], vec![0, 3, 0], vec![1, 3, 0]],
            3,
            2,
            &mut rng,
        );
        assert!(l.coords(&[1, 0, 0]).is_some() && l.coords(&[0, 3, 0]).is_some());
        assert!(l.coords(&[0, 1, 0]).is_none());
        for (h11, rank) in [(5, 2), (10, 4), (60, 12), (200, 30)] {
            let short: Vec<Vec<i32>> = (0..rank)
                .map(|_| (0..h11).map(|_| rng.range(-1, 1) as i32).collect())
                .collect();
            let scaled = |s: &[i32], m: i32| s.iter().map(|&x| m * x).collect::<Vec<_>>();
            let sum = |a: &[i32], b: &[i32]| a.iter().zip(b).map(|(x, y)| x + y).collect();
            // Twice and three times each, and the sums of neighbors and the first, span
            // everything.
            let mut gens: Vec<Vec<i32>> = short.iter().map(|s| scaled(s, 2)).collect();
            gens.extend(short.iter().map(|s| scaled(s, 3)));
            gens.extend((0..rank - 1).map(|i| sum(&short[i], &short[i + 1])));
            gens.push(scaled(&short[0], 5));
            gens.push(scaled(&short[0], 7));
            let lattice = check(&gens, h11, rank, &mut rng);
            let expected = Lattice::spanned(short.iter().map(|v| &v[..]), h11)
                .unwrap()
                .unwrap();
            for j in 0..rank {
                let b: Vec<i32> = lattice.basis_vector(j).iter().map(|&x| x as i32).collect();
                assert!(expected.coords(&b).is_some(), "{h11} {rank}");
            }
        }
    }

    /// Random sublattices given by skewed generators get short bases.
    #[test]
    fn test_reduction() {
        let mut rng = TestRng::new(8);
        for (h11, rank) in [(6, 3), (12, 5), (40, 8), (100, 20), (491, 64)] {
            // Short, nearly orthogonal vectors, then skewed by unimodular row operations.
            let max = |v: &[Vec<i32>]| v.iter().flatten().map(|x| x.abs()).max().unwrap();
            let short: Vec<Vec<i32>> = (0..rank)
                .map(|_| (0..h11).map(|_| rng.range(-1, 1) as i32).collect())
                .collect();
            let mut gens = short.clone();
            while max(&gens) < 50 {
                let (i, j) = (
                    rng.range(0, rank as i64 - 1) as usize,
                    rng.range(0, rank as i64 - 1) as usize,
                );
                if i != j {
                    let m = rng.range(-2, 2) as i32;
                    let gj = gens[j].clone();
                    for (x, y) in gens[i].iter_mut().zip(gj) {
                        *x += m * y;
                    }
                }
            }
            let lattice = check(&gens, h11, rank, &mut rng);
            let reduced = (0..rank)
                .map(|j| {
                    lattice
                        .basis_vector(j)
                        .iter()
                        .map(|&x| x.abs())
                        .max()
                        .unwrap()
                })
                .max()
                .unwrap();
            assert!(reduced <= 2, "{h11} {rank}: {reduced} from {}", max(&gens));
        }
    }
}
