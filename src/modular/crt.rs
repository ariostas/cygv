//! Lifting residues to integers by the Chinese remainder theorem.
//!
//! The passes of the modular pipeline each give the invariants modulo a few primes.
//! [`Residues`] collects them, prime by prime, and lifts each to the integer of least
//! absolute value with those residues. The lift from all primes but the last `check` ones
//! is accepted when it agrees with the residues modulo those as well, at every curve: a lift
//! that was cut short by too few primes would have to agree with each check prime by
//! accident.
//!
//! The lift goes through Garner's mixed-radix digits, $x = v_0 + v_1 p_0 + v_2 p_0 p_1 +
//! \dots$ with $0 \le v_i < p_i$. Whether the least absolute value is negative, and the
//! residue of the lift modulo another prime, follow from the digits in word-sized
//! arithmetic, so only the final integers need a big-number type.
//!
//! The residue sets do not depend on how the primes were split into passes, which is what
//! lets a run be split across machines (`PLAN.md`, section 4.1).

use super::field::{LaneField, Lanes};
use super::key::Key;
use super::table::SparsePoly;
use malachite::base::num::arithmetic::traits::{FloorSqrt, Gcd, UnsignedAbs};
use malachite::{Integer, Rational};
use rayon::prelude::*;

/// The fewest primes [`Residues::fraction`] reconstructs from: 248 bits, enough for a
/// numerator and a denominator of 124 bits each.
const FRACTION_PRIMES: usize = 4;

/// The fraction $a / b$ with $|a|, b \le$ `bound` and $a \equiv b x \pmod m$, for $0 \le x <
/// m$, if there is one in lowest terms (Wang's rational reconstruction: the extended
/// Euclidean algorithm on $m$ and $x$, stopped at the first remainder below the bound).
fn reconstruct(x: &Integer, m: &Integer, bound: &Integer) -> Option<(Integer, Integer)> {
    let (mut r0, mut r1) = (m.clone(), x.clone());
    let (mut t0, mut t1) = (Integer::from(0), Integer::from(1));
    while &r1 > bound {
        let q = &r0 / &r1;
        let r2 = &r0 - &q * &r1;
        (r0, r1) = (r1, r2);
        let t2 = &t0 - &q * &t1;
        (t0, t1) = (t1, t2);
    }
    let (a, b) = if t1 < 0 { (-r1, -t1) } else { (r1, t1) };
    let coprime = (&a).unsigned_abs().gcd((&b).unsigned_abs()) == 1u32;
    (b != 0 && &b <= bound && coprime).then_some((a, b))
}

/// $a b \bmod p$
#[inline]
fn mulmod(a: u64, b: u64, p: u64) -> u64 {
    (a as u128 * b as u128 % p as u128) as u64
}

/// $a^{-1} \bmod p$ for a prime $p$ not dividing $a$.
fn invmod(a: u64, p: u64) -> u64 {
    let (mut x, mut e, mut r) = (a % p, p - 2, 1);
    while e != 0 {
        if e & 1 == 1 {
            r = mulmod(r, x, p);
        }
        x = mulmod(x, x, p);
        e >>= 1;
    }
    r
}

/// The precomputed constants for lifting from the first `n` of a list of distinct primes.
struct Lift<'a> {
    primes: &'a [u64],
    /// `inv[i][j]` $= p_j^{-1} \bmod p_i$ for $j < i$.
    inv: Vec<Vec<u64>>,
}

impl<'a> Lift<'a> {
    fn new(primes: &'a [u64]) -> Self {
        let inv = (0..primes.len())
            .map(|i| (0..i).map(|j| invmod(primes[j], primes[i])).collect())
            .collect();
        Self { primes, inv }
    }

    /// The mixed-radix digits of the residues `r` modulo the first `digits.len()` primes.
    fn digits(&self, r: impl Fn(usize) -> u64, digits: &mut [u64]) {
        for i in 0..digits.len() {
            let p = self.primes[i];
            let mut x = r(i) % p;
            for (&v, &inv) in digits[..i].iter().zip(&self.inv[i]) {
                x = mulmod((x + p - v % p) % p, inv, p);
            }
            digits[i] = x;
        }
    }

    /// Whether the lift with these digits is negative, i.e. above half the product of the
    /// primes. Half the product, rounded down, has the digits $(p_i - 1) / 2$.
    fn is_negative(&self, digits: &[u64]) -> bool {
        for i in (0..digits.len()).rev() {
            let half = (self.primes[i] - 1) / 2;
            if digits[i] != half {
                return digits[i] > half;
            }
        }
        false
    }

    /// The residue of the lift with these digits modulo the prime `q`.
    fn residue(&self, digits: &[u64], q: u64) -> u64 {
        // Horner's rule from the top digit, and the product of the primes along the way.
        let mut x = 0;
        let mut m = 1;
        for i in (0..digits.len()).rev() {
            x = (mulmod(x, self.primes[i] % q, q) + digits[i] % q) % q;
            m = mulmod(m, self.primes[i] % q, q);
        }
        if self.is_negative(digits) {
            (x + q - m) % q
        } else {
            x
        }
    }

    /// The lift with these digits.
    fn integer(&self, digits: &[u64]) -> Integer {
        let mut x = Integer::from(0);
        let mut m = Integer::from(1);
        for i in (0..digits.len()).rev() {
            x = x * Integer::from(self.primes[i]) + Integer::from(digits[i]);
            m *= Integer::from(self.primes[i]);
        }
        if self.is_negative(digits) {
            x - m
        } else {
            x
        }
    }
}

/// The residues of the invariants of a run, at every curve that is nonzero modulo some prime,
/// with `width` values per curve.
#[derive(Debug, Clone)]
pub struct Residues<K> {
    width: usize,
    primes: Vec<u64>,
    /// The curves, in increasing order of degree and then key, as the passes give them.
    keys: Vec<K>,
    degs: Vec<u32>,
    /// Per prime, `width` residues per curve, zero where a pass found nothing.
    columns: Vec<Vec<u64>>,
}

impl<K: Key> Residues<K> {
    /// An empty set, for `width` values per curve.
    pub fn new(width: usize) -> Self {
        Self {
            width,
            primes: Vec::new(),
            keys: Vec::new(),
            degs: Vec::new(),
            columns: Vec::new(),
        }
    }

    /// The primes so far.
    #[inline]
    pub fn primes(&self) -> &[u64] {
        &self.primes
    }

    /// The number of curves so far.
    #[inline]
    pub fn len(&self) -> usize {
        self.keys.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    /// Adds the result of a pass, i.e. the invariants modulo the primes of `ring`.
    ///
    /// # Panics
    ///
    /// Panics if the width does not match or one of the primes was added before.
    pub fn add_pass<const NL: usize>(
        &mut self,
        ring: &LaneField<NL>,
        values: &SparsePoly<K, Lanes<NL>>,
    ) {
        assert_eq!(values.width(), self.width, "the widths do not match");
        for f in ring.fields() {
            assert!(
                !self.primes.contains(&f.p()),
                "the prime {} was added twice",
                f.p()
            );
            self.primes.push(f.p());
        }
        // Both lists are in order of degree and then key: the index of each value among the
        // curves after the merge, which only moves the old curves if the pass has new ones.
        let order = |a: (u32, K), b: (u32, K)| a.cmp(&b);
        let mut at = Vec::with_capacity(values.len());
        let mut new = 0;
        let mut i = 0;
        for v in 0..values.len() {
            let x = (values.degs()[v], values.keys()[v]);
            while i < self.keys.len() && order((self.degs[i], self.keys[i]), x).is_lt() {
                i += 1;
            }
            if i < self.keys.len() && (self.degs[i], self.keys[i]) == x {
                at.push(i + new);
            } else {
                at.push(i + new);
                new += 1;
            }
        }
        if new > 0 {
            let n = self.keys.len() + new;
            let (mut keys, mut degs) = (Vec::with_capacity(n), Vec::with_capacity(n));
            let mut moved = Vec::with_capacity(self.keys.len());
            let mut i = 0;
            for v in 0..values.len() {
                let x = (values.degs()[v], values.keys()[v]);
                while i < self.keys.len() && order((self.degs[i], self.keys[i]), x).is_lt() {
                    moved.push(keys.len());
                    keys.push(self.keys[i]);
                    degs.push(self.degs[i]);
                    i += 1;
                }
                if i < self.keys.len() && (self.degs[i], self.keys[i]) == x {
                    moved.push(keys.len());
                    i += 1;
                }
                keys.push(x.1);
                degs.push(x.0);
            }
            for i in i..self.keys.len() {
                moved.push(keys.len());
                keys.push(self.keys[i]);
                degs.push(self.degs[i]);
            }
            let w = self.width;
            for col in &mut self.columns {
                let mut c = vec![0; n * w];
                for (old, &to) in moved.iter().enumerate() {
                    c[to * w..(to + 1) * w].copy_from_slice(&col[old * w..(old + 1) * w]);
                }
                *col = c;
            }
            (self.keys, self.degs) = (keys, degs);
        }
        let n = self.width * self.keys.len();
        self.columns.extend((0..NL).map(|_| vec![0; n]));
        let first = self.columns.len() - NL;
        for (v, &s) in at.iter().enumerate() {
            let at = s * self.width;
            for (k, &x) in values.values(v).iter().enumerate() {
                for (l, r) in ring.residues(x).into_iter().enumerate() {
                    self.columns[first + l][at + k] = r;
                }
            }
        }
    }

    /// Whether the lift from all primes but the last `check` agrees with the residues modulo
    /// those at every curve. False unless there are more than `check` primes.
    pub fn is_stable(&self, check: usize) -> bool {
        self.unstable(check, 1).is_some_and(|u| u.is_empty())
    }

    /// The first `limit` values, in order of degree, whose lift from all primes but the last
    /// `check` disagrees with the residues modulo those, as indices of the columns, or `None`
    /// unless there are more than `check` primes.
    pub fn unstable(&self, check: usize, limit: usize) -> Option<Vec<usize>> {
        let n = self.primes.len();
        if n <= check {
            return None;
        }
        let lift = Lift::new(&self.primes);
        let len = self.keys.len() * self.width;
        let mut res = Vec::new();
        // In pieces, so that a lift far from stable is not scanned to the end.
        let piece = (len / 16).max(1 << 16);
        for start in (0..len).step_by(piece) {
            let found: Vec<usize> = (start..len.min(start + piece))
                .into_par_iter()
                .map_init(
                    || vec![0; n - check],
                    |digits, i| {
                        lift.digits(|j| self.columns[j][i], digits);
                        let agrees = (n - check..n)
                            .all(|j| lift.residue(digits, self.primes[j]) == self.columns[j][i]);
                        (!agrees).then_some(i)
                    },
                )
                .flatten()
                .collect();
            res.extend(found);
            if res.len() >= limit {
                res.truncate(limit);
                break;
            }
        }
        Some(res)
    }

    /// A value that is a fraction rather than an integer, with its degree: one of `unstable`,
    /// values whose lift from all primes but the last `check` disagrees with those (from
    /// [`Residues::unstable`]), that is a fraction $a / b$ with $b > 1$ modulo all but the last
    /// `check` primes, by rational reconstruction, and agrees with the residues modulo
    /// those. An integer that is only too large for the primes so far passes that test by
    /// accident with a chance of about one in a prime, as for [`Residues::is_stable`]. So
    /// invariants that are not integers, which is what inconsistent input gives, show after
    /// a few primes rather than after all of those allowed.
    pub fn fraction(&self, check: usize, unstable: &[usize]) -> Option<(u32, Rational)> {
        let n = self.primes.len();
        if n < check + FRACTION_PRIMES {
            return None;
        }
        let lift = Lift::new(&self.primes);
        let m = self.primes[..n - check]
            .iter()
            .fold(Integer::from(1), |m, &p| m * Integer::from(p));
        let bound = Integer::from((&m >> 1u32).unsigned_abs().floor_sqrt());
        let mut digits = vec![0; n - check];
        for &i in unstable {
            lift.digits(|j| self.columns[j][i], &mut digits);
            let mut x = lift.integer(&digits);
            if x < 0 {
                x += &m;
            }
            let Some((a, b)) = reconstruct(&x, &m, &bound) else {
                continue;
            };
            let agrees = (n - check..n).all(|j| {
                let p = Integer::from(self.primes[j]);
                let r = Integer::from(self.columns[j][i]);
                ((&a - &b * r) % p) == 0
            });
            if b > 1 && agrees {
                return Some((
                    self.degs[i / self.width],
                    Rational::from(a) / Rational::from(b),
                ));
            }
        }
        None
    }

    /// The lifts from all primes, as (curve, degree, `width` integers) in increasing order
    /// of degree and then key, leaving out the curves where all of them are zero.
    pub fn lift(&self) -> Vec<(K, u32, Vec<Integer>)> {
        let lift = Lift::new(&self.primes);
        // Collected in place, one entry per curve, and then filtered in place.
        let lifted: Vec<Option<(K, u32, Vec<Integer>)>> = (0..self.keys.len())
            .into_par_iter()
            .map_init(
                || vec![0; self.primes.len()],
                |digits, s| {
                    let at = s * self.width;
                    let values: Vec<Integer> = (at..at + self.width)
                        .map(|i| {
                            lift.digits(|j| self.columns[j][i], digits);
                            lift.integer(digits)
                        })
                        .collect();
                    (values.iter().any(|v| *v != 0)).then(|| (self.keys[s], self.degs[s], values))
                },
            )
            .collect();
        // Mapped in place, without a second buffer.
        let mut lifted = lifted;
        lifted.retain(Option::is_some);
        lifted.into_iter().map(Option::unwrap).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modular::field::primes;
    use crate::modular::key::PackedKey;
    use crate::modular::table::SlotMap;
    use crate::modular::TestRng;
    use malachite::base::num::arithmetic::traits::Mod;
    use malachite::base::num::basic::traits::{NegativeOne, One, Zero};

    #[test]
    fn test_lift() {
        let ps: Vec<u64> = primes().take(5).collect();
        let lift = Lift::new(&ps);
        let mut rng = TestRng::new(31);
        for n in 1..=5 {
            let m: Integer = ps[..n].iter().map(|&p| Integer::from(p)).product();
            let half: Integer = (m - Integer::ONE) / Integer::from(2);
            let mut values = vec![
                Integer::ZERO,
                half.clone(),
                -half.clone(),
                Integer::NEGATIVE_ONE,
            ];
            for _ in 0..20 {
                // Random values of up to the full size, of either sign.
                let mut x = Integer::ZERO;
                for _ in 0..n {
                    x = (x << 64u32) + Integer::from(rng.next_u64());
                }
                x %= &half;
                if rng.next_u64().is_multiple_of(2) {
                    x = -x;
                }
                values.push(x);
            }
            let mut digits = vec![0; n];
            for x in values {
                let res = |j: usize| u64::try_from(&(&x).mod_op(Integer::from(ps[j]))).unwrap();
                lift.digits(res, &mut digits);
                assert_eq!(lift.integer(&digits), x);
                for (j, &q) in ps.iter().enumerate() {
                    assert_eq!(lift.residue(&digits, q), res(j), "{x} mod {q}");
                }
            }
        }
    }

    /// The residues of `values` (per key, two components) as a pass over two primes.
    fn pass(lf: &LaneField<2>, values: &[(u128, [i128; 2])]) -> SparsePoly<PackedKey, Lanes<2>> {
        let mut map = SlotMap::new(2, 4);
        for &(k, v) in values {
            let slot = map.get_or_insert(PackedKey(k), k as u32);
            for (x, &y) in slot.iter_mut().zip(&v) {
                *x = Lanes(core::array::from_fn(|l| {
                    let f = &lf.fields()[l];
                    f.encode(y.rem_euclid(f.p() as i128) as u64)
                }));
            }
        }
        map.take()
    }

    #[test]
    fn test_residues() {
        let ps: Vec<u64> = primes().take(6).collect();
        let lanes: Vec<LaneField<2>> = ps.chunks(2).map(|c| LaneField::new([c[0], c[1]])).collect();
        let int = |v: &[i128]| v.iter().map(|&x| Integer::from(x)).collect::<Vec<_>>();

        // Small values: one pass is enough, with a lift from one prime checked by the other.
        let small = [(1, [5, -3]), (2, [0, 7]), (3, [-1, 0])];
        let mut res = Residues::new(2);
        assert!(!res.is_stable(1));
        res.add_pass(&lanes[0], &pass(&lanes[0], &small));
        assert!(res.is_stable(1));
        assert!(!res.is_stable(2));
        let expected = vec![
            (PackedKey(1), 1, int(&[5, -3])),
            (PackedKey(2), 2, int(&[0, 7])),
            (PackedKey(3), 3, int(&[-1, 0])),
        ];
        assert_eq!(res.lift(), expected);

        // A value beyond one prime fails the check with two primes. A multiple of the first
        // two primes vanishes in the first pass, and is only found by the second; it needs
        // a lift from three primes.
        let big = -(i64::MAX as i128) + 5;
        let p01 = ps[0] as i128 * ps[1] as i128;
        let values = [(1, [big, 1]), (4, [0, p01])];
        let mut res = Residues::new(2);
        res.add_pass(&lanes[0], &pass(&lanes[0], &values));
        assert_eq!(res.len(), 1);
        assert!(!res.is_stable(1));
        res.add_pass(&lanes[1], &pass(&lanes[1], &values));
        assert_eq!(res.len(), 2);
        assert!(res.is_stable(1));
        assert!(!res.is_stable(2));
        res.add_pass(&lanes[2], &pass(&lanes[2], &values));
        assert!(res.is_stable(3));
        let expected = vec![
            (PackedKey(1), 1, int(&[big, 1])),
            (PackedKey(4), 4, int(&[0, p01])),
        ];
        assert_eq!(res.lift(), expected);
    }

    #[test]
    fn test_fraction() {
        let ps: Vec<u64> = primes().take(6).collect();
        let lanes: Vec<LaneField<2>> = ps.chunks(2).map(|c| LaneField::new([c[0], c[1]])).collect();
        // -7/3 modulo each prime, at curve 2, beside integers, one of which needs three
        // primes.
        let pass = |lf: &LaneField<2>| {
            let mut map = SlotMap::new(1, 4);
            let fraction = Lanes(core::array::from_fn(|l| {
                let f = &lf.fields()[l];
                f.mul(f.encode(f.p() - 7), f.inv(f.encode(3)).unwrap())
            }));
            let int = |y: i128| {
                Lanes(core::array::from_fn(|l| {
                    let f = &lf.fields()[l];
                    f.encode(y.rem_euclid(f.p() as i128) as u64)
                }))
            };
            map.get_or_insert(PackedKey(1), 1)[0] = int(5);
            map.get_or_insert(PackedKey(2), 2)[0] = fraction;
            map.get_or_insert(PackedKey(3), 3)[0] = int(i128::MAX / 3);
            map.take()
        };
        let mut res = Residues::new(1);
        res.add_pass(&lanes[0], &pass(&lanes[0]));
        res.add_pass(&lanes[1], &pass(&lanes[1]));
        // Four primes are too few, with one of them to check.
        assert!(!res.is_stable(1));
        assert_eq!(res.fraction(1, &res.unstable(1, 64).unwrap()), None);
        res.add_pass(&lanes[2], &pass(&lanes[2]));
        // Only the fraction is unstable, the integer fitting five primes.
        let unstable = res.unstable(1, 64).unwrap();
        assert_eq!(unstable, [1]);
        assert_eq!(
            res.fraction(1, &unstable),
            Some((2, Rational::from_signeds(-7, 3)))
        );
        // Only the values given are tried.
        assert_eq!(res.fraction(1, &[]), None);
        assert_eq!(res.fraction(1, &[0, 2]), None);
        assert_eq!(res.unstable(6, 1), None);
    }
}
