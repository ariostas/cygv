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
use super::table::{SlotMap, SparsePoly};
use malachite::Integer;

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
    /// The curves, numbered in order of insertion.
    curves: SlotMap<K, ()>,
    /// Per prime, `width` residues per curve, zero where a pass found nothing.
    columns: Vec<Vec<u64>>,
}

impl<K: Key> Residues<K> {
    /// An empty set, for `width` values per curve.
    pub fn new(width: usize) -> Self {
        Self {
            width,
            primes: Vec::new(),
            curves: SlotMap::new(1, 1024),
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
        self.curves.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.curves.is_empty()
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
        let n = self.width * self.curves.len();
        self.columns.extend((0..NL).map(|_| vec![0; n]));
        let first = self.columns.len() - NL;
        for i in 0..values.len() {
            let (s, created) = self.curves.insert(values.keys()[i], values.degs()[i]);
            if created {
                for col in &mut self.columns {
                    col.resize(col.len() + self.width, 0);
                }
            }
            let at = s as usize * self.width;
            for (k, &v) in values.values(i).iter().enumerate() {
                for (l, r) in ring.residues(v).into_iter().enumerate() {
                    self.columns[first + l][at + k] = r;
                }
            }
        }
    }

    /// Whether the lift from all primes but the last `check` agrees with the residues modulo
    /// those at every curve. False unless there are more than `check` primes.
    pub fn is_stable(&self, check: usize) -> bool {
        let n = self.primes.len();
        if n <= check {
            return false;
        }
        let lift = Lift::new(&self.primes);
        let mut digits = vec![0; n - check];
        (0..self.curves.len() * self.width).all(|i| {
            lift.digits(|j| self.columns[j][i], &mut digits);
            (n - check..n).all(|j| lift.residue(&digits, self.primes[j]) == self.columns[j][i])
        })
    }

    /// The lifts from all primes, as (curve, degree, `width` integers) in increasing order
    /// of degree and then key, leaving out the curves where all of them are zero.
    pub fn lift(&self) -> Vec<(K, u32, Vec<Integer>)> {
        let lift = Lift::new(&self.primes);
        let mut digits = vec![0; self.primes.len()];
        let mut order: Vec<u32> = (0..self.curves.len() as u32).collect();
        order.sort_unstable_by_key(|&s| (self.curves.deg(s), self.curves.key(s)));
        order
            .into_iter()
            .filter_map(|s| {
                let at = s as usize * self.width;
                let values: Vec<Integer> = (at..at + self.width)
                    .map(|i| {
                        lift.digits(|j| self.columns[j][i], &mut digits);
                        lift.integer(&digits)
                    })
                    .collect();
                (values.iter().any(|v| *v != 0))
                    .then(|| (self.curves.key(s), self.curves.deg(s), values))
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modular::field::primes;
    use crate::modular::key::PackedKey;
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
}
