//! Arithmetic modulo word-sized primes.
//!
//! A pass of the modular pipeline does all of its arithmetic modulo a few primes just below
//! $2^{62}$. [`Field`] is the arithmetic modulo one of them, in Montgomery form. [`Lanes`]
//! holds the residues of one number modulo several primes side by side, and [`LaneField`]
//! computes with them, one prime per lane. The primes are drawn from the fixed sequence of
//! [`primes`].
//!
//! The stages are written against the [`Ring`] trait rather than against [`LaneField`], so that
//! the certifier can run them over upward-rounded bounds instead.

use core::fmt::Debug;

/// The number system the stages compute in.
///
/// An implementation is a context, e.g. the primes for [`LaneField`], and its elements are
/// plain `Copy` values: a residue does not know its modulus, so every operation goes through
/// the context. [`Default`] must give zero, which lets tables hold elements without access to
/// the context.
///
/// This is deliberately not called a field. Its intended second implementation, the
/// certifier's upward-rounded majorants, replaces subtraction by addition and negation by the
/// identity, and has inverses only for the exact values it is ever asked to invert.
pub trait Ring: Sync {
    /// The element type.
    type Elem: Copy + Default + PartialEq + Send + Sync + Debug;

    /// Whether the arithmetic is exact, so that a sum does not depend on the order of its
    /// terms. When it rounds, the stages take every sum in an order that depends neither on
    /// the scheduling nor on the number of threads, at some cost in time and memory, so that
    /// the results are reproducible.
    const EXACT: bool = true;

    /// The multiplicative identity.
    fn one(&self) -> Self::Elem;

    /// The image of an integer.
    fn int(&self, v: i64) -> Self::Elem;

    fn add(&self, a: Self::Elem, b: Self::Elem) -> Self::Elem;

    fn sub(&self, a: Self::Elem, b: Self::Elem) -> Self::Elem;

    fn neg(&self, a: Self::Elem) -> Self::Elem;

    fn mul(&self, a: Self::Elem, b: Self::Elem) -> Self::Elem;

    /// Computes `acc += a * b`.
    #[inline]
    fn fma(&self, acc: &mut Self::Elem, a: Self::Elem, b: Self::Elem) {
        *acc = self.add(*acc, self.mul(a, b));
    }

    /// The multiplicative inverse, or `None` if `a` has none.
    ///
    /// It is only ever applied to exact values: small integers, and the constant term of a
    /// series.
    fn inv(&self, a: Self::Elem) -> Option<Self::Elem>;

    #[inline]
    fn is_zero(&self, a: Self::Elem) -> bool {
        a == Self::Elem::default()
    }

    /// The inverses of $1, \dots, n - 1$, after a zero in place of the inverse of $0$, or
    /// `None` if one of them has no inverse.
    ///
    /// The period tables are built from these. The default takes a single inverse, of
    /// $(n - 1)!$, and recovers the others by multiplication, which is exact in a field; an
    /// implementation that rounds overrides it to invert each integer by itself.
    fn inverses(&self, n: usize) -> Option<Vec<Self::Elem>> {
        let ints: Vec<_> = (0..n).map(|i| self.int(i as i64)).collect();
        // prefix[i] = i!
        let mut prefix = vec![self.one(); n.max(1)];
        for i in 1..n {
            prefix[i] = self.mul(prefix[i - 1], ints[i]);
        }
        let mut inv = vec![Self::Elem::default(); n];
        if n > 1 {
            // Runs from 1/(i!) down to 1/(1!), taking 1/i = (i - 1)! / i! along the way.
            let mut inv_fact = self.inv(prefix[n - 1])?;
            for i in (1..n).rev() {
                inv[i] = self.mul(inv_fact, prefix[i - 1]);
                inv_fact = self.mul(inv_fact, ints[i]);
            }
        }
        Some(inv)
    }
}

/// Arithmetic modulo an odd prime $p < 2^{62}$, in Montgomery form.
///
/// A residue $a$ is stored as $aR \bmod p$ with $R = 2^{64}$, which turns the reduction after
/// a multiplication into two multiplications and a shift ([`Field::reduce`]). The values a
/// `Field` takes and returns are in that form, except for [`Field::encode`] and
/// [`Field::int`], which convert into it, and [`Field::decode`], which converts out of it.
/// Zero is zero in both forms. The bound on $p$ keeps every intermediate sum below $2^{64}$
/// and every product plus correction below $2^{128}$.
///
/// The layout is fixed (`repr(C)`) so that the constants can be handed to other
/// implementations of the same arithmetic, such as a future GPU kernel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(C)]
pub struct Field {
    /// The prime.
    p: u64,
    /// $-p^{-1} \bmod 2^{64}$.
    pinv: u64,
    /// $R^2 \bmod p$.
    r2: u64,
}

impl Field {
    /// Sets up the arithmetic modulo `p`.
    ///
    /// # Panics
    ///
    /// Panics if `p` is not an odd prime below $2^{62}$.
    pub fn new(p: u64) -> Self {
        assert!(
            p > 2 && p < 1 << 62 && is_prime(p),
            "{p} is not an odd prime below 2^62"
        );
        // Newton's iteration for p^{-1} mod 2^64; each step doubles the number of correct
        // bits, starting from the 3 that inv = p always has.
        let mut inv = p;
        for _ in 0..5 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(p.wrapping_mul(inv)));
        }
        let r = (1u128 << 64) % p as u128;
        Self {
            p,
            pinv: inv.wrapping_neg(),
            r2: (r * r % p as u128) as u64,
        }
    }

    /// The prime.
    #[inline]
    pub fn p(&self) -> u64 {
        self.p
    }

    /// Montgomery reduction: $t R^{-1} \bmod p$, for $t < p^2$.
    #[inline]
    pub fn reduce(&self, t: u128) -> u64 {
        let m = (t as u64).wrapping_mul(self.pinv);
        let u = ((t + m as u128 * self.p as u128) >> 64) as u64;
        if u >= self.p {
            u - self.p
        } else {
            u
        }
    }

    #[inline]
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        self.reduce(a as u128 * b as u128)
    }

    #[inline]
    pub fn add(&self, a: u64, b: u64) -> u64 {
        let c = a + b;
        if c >= self.p {
            c - self.p
        } else {
            c
        }
    }

    #[inline]
    pub fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b {
            a - b
        } else {
            a + self.p - b
        }
    }

    #[inline]
    pub fn neg(&self, a: u64) -> u64 {
        if a == 0 {
            0
        } else {
            self.p - a
        }
    }

    /// Converts an integer into Montgomery form.
    #[inline]
    pub fn encode(&self, x: u64) -> u64 {
        self.mul(x % self.p, self.r2)
    }

    /// Converts a signed integer into Montgomery form.
    #[inline]
    pub fn int(&self, v: i64) -> u64 {
        // p < 2^62 fits an i64, and the remainder is non-negative.
        self.encode(v.rem_euclid(self.p as i64) as u64)
    }

    /// Converts out of Montgomery form, giving the residue in `0..p`.
    #[inline]
    pub fn decode(&self, a: u64) -> u64 {
        self.reduce(a as u128)
    }

    #[inline]
    pub fn one(&self) -> u64 {
        self.encode(1)
    }

    /// Computes $a^e$.
    pub fn pow(&self, mut a: u64, mut e: u64) -> u64 {
        let mut r = self.one();
        while e != 0 {
            if e & 1 == 1 {
                r = self.mul(r, a);
            }
            a = self.mul(a, a);
            e >>= 1;
        }
        r
    }

    /// The inverse of `a`, or `None` if `a` is zero.
    pub fn inv(&self, a: u64) -> Option<u64> {
        // Fermat's little theorem. Inverses are rare enough (tables and constant terms) that
        // an extended Euclid would not be worth the extra code.
        (a != 0).then(|| self.pow(a, self.p - 2))
    }
}

/// Whether `n` is prime, by a deterministic Miller–Rabin test.
///
/// The first twelve primes as bases are enough for every `n` below $3.3 \cdot 10^{24}$, so in
/// particular for every `u64`.
pub fn is_prime(n: u64) -> bool {
    const BASES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    let mulmod = |a: u64, b: u64| (a as u128 * b as u128 % n as u128) as u64;
    if n < 2 {
        return false;
    }
    for p in BASES {
        if n.is_multiple_of(p) {
            return n == p;
        }
    }
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;
    'bases: for a in BASES {
        let mut x = 1;
        let (mut b, mut e) = (a, d);
        while e != 0 {
            if e & 1 == 1 {
                x = mulmod(x, b);
            }
            b = mulmod(b, b);
            e >>= 1;
        }
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulmod(x, x);
            if x == n - 1 {
                continue 'bases;
            }
        }
        return false;
    }
    true
}

/// The primes that passes take their lanes from: the primes below $2^{62}$, counting down.
///
/// The sequence is fixed, so a pass is identified by the indices of its primes in it. It is
/// the same sequence as `nth_prime` in cgv (`cgv/gv.c`), so residues
/// can be compared with it lane for lane.
pub fn primes() -> impl Iterator<Item = u64> {
    core::iter::successors(Some((1u64 << 62) - 1), |n| Some(n - 2)).filter(|&n| is_prime(n))
}

/// The residues of one number modulo the primes of a [`LaneField`], in Montgomery form.
///
/// A plain array with a fixed layout, so that it can be handed to other implementations of
/// the same arithmetic as is. The constants it is relative to travel separately, in the
/// [`LaneField`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(transparent)]
pub struct Lanes<const NL: usize>(pub [u64; NL]);

impl<const NL: usize> Default for Lanes<NL> {
    #[inline]
    fn default() -> Self {
        Self([0; NL])
    }
}

/// Arithmetic modulo `NL` primes at once, one per lane.
///
/// A number is known modulo the product of the primes, which is what lets a single pass
/// determine integers up to about $2^{62 \cdot NL}$ in absolute value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(transparent)]
pub struct LaneField<const NL: usize> {
    fields: [Field; NL],
}

impl<const NL: usize> LaneField<NL> {
    /// Sets up the arithmetic modulo the given primes.
    ///
    /// # Panics
    ///
    /// Panics if a prime is not an odd prime below $2^{62}$, or if two of them are equal.
    pub fn new(primes: [u64; NL]) -> Self {
        for (i, p) in primes.iter().enumerate() {
            assert!(!primes[..i].contains(p), "the prime {p} is repeated");
        }
        Self {
            fields: primes.map(Field::new),
        }
    }

    /// The arithmetic of each lane.
    #[inline]
    pub fn fields(&self) -> &[Field; NL] {
        &self.fields
    }

    /// The residue in each lane, out of Montgomery form.
    #[inline]
    pub fn residues(&self, a: Lanes<NL>) -> [u64; NL] {
        core::array::from_fn(|l| self.fields[l].decode(a.0[l]))
    }
}

impl<const NL: usize> Ring for LaneField<NL> {
    type Elem = Lanes<NL>;

    #[inline]
    fn one(&self) -> Lanes<NL> {
        Lanes(core::array::from_fn(|l| self.fields[l].one()))
    }

    #[inline]
    fn int(&self, v: i64) -> Lanes<NL> {
        Lanes(core::array::from_fn(|l| self.fields[l].int(v)))
    }

    #[inline]
    fn add(&self, a: Lanes<NL>, b: Lanes<NL>) -> Lanes<NL> {
        Lanes(core::array::from_fn(|l| self.fields[l].add(a.0[l], b.0[l])))
    }

    #[inline]
    fn sub(&self, a: Lanes<NL>, b: Lanes<NL>) -> Lanes<NL> {
        Lanes(core::array::from_fn(|l| self.fields[l].sub(a.0[l], b.0[l])))
    }

    #[inline]
    fn neg(&self, a: Lanes<NL>) -> Lanes<NL> {
        Lanes(core::array::from_fn(|l| self.fields[l].neg(a.0[l])))
    }

    #[inline]
    fn mul(&self, a: Lanes<NL>, b: Lanes<NL>) -> Lanes<NL> {
        Lanes(core::array::from_fn(|l| self.fields[l].mul(a.0[l], b.0[l])))
    }

    /// The inverse in every lane, or `None` if the value is zero in any of them.
    fn inv(&self, a: Lanes<NL>) -> Option<Lanes<NL>> {
        let mut r = Lanes::default();
        for l in 0..NL {
            r.0[l] = self.fields[l].inv(a.0[l])?;
        }
        Some(r)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modular::TestRng;
    use malachite::base::num::arithmetic::traits::{Mod, ModInverse, ModPow};
    use malachite::base::num::factorization::traits::IsPrime;
    use malachite::{Integer, Natural};

    /// Primes of every size the arithmetic has to handle, from tiny to the largest allowed.
    fn sample_primes() -> Vec<u64> {
        let mut ps = vec![3, 5, 7, 65_537, 2_147_483_647, 4_294_967_311];
        ps.extend(primes().take(3));
        ps
    }

    #[test]
    fn test_constants() {
        for p in sample_primes() {
            let f = Field::new(p);
            assert_eq!(f.pinv.wrapping_mul(p), u64::MAX, "pinv for {p}");
            let r2 = (Natural::from(1u32) << 128u32) % Natural::from(p);
            assert_eq!(Natural::from(f.r2), r2, "r2 for {p}");
        }
    }

    #[test]
    fn test_operations() {
        let mut rng = TestRng::new(1);
        for p in sample_primes() {
            let f = Field::new(p);
            let big_p = Integer::from(p);
            let modp = |x: Integer| u64::try_from(&x.mod_op(&big_p)).unwrap();
            for _ in 0..2000 {
                // Signed inputs, including ones beyond p, reduced on the way in.
                let x = rng.next_u64() as i64;
                let y = rng.range(-1000, 1000);
                let (a, b) = (f.int(x), f.int(y));
                let (bx, by) = (Integer::from(x), Integer::from(y));
                assert_eq!(f.decode(a), modp(bx.clone()));
                assert_eq!(f.decode(f.add(a, b)), modp(bx.clone() + &by));
                assert_eq!(f.decode(f.sub(a, b)), modp(bx.clone() - &by));
                assert_eq!(f.decode(f.sub(b, a)), modp(by.clone() - &bx));
                assert_eq!(f.decode(f.neg(a)), modp(-bx.clone()));
                assert_eq!(f.decode(f.mul(a, b)), modp(bx.clone() * &by));
                let e = rng.next_u64() % 1000;
                let x = modp(bx.clone());
                assert_eq!(f.decode(f.pow(a, e)), x.mod_pow(e, p));
                match x {
                    0 => assert_eq!(f.inv(a), None),
                    _ => assert_eq!(f.decode(f.inv(a).unwrap()), x.mod_inverse(p).unwrap()),
                }
                let u = rng.next_u64();
                assert_eq!(f.decode(f.encode(u)), u % p);
            }
            assert_eq!(f.decode(f.one()), 1);
            assert_eq!(f.int(0), 0);
            assert_eq!(f.int(p as i64), 0);
            assert_eq!(f.decode(f.int(-1)), p - 1);
            assert_eq!(f.decode(f.int(i64::MIN)), modp(Integer::from(i64::MIN)));
            assert_eq!(f.inv(0), None);
            // The largest residues, where the conditional subtractions matter.
            let m = f.encode(p - 1);
            assert_eq!(f.decode(f.add(m, m)), p - 2);
            assert_eq!(f.decode(f.mul(m, m)), 1);
            assert_eq!(f.sub(0, 0), 0);
        }
    }

    #[test]
    #[should_panic(expected = "is not an odd prime")]
    fn test_field_rejects_large_prime() {
        let p = ((1u64 << 62) + 1..).find(|n| n.is_prime()).unwrap();
        Field::new(p);
    }

    #[test]
    #[should_panic(expected = "is not an odd prime")]
    fn test_field_rejects_composite() {
        Field::new(91);
    }

    #[test]
    fn test_is_prime() {
        let agrees = |n: u64| is_prime(n) == n.is_prime();
        for n in 0..5000 {
            assert!(agrees(n), "{n}");
        }
        let mut rng = TestRng::new(2);
        for _ in 0..5000 {
            let n = rng.next_u64() | 1;
            assert!(agrees(n), "{n}");
            assert!(agrees(n >> 2), "{}", n >> 2);
        }
        // Carmichael numbers and strong pseudoprimes to the first few prime bases.
        for n in [561, 41_041, 3_215_031_751, 3_825_123_056_546_413_051] {
            assert!(agrees(n), "{n}");
        }
        assert!(is_prime(u64::MAX - 58));
        assert!(!is_prime(u64::MAX));
    }

    #[test]
    fn test_primes() {
        let expected = (0..1u64 << 62).rev().filter(|n| n.is_prime());
        assert!(primes().take(20).eq(expected.take(20)));
    }

    #[test]
    fn test_lanes() {
        let ps: [u64; 3] = [7, primes().next().unwrap(), 65_537];
        let lf = LaneField::new(ps);
        let mut rng = TestRng::new(3);
        for _ in 0..1000 {
            let (x, y) = (rng.next_u64() as i64, rng.range(-100, 100));
            let (a, b) = (lf.int(x), lf.int(y));
            let mut acc = lf.int(rng.next_u64() as i64);
            let acc0 = acc;
            lf.fma(&mut acc, a, b);
            for (l, f) in lf.fields().iter().enumerate() {
                let (fa, fb) = (f.int(x), f.int(y));
                assert_eq!(a.0[l], fa);
                assert_eq!(lf.add(a, b).0[l], f.add(fa, fb));
                assert_eq!(lf.sub(a, b).0[l], f.sub(fa, fb));
                assert_eq!(lf.neg(a).0[l], f.neg(fa));
                assert_eq!(lf.mul(a, b).0[l], f.mul(fa, fb));
                assert_eq!(acc.0[l], f.add(acc0.0[l], f.mul(fa, fb)));
                assert_eq!(lf.residues(a)[l], f.decode(fa));
            }
        }
        assert_eq!(lf.residues(lf.one()), [1; 3]);
        assert!(lf.is_zero(Lanes::default()));
        assert!(!lf.is_zero(lf.one()));
        let x = lf.int(-12_345);
        assert_eq!(lf.mul(x, lf.inv(x).unwrap()), lf.one());
        // Zero in one lane only: no inverse, but not zero either.
        let seven = lf.int(7);
        assert_eq!(lf.inv(seven), None);
        assert!(!lf.is_zero(seven));
    }

    #[test]
    fn test_inverses() {
        let lf = LaneField::new([primes().next().unwrap(), 65_537]);
        for n in [0, 1, 2, 3, 100] {
            let inv = lf.inverses(n).unwrap();
            assert_eq!(inv.len(), n);
            for (i, &x) in inv.iter().enumerate() {
                let expected = if i == 0 {
                    Lanes::default()
                } else {
                    lf.inv(lf.int(i as i64)).unwrap()
                };
                assert_eq!(x, expected, "1/{i}");
            }
        }
        // 7 has no inverse modulo 7.
        assert_eq!(LaneField::new([7, 11]).inverses(8), None);
        assert!(LaneField::new([7, 11]).inverses(7).is_some());
    }

    #[test]
    #[should_panic(expected = "is repeated")]
    fn test_lanes_reject_repeated_prime() {
        LaneField::new([7, 11, 7]);
    }
}
