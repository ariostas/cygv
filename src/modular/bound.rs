//! Upward-rounded bounds, the number system of the certifier's majorant run (`PLAN.md`
//! phase 8).
//!
//! The majorant run computes the pipeline over nonnegative reals, with subtraction replaced
//! by addition and negation by the identity, so that every value bounds the absolute value of
//! the exact one. [`BoundRing`] does that arithmetic on [`Bound`]s: an `f32` mantissa with a
//! separate `i32` exponent, since the factorial tables leave the range of floating point
//! from $35!$ on, each operation rounded to nearest and then moved up by one ulp, which keeps it
//! an upper bound. That costs one ulp, $2^{-23}$, per operation: a sum of a million terms is
//! at most a third of a bit above its exact value, far below the slack of the bounds
//! themselves, while a bound takes 8 bytes, which halves the memory of the majorant pass's
//! tables. (An `f64` mantissa with an `i64` exponent took 16 bytes; malachite's `Float` rounds
//! up exactly, but is not `Copy` and is about nine times slower.)

use super::field::Ring;
use malachite::base::num::logic::traits::SignificantBits;
use malachite::{Integer, Rational};

/// A nonnegative number $m \, 2^e$, with $m \in [1, 2)$, or zero, with $m = e = 0$.
///
/// The exponent is never anywhere near overflowing: the values are bounds on numbers of at
/// most millions of bits.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Bound {
    m: f32,
    e: i32,
}

impl Bound {
    /// $m \, 2^e$ for a positive normal `m`, normalized, which is exact: it only moves the
    /// exponent of `m` into `e`.
    #[inline]
    fn new(m: f32, e: i32) -> Self {
        const EXP: u32 = 0xff << 23;
        let bits = m.to_bits();
        Self {
            m: f32::from_bits((bits & !EXP) | (127 << 23)),
            e: e + ((bits & EXP) >> 23) as i32 - 127,
        }
    }

    /// A bound on the absolute value of `n`: its leading 24 bits, moved up by one ulp if any
    /// bit below them is set.
    pub fn from_integer(n: &Integer) -> Self {
        let a = n.unsigned_abs_ref();
        let bits = a.significant_bits();
        if bits == 0 {
            return Self::default();
        }
        let shift = bits.saturating_sub(24);
        let top = u32::try_from(&(a >> shift)).expect("24 bits fit");
        let m = if a.trailing_zeros().is_some_and(|z| z < shift) {
            (top as f32).next_up()
        } else {
            top as f32
        };
        Self::new(m, i32::try_from(shift).expect("the exponent fits"))
    }

    /// The base-2 logarithm, approximately, or $-\infty$ for zero.
    pub fn log2(self) -> f64 {
        if self.m == 0.0 {
            f64::NEG_INFINITY
        } else {
            f64::from(self.m).log2() + f64::from(self.e)
        }
    }

    /// The exact value.
    pub fn to_rational(self) -> Rational {
        Rational::try_from(f64::from(self.m)).expect("the mantissa is finite") << self.e
    }
}

/// The arithmetic of [`Bound`]s, rounded up: subtraction is addition, negation the identity,
/// and every result is at least the exact result of the operation on the operands.
///
/// [`Ring::inv`] rounds up as well, so it is only a bound on the inverse of an exact value;
/// the inverse of a bound would round the wrong way.
#[derive(Debug, Clone, Copy, Default)]
pub struct BoundRing;

impl Ring for BoundRing {
    type Elem = Bound;

    const EXACT: bool = false;

    #[inline]
    fn one(&self) -> Bound {
        Bound { m: 1.0, e: 0 }
    }

    /// The absolute value, rounded up above $2^{53}$.
    #[inline]
    fn int(&self, v: i64) -> Bound {
        let a = v.unsigned_abs();
        match a {
            0 => Bound::default(),
            1..0x100_0000 => Bound::new(a as f32, 0),
            _ => Bound::new((a as f32).next_up(), 0),
        }
    }

    #[inline]
    fn add(&self, a: Bound, b: Bound) -> Bound {
        if a.m == 0.0 {
            return b;
        }
        if b.m == 0.0 {
            return a;
        }
        let (a, b) = if a.e >= b.e { (a, b) } else { (b, a) };
        let d = a.e - b.e;
        if d > 30 {
            // b < 2^(b.e + 1) <= 2^(a.e - 30), below the ulp 2^(a.e - 23) of a.
            Bound::new(a.m.next_up(), a.e)
        } else {
            // Exact: b.m times a power of two no smaller than 2^-30.
            let b = b.m * f32::from_bits((127 - d as u32) << 23);
            Bound::new((a.m + b).next_up(), a.e)
        }
    }

    #[inline]
    fn sub(&self, a: Bound, b: Bound) -> Bound {
        self.add(a, b)
    }

    #[inline]
    fn neg(&self, a: Bound) -> Bound {
        a
    }

    #[inline]
    fn mul(&self, a: Bound, b: Bound) -> Bound {
        if a.m == 0.0 || b.m == 0.0 {
            Bound::default()
        } else {
            Bound::new((a.m * b.m).next_up(), a.e + b.e)
        }
    }

    /// A bound on the inverse of the exact value `a`, or `None` for zero.
    #[inline]
    fn inv(&self, a: Bound) -> Option<Bound> {
        (a.m != 0.0).then(|| Bound::new((1.0 / a.m).next_up(), -a.e))
    }

    /// Each inverse taken by itself, since the default's inverse of $(n - 1)!$ would invert a
    /// bound.
    fn inverses(&self, n: usize) -> Option<Vec<Bound>> {
        Some(
            (0..n)
                .map(|i| match i {
                    0 => Bound::default(),
                    _ => self.inv(self.int(i as i64)).unwrap(),
                })
                .collect(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modular::period::Tables;
    use crate::modular::TestRng;

    /// Checks that `got` bounds `exact` from above, within `ulps` units of $2^{-23}$.
    fn check(got: Bound, exact: &Rational, ulps: u32, what: &str) {
        let g = got.to_rational();
        assert!(g >= *exact, "{what}: {got:?} is below the exact value");
        let slack = Rational::from(1) + Rational::from(ulps) / Rational::from(1u64 << 23);
        assert!(g <= exact * slack, "{what}: {got:?} is too far above it");
        if got.m != 0.0 {
            assert!(
                (1.0..2.0).contains(&got.m),
                "{what}: {got:?} is not normalized"
            );
        }
    }

    fn random(rng: &mut TestRng) -> Bound {
        if rng.range(0, 9) == 0 {
            return Bound::default();
        }
        let m = 1.0 + (rng.next_u64() >> 40) as f32 / (1u32 << 24) as f32 * 0.999;
        Bound::new(m, rng.range(-80, 80) as i32)
    }

    #[test]
    fn test_operations() {
        let r = BoundRing;
        let mut rng = TestRng::new(7);
        for _ in 0..10_000 {
            let (a, b) = (random(&mut rng), random(&mut rng));
            let (x, y) = (a.to_rational(), b.to_rational());
            check(r.add(a, b), &(&x + &y), 2, "add");
            assert_eq!(r.sub(a, b), r.add(a, b));
            assert_eq!(r.neg(a), a);
            check(r.mul(a, b), &(&x * &y), 2, "mul");
            let mut acc = a;
            r.fma(&mut acc, b, a);
            check(acc, &(&x + &y * &x), 4, "fma");
            match r.inv(a) {
                Some(i) => check(i, &(Rational::from(1) / &x), 2, "inv"),
                None => assert_eq!(x, 0u32),
            }
        }
        // Far apart: the smaller summand only moves the larger one up by an ulp.
        let (big, small) = (Bound::new(1.5, 100), Bound::new(1.0, 0));
        check(
            r.add(big, small),
            &(big.to_rational() + small.to_rational()),
            2,
            "far",
        );
        check(
            r.add(small, big),
            &r.add(big, small).to_rational(),
            0,
            "far, swapped",
        );
    }

    #[test]
    fn test_int() {
        let r = BoundRing;
        assert_eq!(r.int(0), Bound::default());
        assert_eq!(r.one(), r.int(1));
        for v in [
            1,
            -1,
            6,
            -1000,
            1 << 24,
            (1 << 24) + 1,
            1 << 53,
            (1 << 53) + 1,
            i64::MAX,
            i64::MIN,
        ] {
            let exact = Rational::from(v.unsigned_abs());
            check(r.int(v), &exact, 2, &format!("int {v}"));
        }
        assert_eq!(r.int(-12), r.int(12));
    }

    #[test]
    fn test_from_integer() {
        for n in [
            Integer::from(0),
            Integer::from(-5),
            Integer::from((1u64 << 24) - 1),
            Integer::from(1u64 << 24),
            Integer::from((1u64 << 24) + 1),
            Integer::from((1u64 << 53) - 1),
            Integer::from(1u64 << 53),
            Integer::from((1u64 << 53) + 1),
            Integer::from(u64::MAX),
            -(Integer::from(3) << 200u32),
            (Integer::from(7) << 3000u32) + Integer::from(1),
        ] {
            let exact = Rational::from(n.unsigned_abs_ref());
            check(Bound::from_integer(&n), &exact, 1, &format!("{n}"));
        }
        // Exact when the bits fit.
        let n = Integer::from(12345) << 100u32;
        assert_eq!(Bound::from_integer(&n).to_rational(), Rational::from(n));
    }

    #[test]
    fn test_tables() {
        // Far beyond the range of f64, which ends at 170!.
        let r = BoundRing;
        let n = 2000;
        let tables = Tables::new(&r, n + 1).unwrap();
        let mut fact = Rational::from(1);
        for i in 1..=n {
            fact *= Rational::from(i);
            if i % 250 == 0 || i < 10 {
                let ulps = 4 * i as u32;
                check(tables.fact[i], &fact, ulps, "fact");
                check(
                    tables.inv_fact[i],
                    &(Rational::from(1) / &fact),
                    ulps,
                    "inv_fact",
                );
                check(
                    tables.inv[i],
                    &(Rational::from(1) / Rational::from(i)),
                    2,
                    "inv",
                );
            }
        }
        assert!((tables.fact[n].log2() - 19052.987).abs() < 0.01);
        assert_eq!(tables.inv[0], Bound::default());
    }
}
