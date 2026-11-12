//! The exact rational numbers the reference pipeline computes with.

use crate::reference::factorial::RecipMut;
use crate::reference::polynomial::coefficient::{Assign, IsInteger, PolynomialCoeff};
use core::cmp::Ordering;
use core::fmt;
use core::ops::{AddAssign, DivAssign, MulAssign, SubAssign};
use malachite::base::num::arithmetic::traits::ReciprocalAssign;
use malachite::base::num::conversion::traits::IsInteger as _;
use malachite::{Integer, Rational as Inner};

/// An exact rational number, the coefficient type of the reference pipeline.
///
/// A thin wrapper of [`malachite::Rational`] that adds the in-place operations with primitive
/// integers [`PolynomialCoeff`] asks for. [`Default`] gives zero.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Rational(pub Inner);

impl Rational {
    /// Zero.
    pub fn new() -> Self {
        Self::default()
    }

    /// The wrapped number.
    pub fn into_inner(self) -> Inner {
        self.0
    }

    /// The numerator, which is the number itself when it is an integer.
    pub fn numerator(&self) -> Integer {
        let n = Integer::from(self.0.numerator_ref());
        if self.0 < 0u32 {
            -n
        } else {
            n
        }
    }
}

impl From<Inner> for Rational {
    fn from(x: Inner) -> Self {
        Self(x)
    }
}

impl From<Rational> for Inner {
    fn from(x: Rational) -> Self {
        x.0
    }
}

impl fmt::Display for Rational {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.0.fmt(f)
    }
}

impl Assign<&Rational> for Rational {
    #[inline]
    fn assign(&mut self, x: &Rational) {
        self.0.clone_from(&x.0);
    }
}

impl AddAssign<&Rational> for Rational {
    #[inline]
    fn add_assign(&mut self, x: &Rational) {
        self.0 += &x.0;
    }
}

impl SubAssign<&Rational> for Rational {
    #[inline]
    fn sub_assign(&mut self, x: &Rational) {
        self.0 -= &x.0;
    }
}

impl MulAssign<&Rational> for Rational {
    #[inline]
    fn mul_assign(&mut self, x: &Rational) {
        self.0 *= &x.0;
    }
}

impl DivAssign<&Rational> for Rational {
    #[inline]
    fn div_assign(&mut self, x: &Rational) {
        self.0 /= &x.0;
    }
}

// Malachite stores integers that fit a limb inline, so converting a primitive does not allocate.
macro_rules! primitive {
    ($($t:ty),*) => {$(
        impl From<$t> for Rational {
            #[inline]
            fn from(x: $t) -> Self {
                Self(Inner::from(x))
            }
        }

        impl From<($t, $t)> for Rational {
            /// The fraction $n / d$, which must have $d \neq 0$.
            #[inline]
            fn from((n, d): ($t, $t)) -> Self {
                Self(Inner::from(n) / Inner::from(d))
            }
        }

        impl Assign<$t> for Rational {
            #[inline]
            fn assign(&mut self, x: $t) {
                self.0 = Inner::from(x);
            }
        }

        impl Assign<($t, $t)> for Rational {
            #[inline]
            fn assign(&mut self, x: ($t, $t)) {
                *self = Self::from(x);
            }
        }

        impl AddAssign<$t> for Rational {
            #[inline]
            fn add_assign(&mut self, x: $t) {
                self.0 += Inner::from(x);
            }
        }

        impl SubAssign<$t> for Rational {
            #[inline]
            fn sub_assign(&mut self, x: $t) {
                self.0 -= Inner::from(x);
            }
        }

        impl MulAssign<$t> for Rational {
            #[inline]
            fn mul_assign(&mut self, x: $t) {
                self.0 *= Inner::from(x);
            }
        }

        impl DivAssign<$t> for Rational {
            #[inline]
            fn div_assign(&mut self, x: $t) {
                self.0 /= Inner::from(x);
            }
        }

        impl PartialEq<$t> for Rational {
            #[inline]
            fn eq(&self, x: &$t) -> bool {
                self.0 == *x
            }
        }

        impl PartialOrd<$t> for Rational {
            #[inline]
            fn partial_cmp(&self, x: &$t) -> Option<Ordering> {
                self.0.partial_cmp(x)
            }
        }
    )*};
}

primitive!(i32, i64, u32, u64);

impl RecipMut for Rational {
    #[inline]
    fn recip_mut(&mut self) {
        self.0.reciprocal_assign();
    }
}

impl IsInteger for Rational {
    #[inline]
    fn is_integer(&self) -> bool {
        self.0.is_integer()
    }
}

impl PolynomialCoeff<Rational> for Rational {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_operations() {
        let mut x = Rational::from((3, 4));
        x += 1;
        assert_eq!(x, Rational::from((7, 4)));
        x *= 4u64;
        assert_eq!(x, 7);
        assert!(x.is_integer());
        x /= -14i64;
        assert_eq!(x, Rational::from((-1, 2)));
        assert!(x < 0 && x > -1);
        x.recip_mut();
        assert_eq!(x, -2);
        x -= &Rational::from((1, 3));
        assert_eq!(x.to_string(), "-7/3");
        assert_eq!(x.numerator(), -7);
        x.assign(5u32);
        assert_eq!(x.numerator(), 5);
        assert_eq!(x.to_string(), "5");
    }
}
