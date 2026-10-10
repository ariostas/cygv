//! Polynomial coefficients

use crate::reference::factorial::RecipMut;
use core::ops::{AddAssign, DivAssign, MulAssign, SubAssign};
use rug::{Assign, Rational};
use std::fmt::Display;

/// Tells whether a number is an integer.
pub trait IsInteger {
    fn is_integer(&self) -> bool;
}

impl IsInteger for Rational {
    #[inline]
    fn is_integer(&self) -> bool {
        Rational::is_integer(self)
    }
}

/// A compound trait of all the traits a reasonable polynomial coefficient should have.
///
/// Coefficients are exact, and [`Default`] must give zero.
pub trait PolynomialCoeff<T>:
    for<'a> Assign<&'a T>
    + for<'a> AddAssign<&'a T>
    + for<'a> SubAssign<&'a T>
    + for<'a> MulAssign<&'a T>
    + for<'a> DivAssign<&'a T>
    + Assign<i32>
    + Assign<i64>
    + Assign<u32>
    + Assign<u64>
    + AddAssign<i32>
    + AddAssign<i64>
    + AddAssign<u32>
    + AddAssign<u64>
    + SubAssign<i32>
    + SubAssign<i64>
    + SubAssign<u32>
    + SubAssign<u64>
    + MulAssign<i32>
    + MulAssign<i64>
    + MulAssign<u32>
    + MulAssign<u64>
    + DivAssign<i32>
    + DivAssign<i64>
    + DivAssign<u32>
    + DivAssign<u64>
    + PartialEq<i32>
    + PartialOrd<T>
    + PartialOrd<i32>
    + PartialOrd<i64>
    + PartialOrd<u32>
    + PartialOrd<u64>
    + Display
    + Clone
    + Default
    + Send
    + Sync
    + RecipMut
    + IsInteger
{
}

impl PolynomialCoeff<Rational> for Rational {}
