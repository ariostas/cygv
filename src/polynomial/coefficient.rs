//! Polynomial coefficients

use crate::factorial::RecipMut;
use core::ops::{AddAssign, DivAssign, MulAssign, SubAssign};
use rug::{Assign, Float, Rational};
use std::fmt::Display;

/// Rounds a number to the closest integer in place.
pub trait RoundMut {
    fn round_mut(&mut self);
}

impl RoundMut for Rational {
    #[inline]
    fn round_mut(&mut self) {
        Rational::round_mut(self);
    }
}

impl RoundMut for Float {
    #[inline]
    fn round_mut(&mut self) {
        Float::round_mut(self);
    }
}

/// Takes the absolute value in place.
pub trait AbsMut {
    fn abs_mut(&mut self);
}

impl AbsMut for Rational {
    #[inline]
    fn abs_mut(&mut self) {
        Rational::abs_mut(self);
    }
}

impl AbsMut for Float {
    #[inline]
    fn abs_mut(&mut self) {
        Float::abs_mut(self);
    }
}

/// A compound trait of all the traits a reasonable polynomial coefficient should reasonably have.
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
    + PartialOrd<f32>
    + PartialOrd<f64>
    + Display
    + Clone
    + Send
    + Sync
    + RecipMut
    + RoundMut
    + AbsMut
{
}

impl PolynomialCoeff<Rational> for Rational {}

impl PolynomialCoeff<Float> for Float {}
