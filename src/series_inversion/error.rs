//! A module for errors that can happen in series inversions.

use crate::polynomial::error::PolynomialError;
use core::fmt;

/// An error structure for semigroup errors.
#[derive(Debug, Clone, PartialEq)]
pub enum SeriesInversionError {
    NonIntegerGVError,
    PolynomialError(PolynomialError),
}

/// Implement Display trait for SeriesInversionError.
impl fmt::Display for SeriesInversionError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            SeriesInversionError::NonIntegerGVError => {
                write!(f, "A non-integer GV invariant was found")
            }
            SeriesInversionError::PolynomialError(e) => e.fmt(f),
        }
    }
}

impl From<PolynomialError> for SeriesInversionError {
    fn from(e: PolynomialError) -> Self {
        SeriesInversionError::PolynomialError(e)
    }
}
