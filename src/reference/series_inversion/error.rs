//! A module for errors that can happen in series inversions.

use crate::reference::polynomial::error::PolynomialError;
use core::fmt;

/// An error structure for semigroup errors.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SeriesInversionError {
    NonIntegerGVError,
    PolynomialError(PolynomialError),
}

/// Implement Display trait for SeriesInversionError.
impl fmt::Display for SeriesInversionError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            SeriesInversionError::NonIntegerGVError => {
                write!(
                    f,
                    "a GV invariant is not an integer; {}",
                    crate::NON_INTEGER_CAUSES
                )
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
