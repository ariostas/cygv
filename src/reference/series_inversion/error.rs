//! A module for errors that can happen in series inversions.

use crate::reference::polynomial::error::PolynomialError;
use core::fmt;

/// An error structure for semigroup errors.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SeriesInversionError {
    NonIntegerGVError,
    /// The instanton corrections of a threefold do not give the same invariant for a curve
    /// class, which only an inconsistent geometry does.
    Inconsistent,
    PolynomialError(PolynomialError),
    /// The computation was interrupted ([`crate::Interrupt`]).
    Interrupted,
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
            SeriesInversionError::Inconsistent => {
                write!(
                    f,
                    "the instanton corrections do not agree on an invariant; {}",
                    crate::NON_INTEGER_CAUSES
                )
            }
            SeriesInversionError::PolynomialError(e) => e.fmt(f),
            SeriesInversionError::Interrupted => write!(f, "the computation was interrupted"),
        }
    }
}

impl From<PolynomialError> for SeriesInversionError {
    fn from(e: PolynomialError) -> Self {
        match e {
            PolynomialError::Interrupted => SeriesInversionError::Interrupted,
            e => SeriesInversionError::PolynomialError(e),
        }
    }
}
