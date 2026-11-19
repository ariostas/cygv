//! A module for errors that can happen in polynomial operations.

use core::fmt;

/// An error enum for polynomial errors.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PolynomialError {
    ZeroConstantTermError,
    NonZeroConstantTermError,
    /// The computation was interrupted ([`crate::Interrupt`]).
    Interrupted,
}

/// Implement Display trait for PolynomialError.
impl fmt::Display for PolynomialError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let message = match self {
            PolynomialError::ZeroConstantTermError => {
                "this operation requires a non-zero constant term"
            }
            PolynomialError::NonZeroConstantTermError => {
                "this operation requires a zero constant term"
            }
            PolynomialError::Interrupted => "the computation was interrupted",
        };
        write!(f, "{message}")
    }
}
