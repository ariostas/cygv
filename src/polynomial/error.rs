//! A module for errors that can happen in polynomial operations.

use core::fmt;

/// An error enum for polynomial errors.
#[derive(Debug, Clone)]
pub enum PolynomialError {
    ZeroConstantTermError,
    NonZeroConstantTermError,
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
        };
        write!(f, "{message}")
    }
}
