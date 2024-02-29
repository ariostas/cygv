//! A module for errors that can happen in semigroup constructions.

use core::fmt;

/// An error structure for non-positive degrees of non-identity elements.
#[derive(Debug, Clone, PartialEq)]
pub enum SemigroupError {
    NonPositiveDegreeError,
    MissingIdentityError,
}

/// Implement Display trait for SemigroupError.
impl fmt::Display for SemigroupError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let message = match self {
            SemigroupError::NonPositiveDegreeError => {
                "the degree of all non-zero elements must be positive"
            }
            SemigroupError::MissingIdentityError => "the identity element is missing",
        };
        write!(f, "{message}")
    }
}
