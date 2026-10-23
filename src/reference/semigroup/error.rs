//! A module for errors that can happen in semigroup constructions.

use core::fmt;

/// An error structure for semigroup errors.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SemigroupError {
    NonPositiveDegreeError,
    MissingIdentityError,
    DimensionMismatchError,
    DuplicateElementError,
    DegreeTooLargeError,
}

/// Implement Display trait for SemigroupError.
impl fmt::Display for SemigroupError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let message = match self {
            SemigroupError::NonPositiveDegreeError => {
                "the degree of all non-zero elements must be positive"
            }
            SemigroupError::MissingIdentityError => "the identity element is missing",
            SemigroupError::DimensionMismatchError => {
                "the dimensions of the semigroup data do not match"
            }
            SemigroupError::DuplicateElementError => "an element is given more than once",
            SemigroupError::DegreeTooLargeError => "a degree does not fit in 32 bits",
        };
        write!(f, "{message}")
    }
}
