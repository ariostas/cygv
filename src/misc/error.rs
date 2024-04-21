//! A module for errors that can happen in miscellaneous computations.

use core::fmt;

/// An error enum for miscellaneous errors.
#[derive(Debug, Clone)]
pub enum MiscError {
    EmptyIntNums,
    WrongDimIntNums,
    NegativeIndex,
    RepeatedIdxIntNums,
}

/// Implement Display trait for PolynomialError.
impl fmt::Display for MiscError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let message = match self {
            MiscError::EmptyIntNums => {
                "there must be at least one nonzero intersection number"
            }
            MiscError::WrongDimIntNums => {
                "the intersection number must be given as a triplet of indices followed by the intersection number"
            }
            MiscError::NegativeIndex => {
                "indices of intersection numbers must be positive"
            }
            MiscError::RepeatedIdxIntNums => {
                "the input intersection numbers must be uniquely specified"
            }
        };
        write!(f, "{message}")
    }
}
