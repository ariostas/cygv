//! A module for errors that can happen in the computation of the fundamental period.

use core::fmt;

/// An error structure for non-positive degrees of non-identity elements.
#[derive(Debug, Clone, PartialEq)]
pub enum FundamentalPeriodError {
    CYDimLessThanThree,
    InconsistentNefPartition,
    /// A curve has negative anticanonical degree but fewer than two negative intersections
    /// with the toric divisors.
    NegativeAnticanonicalDegree,
    /// The strata do not fit the geometry, for the given reason.
    InvalidStrata(&'static str),
}

/// Implement Display trait for SemigroupError.
impl fmt::Display for FundamentalPeriodError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let message = match self {
            FundamentalPeriodError::CYDimLessThanThree => {
                "the dimension of the CY must be at least three"
            }
            FundamentalPeriodError::InconsistentNefPartition => {
                "the nef partition data is inconsistent"
            }
            FundamentalPeriodError::NegativeAnticanonicalDegree => {
                "a curve has negative anticanonical degree but fewer than two negative \
                 intersections with the toric divisors"
            }
            FundamentalPeriodError::InvalidStrata(reason) => reason,
        };
        write!(f, "{message}")
    }
}
