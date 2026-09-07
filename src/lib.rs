#![doc = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/README.md"))]

pub mod factorial;
pub mod fundamental_period;
pub mod hkty;
pub mod instanton;
#[cfg(feature = "cli")]
pub mod io;
pub mod misc;
pub mod polynomial;
pub mod pool;
#[cfg(feature = "python")]
pub mod python;
pub mod semigroup;
pub mod series_inversion;

/// Which invariants to compute.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum InvariantKind {
    /// Gopakumar-Vafa invariants, which are integers.
    GV,
    /// Gromov-Witten invariants, which are rational.
    GW,
}

impl InvariantKind {
    /// The name of these invariants in the input and output documents.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::GV => "gv",
            Self::GW => "gw",
        }
    }
}

/// The kind of Calabi-Yau that the invariants are computed for.
///
/// Threefolds are treated separately throughout. Their intersection numbers are
/// canonicalized as unordered triplets, and their invariants are labeled by a curve
/// class alone. In higher dimensions the first index of the intersection numbers is
/// a reference surface, which labels the resulting invariants alongside the curve
/// class.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CYKind {
    /// A CY threefold.
    Threefold,
    /// A CY of dimension greater than three.
    Nfold,
}

impl CYKind {
    /// Selects the variant from the `is_threefold` flag of the input formats.
    pub fn from_is_threefold(is_threefold: bool) -> Self {
        if is_threefold {
            Self::Threefold
        } else {
            Self::Nfold
        }
    }

    /// Whether the CY is a threefold.
    pub fn is_threefold(self) -> bool {
        matches!(self, Self::Threefold)
    }
}

// Re-export the main structs
#[doc(inline)]
pub use polynomial::{
    prettyprint::PrettyPrintPolynomial, properties::PolynomialProperties, Polynomial,
};
#[doc(inline)]
pub use pool::NumberPool;
#[doc(inline)]
pub use semigroup::Semigroup;

// Re-export main trait
#[doc(inline)]
pub use polynomial::coefficient::PolynomialCoeff;

// Re-export python module
#[cfg(feature = "python")]
#[doc(inline)]
pub use python::cygv;

// Re-export main functions
#[doc(inline)]
pub use hkty::*;
