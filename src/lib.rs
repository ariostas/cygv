#![doc = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/README.md"))]

pub mod factorial;
pub mod fundamental_period;
pub mod hkty;
pub mod instanton;
pub mod misc;
pub mod polynomial;
pub mod python;
pub mod semigroup;
pub mod series_inversion;

// Re-export the main structs
#[doc(inline)]
pub use polynomial::{
    prettyprint::PrettyPrintPolynomial, properties::PolynomialProperties, Polynomial,
};
#[doc(inline)]
pub use semigroup::Semigroup;

// Re-export main trait
#[doc(inline)]
pub use polynomial::coefficient::PolynomialCoeff;

// Re-export python module
#[doc(inline)]
pub use python::cygv;

// Re-export main functions
#[doc(inline)]
pub use hkty::*;
