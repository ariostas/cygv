pub mod factorial;
pub mod fundamental_period;
pub mod instanton;
pub mod misc;
pub mod polynomial;
pub mod pool;
pub mod semigroup;
pub mod series_inversion;

// Re-export the main structs
pub use polynomial::{
    prettyprint::PrettyPrintPolynomial, properties::PolynomialProperties, Polynomial,
};
pub use pool::NumberPool;
pub use semigroup::Semigroup;

// Re-export main trait
pub use polynomial::coefficient::PolynomialCoeff;
