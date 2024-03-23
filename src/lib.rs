pub mod factorial;
pub mod polynomial;
pub mod pool;
pub mod semigroup;

// Re-export the main structs
pub use polynomial::{
    prettyprint::PrettyPrintPolynomial, properties::PolynomialProperties, Polynomial,
};
pub use pool::NumberPool;
pub use semigroup::Semigroup;
