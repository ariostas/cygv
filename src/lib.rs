/*!
# GV and GW invariants of CY manifolds

[![Rust CI](https://github.com/ariostas/cygv/actions/workflows/rust.yml/badge.svg)](https://github.com/ariostas/cygv/actions/workflows/rust.yml)

**Warning**: This project is still in early stages. The code and documentation are under active development and may change significantly.

This project implements an efficient algorithm to perform the HKTY procedure [[1], [2]] to compute Gopakumar-Vafa (GV) and Gromov-Witten (GW) invariants of Calabi-Yau (CY) manifolds obtained as hypersurfaces or complete intersections in toric varieties. This project is based on the work presented in the paper "[Computational Mirror Symmetry]", but written in the Rust programming language and with some additional improvements.

## License

Licensed under either of

 * Apache License, Version 2.0
   (<http://www.apache.org/licenses/LICENSE-2.0>)
 * MIT license
   (<http://opensource.org/licenses/MIT>)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.

[1]: https://arxiv.org/abs/hep-th/9308122
[2]: https://arxiv.org/abs/hep-th/9406055
[Computational Mirror Symmetry]: https://arxiv.org/abs/2303.00757
*/

pub mod factorial;
pub mod fundamental_period;
pub mod hkty;
pub mod instanton;
pub mod misc;
pub mod polynomial;
pub mod pool;
pub mod python;
pub mod semigroup;
pub mod series_inversion;

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
#[doc(inline)]
pub use python::cygv;

// Re-export main functions
#[doc(inline)]
pub use hkty::*;
