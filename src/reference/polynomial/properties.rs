//! A module for properties of polynomials.

use crate::reference::polynomial::coefficient::PolynomialCoeff;
use crate::reference::semigroup::Semigroup;
use nalgebra::DVectorView;
use rustc_hash::FxHashMap as HashMap;

/// A structure containing properties of a polynomial.
///
/// This structure contains the data common to all polynomials and is needed
/// for most operations.
#[derive(Clone, Debug)]
pub struct PolynomialProperties<'a, T> {
    pub semigroup: &'a Semigroup,
    pub monomial_map: HashMap<DVectorView<'a, i32>, usize>,
    /// A zero coefficient, kept around as the template new coefficients are
    /// cloned from.
    pub zero: T,
}

impl<'a, T: PolynomialCoeff<T>> PolynomialProperties<'a, T> {
    /// Create a new PolynomialProperties structure.
    pub fn new(semigroup: &'a Semigroup) -> Self {
        let mut poly_props = Self {
            semigroup,
            monomial_map: HashMap::default(),
            zero: T::default(),
        };

        for (i, c) in poly_props.semigroup.elements.column_iter().enumerate() {
            poly_props.monomial_map.insert(c, i);
        }

        poly_props
    }
}
