//! A module for properties of polynomials.

use crate::polynomial::coefficient::PolynomialCoeff;
use crate::semigroup::Semigroup;
use nalgebra::DVectorView;
use std::collections::HashMap;

/// A structure containing properties of a polynomial.
///
/// This structure contains the data common to all polynomials and is needed
/// for most operations.
#[derive(Clone, Debug)]
pub struct PolynomialProperties<'a, T> {
    pub semigroup: &'a Semigroup,
    pub monomial_map: HashMap<DVectorView<'a, i32>, usize>,
    pub zero_cutoff: T,
    /// A zero coefficient, kept around as the template new coefficients are
    /// cloned from. For [`rug::Float`] coefficients it carries the precision
    /// they all have to be created with.
    pub zero: T,
}

impl<'a, T: PolynomialCoeff<T>> PolynomialProperties<'a, T> {
    /// Create a new PolynomialProperties structure.
    pub fn new(semigroup: &'a Semigroup, zero_cutoff: &T) -> Self {
        let mut zero = zero_cutoff.clone();
        zero.assign(0i32);
        let mut poly_props = Self {
            semigroup,
            monomial_map: HashMap::new(),
            zero_cutoff: zero_cutoff.clone(),
            zero,
        };

        for (i, c) in poly_props.semigroup.elements.column_iter().enumerate() {
            poly_props.monomial_map.insert(c, i);
        }

        poly_props
    }
}
