//! A module for properties of polynomials.

use crate::semigroup::Semigroup;
use nalgebra::DVectorView;
use std::collections::HashMap;

/// A structure containing properties of a polynomial.
///
/// This structure contains the data common to all polynomaials and is needed
/// for most operations.
#[derive(Clone, Debug)]
pub struct PolynomialProperties<'a, T> {
    pub semigroup: &'a Semigroup,
    pub monomial_map: HashMap<DVectorView<'a, i32>, usize>,
    pub zero_cutoff: T,
}

impl<'a, T: Clone> PolynomialProperties<'a, T> {
    /// Create a new PolynomialProperties structure.
    pub fn new(semigroup: &'a Semigroup, zero_cutoff: &T) -> Self {
        let mut poly_props = Self {
            semigroup,
            monomial_map: HashMap::new(),
            zero_cutoff: zero_cutoff.clone(),
        };

        for (i, c) in poly_props.semigroup.elements.column_iter().enumerate() {
            poly_props.monomial_map.insert(c, i);
        }

        poly_props
    }
}
