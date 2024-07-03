//! A module for properties of polynomials.

use std::vec;

use crate::semigroup::Semigroup;
use nalgebra::DVectorView;
use boomphf::hashmap::BoomHashMap;

/// A structure containing properties of a polynomial.
///
/// This structure contains the data common to all polynomaials and is needed
/// for most operations.
#[derive(Clone, Debug)]
pub struct PolynomialProperties<'a, T> {
    pub semigroup: &'a Semigroup,
    pub monomial_map: BoomHashMap<DVectorView<'a, i32>, u32>,
    pub zero_cutoff: T,
}

impl<'a, T: Clone> PolynomialProperties<'a, T> {
    /// Create a new PolynomialProperties structure.
    pub fn new(semigroup: &'a Semigroup, zero_cutoff: &T) -> Self {
        let mut poly_props = PolynomialProperties {
            semigroup,
            monomial_map: BoomHashMap::new(vec![], vec![]),
            zero_cutoff: zero_cutoff.clone(),
        };
        
        let monomial_views: Vec<_> = poly_props.semigroup.elements.column_iter().collect();
        let indices: Vec<u32> = (0..poly_props.semigroup.elements.ncols() as u32).collect();
        println!("Generating phm");
        poly_props.monomial_map = BoomHashMap::new(monomial_views, indices);
        println!("Done with phm");

        poly_props
    }
}
