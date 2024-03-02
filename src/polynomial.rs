//! Polynomial
//!
//! This module provides tools to work with polynomials containing a fixed set
//! of monomials.

pub mod error;
pub mod properties;

use crate::cache::NumCache;
use core::ops::{AddAssign, DivAssign, MulAssign, SubAssign};
use error::PolynomialError;
use nalgebra::DVector;
use properties::PolynomialProperties;
use rug::Assign;
use std::collections::HashMap;

/// A polynomial structure.
///
/// The structure contains a hash map from monomial indices to the value of
/// the coefficient. It also contains a vector listing the sorted indices of
/// its nonzero coefficients. Even though this is redundant information, it is
/// useful to have this list to make some operations more efficient.
#[derive(Clone, Debug, PartialEq)]
pub struct Polynomial<T> {
    pub coeffs: HashMap<u32, T>,
    pub nonzero: Vec<u32>,
}

impl<T> Polynomial<T>
where
    for<'a> T: Clone
        + AddAssign<&'a T>
        + Assign<i32>
        + Assign<&'a T>
        + DivAssign<&'a T>
        + DivAssign<u32>
        + MulAssign<&'a T>
        + MulAssign<i32>
        + MulAssign<u32>
        + PartialEq<i32>
        + PartialOrd<T>
        + SubAssign<&'a T>,
{
    /// Create a new polynomial.
    pub fn new() -> Self {
        Self {
            coeffs: HashMap::new(),
            nonzero: Vec::new(),
        }
    }

    /// Create a new polynomial with one as its constant term.
    pub fn one(coeff_cache: &mut NumCache<T>) -> Self {
        let mut coeffs = HashMap::new();
        let mut one = coeff_cache.pop();
        one.assign(1i32);
        coeffs.insert(0, one);
        Self {
            coeffs,
            nonzero: vec![0],
        }
    }

    /// Drop the polynomial, while moving the allocated coefficients to the
    /// cache.
    pub fn drop(self, coeff_cache: &mut NumCache<T>) {
        for c in self.coeffs.into_values() {
            coeff_cache.push(c);
        }
    }

    /// Clean up the polynomial, removing monomials whose coefficient are
    /// below a given cutoff.
    pub fn clean_up(
        &mut self,
        poly_props: &PolynomialProperties<T>,
        coeff_cache: &mut NumCache<T>,
    ) {
        let mut new_nonzero = Vec::new();
        let mut to_delete = Vec::new();
        let mut pos_cutoff = coeff_cache.pop();
        let mut neg_cutoff = coeff_cache.pop();
        pos_cutoff.assign(&poly_props.zero_cutoff);
        neg_cutoff.assign(&poly_props.zero_cutoff);
        neg_cutoff *= -1;
        for (i, c) in self.coeffs.iter() {
            if *c <= pos_cutoff && *c >= neg_cutoff {
                to_delete.push(*i);
            } else {
                new_nonzero.push(*i);
            }
        }
        for i in to_delete.iter() {
            self.coeffs.remove(i);
        }
        self.nonzero = new_nonzero;
        coeff_cache.push(pos_cutoff);
        coeff_cache.push(neg_cutoff);
    }

    /// Find the minimum degree of the polynomial.
    pub fn min_degree(&self, poly_props: &PolynomialProperties<T>) -> u32 {
        for &i in self.nonzero.iter() {
            return poly_props.semigroup.degrees[i as usize];
        }
        poly_props.semigroup.max_degree + 1
    }

    /// Clone the polynomial, but truncated to a given degree.
    pub fn truncated(
        &self,
        max_deg: u32,
        poly_props: &PolynomialProperties<T>,
        coeff_cache: &mut NumCache<T>,
    ) -> Self {
        let mut res = Self::new();
        for &i in self.nonzero.iter() {
            if poly_props.semigroup.degrees[i as usize] > max_deg {
                break;
            }
            res.nonzero.push(i);
            let mut coeff = coeff_cache.pop();
            coeff.assign(self.coeffs.get(&i).unwrap());
            res.coeffs.insert(i, coeff);
        }
        res
    }

    /// Add polynomials in place.
    pub fn add_assign(&mut self, rhs: &Self, coeff_cache: &mut NumCache<T>) {
        let mut resort = false;
        for (k, v) in rhs.coeffs.iter() {
            let c = self.coeffs.entry(*k).or_insert_with(|| {
                resort = true;
                let mut new_coeff = coeff_cache.pop();
                new_coeff.assign(0i32);
                new_coeff
            });
            *c += v;
        }
        if resort {
            self.nonzero = self.coeffs.keys().cloned().collect();
            self.nonzero.sort_unstable();
        }
    }

    /// Subtract polynomials in place.
    pub fn sub_assign(&mut self, rhs: &Self, coeff_cache: &mut NumCache<T>) {
        let mut resort = false;
        for (k, v) in rhs.coeffs.iter() {
            let c = self.coeffs.entry(*k).or_insert_with(|| {
                resort = true;
                let mut new_coeff = coeff_cache.pop();
                new_coeff.assign(0i32);
                new_coeff
            });
            *c -= v;
        }
        if resort {
            self.nonzero = self.coeffs.keys().cloned().collect();
            self.nonzero.sort_unstable();
        }
    }

    /// Multiply polynomial by a scalar in place.
    pub fn mul_scalar_assign<U>(&mut self, scalar: U)
    where
        T: MulAssign<U>,
        U: Copy,
    {
        for c in self.coeffs.values_mut() {
            *c *= scalar;
        }
    }

    /// Divide polynomial by a scalar in place.
    pub fn div_scalar_assign<U>(&mut self, scalar: U)
    where
        T: DivAssign<U>,
        U: Copy,
    {
        for c in self.coeffs.values_mut() {
            *c /= scalar;
        }
    }

    /// Multiply polynomials.
    pub fn mul(
        &self,
        rhs: &Self,
        poly_props: &PolynomialProperties<T>,
        coeff_cache: &mut NumCache<T>,
    ) -> Self {
        let mut res = Self::new();
        let mut deg1;
        let mut deg2;
        let max_deg = poly_props.semigroup.max_degree;
        let mut tmp_vec = DVector::zeros(poly_props.semigroup.elements.nrows() as usize);
        let (pshort, plong) = if self.nonzero.len() < rhs.nonzero.len() {
            (self, rhs)
        } else {
            (rhs, self)
        };
        let mut tmp_var = coeff_cache.pop();
        for &i in pshort.nonzero.iter() {
            deg1 = poly_props.semigroup.degrees[i as usize];
            for &j in plong.nonzero.iter() {
                deg2 = poly_props.semigroup.degrees[j as usize];
                if deg1 + deg2 > max_deg {
                    continue;
                }
                poly_props
                    .semigroup
                    .elements
                    .column(i as usize)
                    .iter()
                    .zip(poly_props.semigroup.elements.column(j as usize).iter())
                    .zip(tmp_vec.iter_mut())
                    .for_each(|((x, y), z)| *z = *x + *y);
                let Some(mon) = poly_props.monomial_map.get(&tmp_vec.as_view()) else {
                    continue;
                };
                let c = res.coeffs.entry(*mon).or_insert_with(|| {
                    res.nonzero.push(*mon);
                    let mut coeff = coeff_cache.pop();
                    coeff.assign(0i32);
                    coeff
                });
                tmp_var.assign(pshort.coeffs.get(&i).unwrap());
                tmp_var *= plong.coeffs.get(&j).unwrap();
                *c += &tmp_var;
            }
        }
        coeff_cache.push(tmp_var);
        res
    }

    /// Multiply polynomial by a scalar, and then add in place.
    pub fn mul_scalar_add_assign<U>(&mut self, scalar: U, rhs: &Self, coeff_cache: &mut NumCache<T>)
    where
        T: MulAssign<U>,
        U: Copy,
    {
        let mut tmp_var = coeff_cache.pop();
        let mut resort = false;
        for (k, v) in rhs.coeffs.iter() {
            let Some(c) = self.coeffs.get_mut(k) else {
                let mut new_var = coeff_cache.pop();
                new_var.assign(v);
                new_var *= scalar;
                self.nonzero.push(*k);
                self.coeffs.insert(*k, new_var);
                resort = true;
                continue;
            };
            tmp_var.assign(v);
            tmp_var *= scalar;
            *c += &tmp_var;
        }
        if resort {
            self.nonzero = self.coeffs.keys().cloned().collect();
            self.nonzero.sort_unstable();
        }
        coeff_cache.push(tmp_var);
    }

    /// Divide polynomial by a scalar, and then add in place.
    pub fn div_scalar_add_assign<U>(&mut self, scalar: U, rhs: &Self, coeff_cache: &mut NumCache<T>)
    where
        T: DivAssign<U>,
        U: Copy,
    {
        let mut tmp_var = coeff_cache.pop();
        let mut resort = false;
        for (k, v) in rhs.coeffs.iter() {
            let Some(c) = self.coeffs.get_mut(k) else {
                let mut new_var = coeff_cache.pop();
                new_var.assign(v);
                new_var /= scalar;
                self.nonzero.push(*k);
                self.coeffs.insert(*k, new_var);
                resort = true;
                continue;
            };
            tmp_var.assign(v);
            tmp_var /= scalar;
            *c += &tmp_var;
        }
        if resort {
            self.nonzero = self.coeffs.keys().cloned().collect();
            self.nonzero.sort_unstable();
        }
        coeff_cache.push(tmp_var);
    }

    /// Clone the polynomial.
    pub fn clone(&self, coeff_cache: &mut NumCache<T>) -> Self {
        let mut res = Self::new();
        for (i, c) in self.coeffs.iter() {
            res.nonzero.push(*i);
            let mut coeff = coeff_cache.pop();
            coeff.assign(c);
            res.coeffs.insert(*i, coeff);
        }
        res
    }

    /// Move all data into another polynomial.
    pub fn move_into(self, other: &mut Self, coeff_cache: &mut NumCache<T>) {
        for (_, c) in other.coeffs.drain() {
            coeff_cache.push(c);
        }
        other.coeffs = self.coeffs;
        other.nonzero = self.nonzero;
    }

    /// Compute the reciprocal of the polynomial.
    ///
    /// Only works for polynomials with a nonzero constant term.
    pub fn recipr(
        &self,
        poly_props: &PolynomialProperties<T>,
        coeff_cache: &mut NumCache<T>,
    ) -> Result<Self, PolynomialError> {
        let Some(c) = self.coeffs.get(&0) else {
            return Err(PolynomialError::ZeroConstantTermError);
        };
        if *c == 0 {
            return Err(PolynomialError::ZeroConstantTermError);
        }
        let mut res = Self::one(coeff_cache);
        let max_deg = poly_props.semigroup.max_degree;
        let mut p0 = self.clone(coeff_cache);
        p0.nonzero.remove(0);
        let const_term = p0.coeffs.remove(&0).unwrap();
        p0.div_scalar_assign(&const_term);
        let mut tmp_poly = Self::one(coeff_cache);
        let min_deg = p0.min_degree(poly_props);
        for i in 1..=max_deg / min_deg {
            let tmp_poly2 = tmp_poly.mul(&tmp_poly, poly_props, coeff_cache);
            tmp_poly2.move_into(&mut tmp_poly, coeff_cache);
            if i % 2 == 0 {
                res.add_assign(&tmp_poly, coeff_cache)
            } else {
                res.mul_scalar_add_assign(-1, &tmp_poly, coeff_cache);
            }
        }
        res.div_scalar_assign(&const_term);
        tmp_poly.drop(coeff_cache);
        Ok(res)
    }

    /// Compute powers of polynomials.
    pub fn pow(
        &self,
        n: i32,
        poly_props: &PolynomialProperties<T>,
        coeff_cache: &mut NumCache<T>,
    ) -> Result<Self, PolynomialError> {
        if n == 0 {
            return Ok(Self::one(coeff_cache));
        } else if n == 1 {
            return Ok(self.clone(coeff_cache));
        }
        let mut res = Self::one(coeff_cache);
        let max_deg = poly_props.semigroup.max_degree;
        let min_deg = self.min_degree(poly_props);
        let invert = n < 0;
        let mut n = n.abs() as u32;
        let mut tmp_poly = if invert {
            self.truncated(max_deg - (n - 1) * min_deg, poly_props, coeff_cache)
        } else {
            let tmp_poly2 = self.recipr(&poly_props, coeff_cache)?;
            let tmp_poly3 =
                tmp_poly2.truncated(max_deg - (n - 1) * min_deg, poly_props, coeff_cache);
            tmp_poly2.drop(coeff_cache);
            tmp_poly3
        };
        loop {
            if n % 2 != 0 {
                let tmp_poly2 = res.mul(&tmp_poly, poly_props, coeff_cache);
                tmp_poly2.move_into(&mut res, coeff_cache);
            }
            n >>= 1;
            if n == 0 {
                break;
            }
            let tmp_poly2 = tmp_poly.mul(&tmp_poly, poly_props, coeff_cache);
            tmp_poly2.move_into(&mut tmp_poly, coeff_cache);
        }
        tmp_poly.drop(coeff_cache);
        Ok(res)
    }

    /// For a polynomial $p$, compute $\exp(p)$ and $\exp(-p)$.
    pub fn exp_pos_neg(
        &self,
        poly_props: &PolynomialProperties<T>,
        coeff_cache: &mut NumCache<T>,
    ) -> Result<(Self, Self), PolynomialError> {
        if self.coeffs.get(&0).is_some() {
            return Err(PolynomialError::NonZeroConstantTermError);
        }
        let mut res_pos = Self::one(coeff_cache);
        let mut res_neg = Self::one(coeff_cache);
        let mut tmp_poly = Self::one(coeff_cache);
        let min_deg = self.min_degree(poly_props);
        let max_deg = poly_props.semigroup.max_degree;
        let mut tmp_var = coeff_cache.pop();
        tmp_var.assign(1);
        for i in 1..=max_deg / min_deg {
            let tmp_poly2 = self.mul(&tmp_poly, poly_props, coeff_cache);
            tmp_poly2.move_into(&mut tmp_poly, coeff_cache);
            tmp_var *= i;
            res_pos.div_scalar_add_assign(&tmp_var, &tmp_poly, coeff_cache);
            if i % 2 != 0 {
                tmp_var *= -1;
            }
            res_neg.div_scalar_add_assign(&tmp_var, &tmp_poly, coeff_cache);
            if i % 2 != 0 {
                tmp_var *= -1;
            }
        }
        tmp_poly.drop(coeff_cache);
        Ok((res_pos, res_neg))
    }

    /// Compute the dilogarithm of a polynomial.
    pub fn li_2(
        &self,
        poly_props: &PolynomialProperties<T>,
        coeff_cache: &mut NumCache<T>,
    ) -> Result<Self, PolynomialError> {
        if self.coeffs.get(&0).is_some() {
            return Err(PolynomialError::NonZeroConstantTermError);
        }
        let mut res = self.clone(coeff_cache);
        let mut tmp_poly = self.clone(coeff_cache);
        let min_deg = self.min_degree(poly_props);
        let max_deg = poly_props.semigroup.max_degree;
        for i in 2..=max_deg / min_deg {
            let tmp_poly2 = self.mul(&tmp_poly, poly_props, coeff_cache);
            tmp_poly2.move_into(&mut tmp_poly, coeff_cache);
            res.div_scalar_add_assign(i * i, &tmp_poly, coeff_cache);
        }
        tmp_poly.drop(coeff_cache);
        Ok(res)
    }
}

#[macro_export]
macro_rules! polynomial {
    ( $t:expr, $( ($k:expr, $v:expr) ),* ) => {
        {
            let mut nonzero = Vec::new();
            let mut coeffs = HashMap::new();
            $(
                let mut c = $t.clone();
                c.assign($v);
                coeffs.insert($k, c);
                nonzero.push($k);
            )*
            Polynomial { coeffs, nonzero }}
    };
}

#[cfg(test)]
mod tests {
    use crate::semigroup;

    use super::*;
    use nalgebra::{DMatrix, RowDVector};
    use rug::Rational;
    use semigroup::Semigroup;

    // Create a simple semigroup for testing.
    fn example_semigroup() -> Semigroup {
        #[rustfmt::skip]
        let exponents = DMatrix::from_column_slice(2, 6,
            &[
                0, 0, // 1
                1, 0, // x
                0, 1, // y
                2, 0, // x^2
                1, 1, // x*y
                0, 2  // y^2
            ]
        );
        let grading_vector = RowDVector::from_row_slice(&[1, 1]);
        Semigroup::from_data(exponents, grading_vector).unwrap()
    }

    #[test]
    fn test_polynomial() {
        let tmp_rational = Rational::new();
        let semigroup = example_semigroup();
        let poly_props = PolynomialProperties::new(&semigroup, &tmp_rational);

        let mut coeff_cache = NumCache::new(tmp_rational.clone(), 100);

        #[rustfmt::skip]
        let p1 = polynomial!(
            tmp_rational,
            (0, 1),    // 1
            (1, 2),    // 2*x
            (2, 3),    // 3*y
            (3, 4),    // 4*x^2
            (4, 5),    // 5*x*y
            (5, 6)     // 6*y^2
        );

        // Test multiplication
        let p2 = p1.mul(&p1, &poly_props, &mut coeff_cache);
        let p2_res = polynomial!(
            tmp_rational,
            (0, 1),  // 1
            (1, 4),  // 4*x
            (2, 6),  // 6*y
            (3, 12), // 12*x^2
            (4, 22), // 22*x*y
            (5, 21)  // 21*y^2
        );
        assert_eq!(p2, p2_res);

        // TODO: add more tests
    }
}
