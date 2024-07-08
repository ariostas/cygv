//! Polynomials containing a fixed set of monomials.

pub mod coefficient;
pub mod error;
pub mod prettyprint;
pub mod properties;

use coefficient::PolynomialCoeff;
use core::ops::{DivAssign, MulAssign};
use error::PolynomialError;
use properties::PolynomialProperties;
use std::collections::HashMap;

/// A polynomial structure.
///
/// The structure contains a hash map from monomial indices to the value of
/// the coefficient. It also contains a vector listing the sorted indices of
/// its nonzero coefficients. Even though this is redundant information, it is
/// useful to have this list to make some operations more efficient.
#[derive(Debug, PartialEq)]
pub struct Polynomial<T> {
    pub coeffs: HashMap<usize, T>,
    pub nonzero: Vec<usize>,
}

impl<T> Default for Polynomial<T>
where
    T: PolynomialCoeff<T>,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T> Clone for Polynomial<T>
where
    T: PolynomialCoeff<T>,
{
    fn clone(&self) -> Self {
        Self {
            coeffs: self.coeffs.clone(),
            nonzero: self.nonzero.clone(),
        }
    }
}

impl<T> Polynomial<T>
where
    T: PolynomialCoeff<T>,
{
    /// Create a new polynomial.
    pub fn new() -> Self {
        Self {
            coeffs: HashMap::new(),
            nonzero: Vec::new(),
        }
    }

    /// Create a new polynomial equal to zero.
    pub fn zero() -> Self {
        Self::new()
    }

    /// Create a new polynomial with one as its constant term.
    pub fn one(poly_props: &PolynomialProperties<T>) -> Self {
        let mut coeffs = HashMap::new();
        let one = poly_props.one();
        coeffs.insert(0, one);
        Self {
            coeffs,
            nonzero: vec![0],
        }
    }

    /// Clear the polynomial, while moving the allocated coefficients to the
    /// pool.
    pub fn clear(&mut self) {
        self.coeffs.clear();
        self.nonzero.clear();
    }

    /// Clean up the polynomial, removing monomials whose coefficient are
    /// below a given cutoff.
    pub fn clean_up(&mut self, poly_props: &PolynomialProperties<T>) {
        let mut new_nonzero = Vec::new();
        let mut to_delete = Vec::new();
        let pos_cutoff = poly_props.zero_cutoff.clone();
        let neg_cutoff = -poly_props.zero_cutoff.clone();
        for i in self.nonzero.drain(..) {
            let c = self.coeffs.get(&i).unwrap();
            if *c <= pos_cutoff && *c >= neg_cutoff {
                to_delete.push(i);
            } else {
                new_nonzero.push(i);
            }
        }
        for i in to_delete.iter() {
            self.coeffs.remove(i);
        }
        self.nonzero = new_nonzero;
    }

    /// Find the minimum degree of the polynomial.
    pub fn min_degree(&self, poly_props: &PolynomialProperties<T>) -> u32 {
        let Some(idx) = self.nonzero.first() else {
            return poly_props.semigroup.max_degree + 1;
        };
        poly_props.semigroup.degrees[*idx]
    }

    /// Clone the polynomial, but truncated to a given degree.
    pub fn truncated(&self, max_deg: u32, poly_props: &PolynomialProperties<T>) -> Self {
        let mut res = Self::new();
        for &i in self.nonzero.iter() {
            if poly_props.semigroup.degrees[i] > max_deg {
                break;
            }
            res.nonzero.push(i);
            let coeff = self.coeffs.get(&i).unwrap().clone();
            res.coeffs.insert(i, coeff);
        }
        res
    }

    /// Add polynomials in place.
    pub fn add_assign(&mut self, rhs: &Self, poly_props: &PolynomialProperties<T>) {
        let mut resort = false;
        for (k, v) in rhs.coeffs.iter() {
            let c = self.coeffs.entry(*k).or_insert_with(|| {
                resort = true;
                poly_props.zero()
            });
            *c += v;
        }
        if resort {
            self.nonzero = self.coeffs.keys().cloned().collect();
            self.nonzero.sort_unstable();
        }
    }

    /// Subtract polynomials in place.
    pub fn sub_assign(&mut self, rhs: &Self, poly_props: &PolynomialProperties<T>) {
        let mut resort = false;
        for (k, v) in rhs.coeffs.iter() {
            let c = self.coeffs.entry(*k).or_insert_with(|| {
                resort = true;
                poly_props.zero()
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
    pub fn mul(&self, rhs: &Self, poly_props: &PolynomialProperties<T>) -> Self {
        let mut res = Self::new();
        let max_deg = poly_props.semigroup.max_degree;
        let (pshort, plong) = if self.nonzero.len() < rhs.nonzero.len() {
            (self, rhs)
        } else {
            (rhs, self)
        };
        for &i in pshort.nonzero.iter() {
            let deg1 = poly_props.semigroup.degrees[i];
            for &j in plong.nonzero.iter() {
                let deg2 = poly_props.semigroup.degrees[j];
                if deg1 + deg2 > max_deg {
                    break;
                }
                let tmp_vec = poly_props.semigroup.elements.column(i)
                    + poly_props.semigroup.elements.column(j);
                let Some(mon) = poly_props.monomial_map.get(&tmp_vec.as_view()) else {
                    continue;
                };
                let c = res.coeffs.entry(*mon).or_insert_with(|| {
                    res.nonzero.push(*mon);
                    poly_props.zero()
                });
                let mut tmp_var = pshort.coeffs.get(&i).unwrap().clone();
                tmp_var *= plong.coeffs.get(&j).unwrap();
                *c += &tmp_var;
            }
        }
        res.nonzero.sort_unstable();
        res
    }

    /// Multiply polynomial by a scalar, and then add in place.
    pub fn mul_scalar_add_assign<U>(
        &mut self,
        scalar: U,
        rhs: &Self,
        poly_props: &PolynomialProperties<T>,
    ) where
        T: MulAssign<U>,
        U: Copy,
    {
        let mut resort = false;
        for (k, v) in rhs.coeffs.iter() {
            let c = self.coeffs.entry(*k).or_insert_with(|| {
                resort = true;
                poly_props.zero()
            });
            let mut tmp_var = v.clone();
            tmp_var *= scalar;
            *c += &tmp_var;
        }
        if resort {
            self.nonzero = self.coeffs.keys().cloned().collect();
            self.nonzero.sort_unstable();
        }
    }

    /// Divide polynomial by a scalar, and then add in place.
    pub fn div_scalar_add_assign<U>(
        &mut self,
        scalar: U,
        rhs: &Self,
        poly_props: &PolynomialProperties<T>,
    ) where
        T: DivAssign<U>,
        U: Copy,
    {
        let mut resort = false;
        for (k, v) in rhs.coeffs.iter() {
            let c = self.coeffs.entry(*k).or_insert_with(|| {
                resort = true;
                poly_props.zero()
            });
            let mut tmp_var = v.clone();
            tmp_var /= scalar;
            *c += &tmp_var;
        }
        if resort {
            self.nonzero = self.coeffs.keys().cloned().collect();
            self.nonzero.sort_unstable();
        }
    }

    /// Compute the reciprocal of the polynomial.
    ///
    /// Only works for polynomials with a nonzero constant term.
    pub fn recipr(&self, poly_props: &PolynomialProperties<T>) -> Result<Self, PolynomialError> {
        let Some(c) = self.coeffs.get(&0) else {
            return Err(PolynomialError::ZeroConstantTermError);
        };
        if *c == 0 {
            return Err(PolynomialError::ZeroConstantTermError);
        }
        let mut res = Self::one(poly_props);
        let max_deg = poly_props.semigroup.max_degree;
        let mut p0 = self.clone();
        p0.nonzero.remove(0);
        let const_term = p0.coeffs.remove(&0).unwrap();
        p0.div_scalar_assign(&const_term);
        let mut tmp_poly = Self::one(poly_props);
        let min_deg = p0.min_degree(poly_props);
        for i in 1..=max_deg / min_deg {
            tmp_poly = tmp_poly.mul(&p0, poly_props);
            if i % 2 == 0 {
                res.add_assign(&tmp_poly, poly_props)
            } else {
                res.mul_scalar_add_assign(-1, &tmp_poly, poly_props);
            }
        }
        res.div_scalar_assign(&const_term);
        Ok(res)
    }

    /// Compute powers of polynomials.
    pub fn pow(
        &self,
        n: i32,
        poly_props: &PolynomialProperties<T>,
    ) -> Result<Self, PolynomialError> {
        if n == 0 {
            return Ok(Self::one(poly_props));
        } else if n == 1 {
            return Ok(self.clone());
        }
        let mut res = Self::one(poly_props);
        let max_deg = poly_props.semigroup.max_degree;
        let min_deg = self.min_degree(poly_props);
        let invert = n < 0;
        let mut n = n.unsigned_abs();
        let mut tmp_poly = if invert {
            let tmp_poly2 = self.recipr(poly_props)?;
            tmp_poly2.truncated(max_deg - (n - 1) * min_deg, poly_props)
        } else {
            self.truncated(max_deg - (n - 1) * min_deg, poly_props)
        };
        loop {
            if n % 2 != 0 {
                res = res.mul(&tmp_poly, poly_props);
            }
            n >>= 1;
            if n == 0 {
                break;
            }
            tmp_poly = tmp_poly.mul(&tmp_poly, poly_props);
        }
        Ok(res)
    }

    /// For a polynomial $p$, compute $\exp(p)$ and $\exp(-p)$.
    pub fn exp_pos_neg(
        &self,
        poly_props: &PolynomialProperties<T>,
    ) -> Result<(Self, Self), PolynomialError> {
        if self.coeffs.contains_key(&0) {
            return Err(PolynomialError::NonZeroConstantTermError);
        }
        let mut res_pos = Self::one(poly_props);
        let mut res_neg = Self::one(poly_props);
        let mut tmp_poly = Self::one(poly_props);
        let min_deg = self.min_degree(poly_props);
        let max_deg = poly_props.semigroup.max_degree;
        let mut tmp_var = poly_props.one();
        for i in 1..=max_deg / min_deg {
            tmp_poly = self.mul(&tmp_poly, poly_props);
            tmp_var *= i;
            res_pos.div_scalar_add_assign(&tmp_var, &tmp_poly, poly_props);
            if i % 2 != 0 {
                tmp_var *= -1;
            }
            res_neg.div_scalar_add_assign(&tmp_var, &tmp_poly, poly_props);
            if i % 2 != 0 {
                tmp_var *= -1;
            }
        }
        Ok((res_pos, res_neg))
    }

    /// Compute the dilogarithm of a polynomial.
    pub fn li_2(&self, poly_props: &PolynomialProperties<T>) -> Result<Self, PolynomialError> {
        if self.coeffs.contains_key(&0) {
            return Err(PolynomialError::NonZeroConstantTermError);
        }
        let mut res = self.clone();
        let mut tmp_poly = self.clone();
        let min_deg = self.min_degree(poly_props);
        let max_deg = poly_props.semigroup.max_degree;
        for i in 2..=max_deg / min_deg {
            tmp_poly = self.mul(&tmp_poly, poly_props);
            res.div_scalar_add_assign(i * i, &tmp_poly, poly_props);
        }
        Ok(res)
    }
}

#[macro_export]
/// Construct a polynomial directly from data.
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
    use super::*;
    use crate::semigroup;
    use nalgebra::{DMatrix, RowDVector};
    use rug::{Assign, Rational};
    use semigroup::Semigroup;

    // Construct simple data for testing
    fn example_semigroup() -> Semigroup {
        let elements = DMatrix::from_column_slice(
            2,
            6,
            &[
                0, 0, // 1
                1, 0, // x
                0, 1, // y
                2, 0, // x^2
                1, 1, // x*y
                0, 2, // y^2
            ],
        );
        let grading_vector = RowDVector::from_row_slice(&[1, 1]);
        Semigroup::from_data(elements, grading_vector).unwrap()
    }

    #[test]
    fn test_identity() {
        let tmp_rational = Rational::new();
        let semigroup = example_semigroup();
        let poly_props = PolynomialProperties::new(&semigroup, &tmp_rational);

        let p = Polynomial::one(&poly_props);
        let p_res = polynomial!(
            tmp_rational,
            (0, 1) // 1
        );
        assert_eq!(p, p_res);
    }

    #[test]
    fn test_cleanup() {
        let tmp_rational = Rational::new();
        let semigroup = example_semigroup();
        let poly_props = PolynomialProperties::new(&semigroup, &tmp_rational);

        let mut p = polynomial!(
            tmp_rational,
            (0, 1), // 1
            (1, 0), // 0*x
            (2, 1), // 1*y
            (3, 0), // 0*x^2
            (4, 1), // 1*x*y
            (5, 0)  // 0*y^2
        );
        let p_res = polynomial!(
            tmp_rational,
            (0, 1), // 1
            (2, 1), // 1*y
            (4, 1)  // 1*x*y
        );
        p.clean_up(&poly_props);
        assert_eq!(p, p_res);
    }

    #[test]
    fn test_min_degree() {
        let tmp_rational = Rational::new();
        let semigroup = example_semigroup();
        let poly_props = PolynomialProperties::new(&semigroup, &tmp_rational);

        let p0 = Polynomial::one(&poly_props);
        let p1 = polynomial!(
            tmp_rational,
            (2, 3), // 3*y
            (3, 4), // 4*x^2
            (4, 5), // 5*x*y
            (5, 6)  // 6*y^2
        );
        let p2 = polynomial!(
            tmp_rational,
            (5, 6) // 6*y^2
        );
        let p3: Polynomial<Rational> = Polynomial::new();

        assert_eq!(p0.min_degree(&poly_props), 0);
        assert_eq!(p1.min_degree(&poly_props), 1);
        assert_eq!(p2.min_degree(&poly_props), 2);
        assert_eq!(p3.min_degree(&poly_props), 3);
    }

    #[test]
    fn test_truncated() {
        let tmp_rational = Rational::new();
        let semigroup = example_semigroup();
        let poly_props = PolynomialProperties::new(&semigroup, &tmp_rational);

        let p = polynomial!(
            tmp_rational,
            (0, 1), // 1
            (1, 0), // 0*x
            (2, 1), // 1*y
            (3, 0), // 0*x^2
            (4, 1), // 1*x*y
            (5, 0)  // 0*y^2
        );
        let p_res = polynomial!(
            tmp_rational,
            (0, 1), // 1
            (1, 0), // 0*x
            (2, 1)  // 1*y
        );
        let p = p.truncated(1, &poly_props);
        assert_eq!(p, p_res);
    }

    #[test]
    fn test_add_assign() {
        let tmp_rational = Rational::new();
        let semigroup = example_semigroup();
        let poly_props = PolynomialProperties::new(&semigroup, &tmp_rational);

        let mut p = polynomial!(
            tmp_rational,
            (0, 1), // 1
            (1, 2), // 2*x
            (2, 3), // 3*y
            (3, 4), // 4*x^2
            (4, 5), // 5*x*y
            (5, 6)  // 6*y^2
        );
        let p1 = p.clone();
        p.add_assign(&p1, &poly_props);
        let p_res = polynomial!(
            tmp_rational,
            (0, 2),  // 2
            (1, 4),  // 4*x
            (2, 6),  // 6*y
            (3, 8),  // 8*x^2
            (4, 10), // 10*x*y
            (5, 12)  // 12*y^2
        );
        assert_eq!(p, p_res);
    }

    #[test]
    fn test_sub_assign() {
        let tmp_rational = Rational::new();
        let semigroup = example_semigroup();
        let poly_props = PolynomialProperties::new(&semigroup, &tmp_rational);

        let mut p = polynomial!(
            tmp_rational,
            (0, 1), // 1
            (1, 2), // 2*x
            (2, 3), // 3*y
            (3, 4), // 4*x^2
            (4, 5), // 5*x*y
            (5, 6)  // 6*y^2
        );
        let p1 = polynomial!(
            tmp_rational,
            (0, 2),  // 2
            (1, 4),  // 4*x
            (2, 6),  // 6*y
            (3, 8),  // 8*x^2
            (4, 10), // 10*x*y
            (5, 12)  // 12*y^2
        );
        p.sub_assign(&p1, &poly_props);
        let p_res = polynomial!(
            tmp_rational,
            (0, -1), // -1
            (1, -2), // -2*x
            (2, -3), // -3*y
            (3, -4), // -4*x^2
            (4, -5), // -5*x*y
            (5, -6)  // -6*y^2
        );
        assert_eq!(p, p_res);
    }

    #[test]
    fn test_mul_scalar_assign() {
        let tmp_rational = Rational::from((2, 1));

        let mut p = polynomial!(
            tmp_rational,
            (0, 1), // 1
            (1, 2), // 2*x
            (2, 3), // 3*y
            (3, 4), // 4*x^2
            (4, 5), // 5*x*y
            (5, 6)  // 6*y^2
        );
        p.mul_scalar_assign(&tmp_rational);
        let p_res = polynomial!(
            tmp_rational,
            (0, 2),  // 2
            (1, 4),  // 4*x
            (2, 6),  // 6*y
            (3, 8),  // 8*x^2
            (4, 10), // 10*x*y
            (5, 12)  // 12*y^2
        );
        assert_eq!(p, p_res);
    }

    #[test]
    fn test_div_scalar_assign() {
        let tmp_rational = Rational::from((1, 2));

        let mut p = polynomial!(
            tmp_rational,
            (0, 1), // 1
            (1, 2), // 2*x
            (2, 3), // 3*y
            (3, 4), // 4*x^2
            (4, 5), // 5*x*y
            (5, 6)  // 6*y^2
        );
        p.div_scalar_assign(&tmp_rational);
        let p_res = polynomial!(
            tmp_rational,
            (0, 2),  // 2
            (1, 4),  // 4*x
            (2, 6),  // 6*y
            (3, 8),  // 8*x^2
            (4, 10), // 10*x*y
            (5, 12)  // 12*y^2
        );
        assert_eq!(p, p_res);
    }

    #[test]
    fn test_mul() {
        let tmp_rational = Rational::new();
        let semigroup = example_semigroup();
        let poly_props = PolynomialProperties::new(&semigroup, &tmp_rational);

        let p1 = polynomial!(
            tmp_rational,
            (0, 1), // 1
            (1, 2), // 2*x
            (2, 3), // 3*y
            (3, 4), // 4*x^2
            (4, 5), // 5*x*y
            (5, 6)  // 6*y^2
        );
        let p = p1.mul(&p1, &poly_props);
        let p_res = polynomial!(
            tmp_rational,
            (0, 1),  // 1
            (1, 4),  // 4*x
            (2, 6),  // 6*y
            (3, 12), // 12*x^2
            (4, 22), // 22*x*y
            (5, 21)  // 21*y^2
        );
        assert_eq!(p, p_res);
    }

    #[test]
    fn test_mul_scalar_add_assign() {
        let mut tmp_rational = Rational::new();
        let semigroup = example_semigroup();
        let poly_props = PolynomialProperties::new(&semigroup, &tmp_rational);
        tmp_rational.assign(2);

        let mut p = polynomial!(
            tmp_rational,
            (0, 1), // 1
            (1, 2), // 2*x
            (2, 3), // 3*y
            (3, 4), // 4*x^2
            (4, 5), // 5*x*y
            (5, 6)  // 6*y^2
        );
        let p1 = p.clone();
        p.mul_scalar_add_assign(&tmp_rational, &p1, &poly_props);
        let p_res = polynomial!(
            tmp_rational,
            (0, 3),  // 3
            (1, 6),  // 6*x
            (2, 9),  // 9*y
            (3, 12), // 12*x^2
            (4, 15), // 15*x*y
            (5, 18)  // 18*y^2
        );
        assert_eq!(p, p_res);
    }

    #[test]
    fn test_div_scalar_add_assign() {
        let mut tmp_rational = Rational::new();
        let semigroup = example_semigroup();
        let poly_props = PolynomialProperties::new(&semigroup, &tmp_rational);
        tmp_rational.assign((1, 2));

        let mut p = polynomial!(
            tmp_rational,
            (0, 1), // 1
            (1, 2), // 2*x
            (2, 3), // 3*y
            (3, 4), // 4*x^2
            (4, 5), // 5*x*y
            (5, 6)  // 6*y^2
        );
        let p1 = p.clone();
        p.div_scalar_add_assign(&tmp_rational, &p1, &poly_props);
        let p_res = polynomial!(
            tmp_rational,
            (0, 3),  // 3
            (1, 6),  // 6*x
            (2, 9),  // 9*y
            (3, 12), // 12*x^2
            (4, 15), // 15*x*y
            (5, 18)  // 18*y^2
        );
        assert_eq!(p, p_res);
    }

    #[test]
    fn test_recipr() {
        let tmp_rational = Rational::new();
        let semigroup = example_semigroup();
        let poly_props = PolynomialProperties::new(&semigroup, &tmp_rational);

        let p1 = polynomial!(
            tmp_rational,
            (0, 1), // 1
            (1, 2), // 2*x
            (2, 3), // 3*y
            (3, 4), // 4*x^2
            (4, 5), // 5*x*y
            (5, 6)  // 6*y^2
        );
        let mut p = p1.recipr(&poly_props).unwrap();
        p.clean_up(&poly_props);
        let p_res = polynomial!(
            tmp_rational,
            (0, 1),  // 1
            (1, -2), // -2*x
            (2, -3), // -3*y
            (4, 7),  // 7*x*y
            (5, 3)   // 3*y^2
        );
        assert_eq!(p, p_res);
    }

    #[test]
    fn test_pow() {
        let tmp_rational = Rational::new();
        let semigroup = example_semigroup();
        let poly_props = PolynomialProperties::new(&semigroup, &tmp_rational);

        let p = polynomial!(
            tmp_rational,
            (0, 1), // 1
            (1, 2), // 2*x
            (2, 3), // 3*y
            (3, 4), // 4*x^2
            (4, 5), // 5*x*y
            (5, 6)  // 6*y^2
        );

        let p_0 = p.pow(0, &poly_props).unwrap();
        let p_0_res = Polynomial::one(&poly_props);
        assert_eq!(p_0, p_0_res);

        let p_m1 = p.pow(-1, &poly_props).unwrap();
        let p_m1_res = p.recipr(&poly_props).unwrap();
        assert_eq!(p_m1, p_m1_res);

        let p_1 = p.pow(1, &poly_props).unwrap();
        assert_eq!(p, p_1);

        let p_2 = p.pow(2, &poly_props).unwrap();
        let p_2_res = p.mul(&p, &poly_props);
        assert_eq!(p_2, p_2_res);

        let p_3 = p.pow(3, &poly_props).unwrap();
        let p_3_res = p.mul(&p_2, &poly_props);
        assert_eq!(p_3, p_3_res);
    }

    #[test]
    fn test_exp_pos_neg() {
        let tmp_rational = Rational::new();
        let semigroup = example_semigroup();
        let poly_props = PolynomialProperties::new(&semigroup, &tmp_rational);

        let p = polynomial!(
            tmp_rational,
            (1, 2), // 2*x
            (2, 3), // 3*y
            (3, 4), // 4*x^2
            (4, 5), // 5*x*y
            (5, 6)  // 6*y^2
        );

        let (p_exp_pos, p_exp_neg) = p.exp_pos_neg(&poly_props).unwrap();
        let p_exp_pos_res = polynomial!(
            tmp_rational,
            (0, 1),       // 1
            (1, 2),       // 2*x
            (2, 3),       // 3*y
            (3, 6),       // 6*x^2
            (4, 11),      // 11*x*y
            (5, (21, 2))  // (21/2)*y^2
        );
        let p_exp_neg_res = polynomial!(
            tmp_rational,
            (0, 1),       // 1
            (1, -2),      // -2*x
            (2, -3),      // -3*y
            (3, -2),      // -2*x^2
            (4, 1),       // 1*x*y
            (5, (-3, 2))  // (-3/2)*y^2
        );
        assert_eq!(p_exp_pos, p_exp_pos_res);
        assert_eq!(p_exp_neg, p_exp_neg_res);
    }

    #[test]
    fn test_li_2() {
        let tmp_rational = Rational::new();
        let semigroup = example_semigroup();
        let poly_props = PolynomialProperties::new(&semigroup, &tmp_rational);

        let p = polynomial!(
            tmp_rational,
            (1, 2), // 2*x
            (2, 3), // 3*y
            (3, 4), // 4*x^2
            (4, 5), // 5*x*y
            (5, 6)  // 6*y^2
        );

        let p_li_2 = p.li_2(&poly_props).unwrap();
        let p_li_2_res = polynomial!(
            tmp_rational,
            (1, 2),       // 2*x
            (2, 3),       // 3*y
            (3, 5),       // 5*x^2
            (4, 8),       // 8*x*y
            (5, (33, 4))  // (33/4)*y^2
        );
        assert_eq!(p_li_2, p_li_2_res);
    }
}
