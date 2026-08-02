//! Pretty-print polynomials.

use super::Polynomial;
use super::PolynomialProperties;
use core::fmt;

/// A pretty-print struct for polynomials.
#[derive(Debug, Clone)]
pub struct PrettyPrintPolynomial<'a, T> {
    pub polynomial: &'a Polynomial<T>,
    pub properties: &'a PolynomialProperties<'a, T>,
}

/// Implement Display trait for PrettyPrintPolynomial.
impl<'a, T: fmt::Display> fmt::Display for PrettyPrintPolynomial<'a, T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let variables = ["x", "y", "z", "w", "v", "u"];
        let n_variables = self.properties.semigroup.elements.nrows();
        let mut message = String::new();
        for (i_ind, i) in self.polynomial.nonzero.iter().enumerate() {
            let coeff = self.polynomial.coeffs.get(i).unwrap();
            if i_ind != 0 {
                message += " + ";
            }
            message += &format!("{coeff}");
            for (v_ind, v) in self
                .properties
                .semigroup
                .elements
                .column(*i)
                .iter()
                .enumerate()
            {
                if *v == 0 {
                    continue;
                }
                if n_variables <= 6 {
                    message += &format!("*{}", variables[v_ind]);
                } else {
                    message += &format!("*x_{}", v_ind);
                }
                if *v != 1 {
                    message += &format!("^{}", v);
                }
            }
        }
        write!(f, "{message}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pool::NumberPool;
    use crate::semigroup::Semigroup;
    use nalgebra::{DMatrix, DVector, RowDVector};
    use rug::{Assign, Rational};

    #[test]
    fn test_display_named_variables() {
        // Two variables, but more than six monomials: the variables must still be
        // printed with their names, since the count is the number of rows.
        #[rustfmt::skip]
        let elements = DMatrix::from_column_slice(2, 7,
            &[
                0, 0,
                1, 0,
                0, 1,
                2, 0,
                1, 1,
                0, 2,
                3, 0,
            ]
        );
        let grading_vector = RowDVector::from_row_slice(&[1, 1]);
        let sg = Semigroup::from_data(elements, grading_vector).unwrap();
        let tmp_rational = Rational::new();
        let poly_props = PolynomialProperties::new(&sg, &tmp_rational);
        let mut coeff_pool = NumberPool::new(tmp_rational.clone(), 10);

        let index_of = |monomial: &[i32]| {
            let monomial = DVector::from_column_slice(monomial);
            poly_props.monomial_map[&monomial.as_view()]
        };
        let display = |p: &Polynomial<Rational>| {
            format!(
                "{}",
                PrettyPrintPolynomial {
                    polynomial: p,
                    properties: &poly_props,
                }
            )
        };

        let mut p = Polynomial::one(&mut coeff_pool);
        assert_eq!(display(&p), "1");

        let mut coeff = coeff_pool.pop();
        coeff.assign(2);
        p.coeffs.insert(index_of(&[1, 1]), coeff);
        p.nonzero.push(index_of(&[1, 1]));
        assert_eq!(display(&p), "1 + 2*x*y");

        let mut coeff = coeff_pool.pop();
        coeff.assign(3);
        p.coeffs.insert(index_of(&[3, 0]), coeff);
        p.nonzero.push(index_of(&[3, 0]));
        assert_eq!(display(&p), "1 + 2*x*y + 3*x^3");
    }
}
