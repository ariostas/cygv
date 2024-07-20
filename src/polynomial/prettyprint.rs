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
        let n_variables = self.properties.semigroup.elements.ncols();
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
