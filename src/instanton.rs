//! Instanton corrections.

use crate::fundamental_period::FundamentalPeriod;
use crate::polynomial::{
    coefficient::PolynomialCoeff, error::PolynomialError, properties::PolynomialProperties,
    Polynomial,
};
use crate::CYKind;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

pub struct InstantonData<T> {
    pub inst: Vec<Polynomial<T>>,
    pub expalpha: Vec<(Polynomial<T>, Polynomial<T>)>,
}

/// Computes one instanton correction, i.e. one entry of `InstantonData::inst`.
fn compute_inst<T>(
    t: usize,
    f_poly: &HashMap<(usize, usize), Polynomial<T>>,
    poly_props: &PolynomialProperties<T>,
    intnum_dict: &HashMap<(usize, usize, usize), i32>,
    cy_kind: CYKind,
) -> Polynomial<T>
where
    T: PolynomialCoeff<T>,
{
    let h11 = poly_props.semigroup.elements.nrows();
    let mut intnum_ind = [0_usize; 3];
    let mut tmp_num = poly_props.zero.clone();
    let mut p = Polynomial::new();
    for a in 0..h11 {
        for b in a..h11 {
            intnum_ind[0] = t;
            intnum_ind[1] = a;
            intnum_ind[2] = b;
            if cy_kind.is_threefold() {
                intnum_ind.sort_unstable();
            }
            let Some(x) = intnum_dict.get(&(intnum_ind[0], intnum_ind[1], intnum_ind[2])) else {
                continue;
            };
            let mut tmp_poly = f_poly[&(a, b)].clone(&poly_props.zero);
            if a != b {
                tmp_poly.mul_scalar_assign(*x);
            } else {
                tmp_num.assign(*x);
                tmp_num /= 2;
                tmp_poly.mul_scalar_assign(&tmp_num);
            }
            p.add_assign(&tmp_poly, &poly_props.zero);
        }
    }
    p.clean_up(poly_props);
    p
}

/// Compute instanton corrections, as well as other objects needed for the series inversion.
pub fn compute_instanton_data<T>(
    fp: FundamentalPeriod<T>,
    poly_props: &PolynomialProperties<T>,
    intnum_idxpairs: &HashSet<(usize, usize)>,
    n_indices: usize,
    intnum_dict: &HashMap<(usize, usize, usize), i32>,
    cy_kind: CYKind,
) -> Result<InstantonData<T>, PolynomialError>
where
    T: PolynomialCoeff<T>,
{
    let h11 = poly_props.semigroup.elements.nrows();

    // Compute alpha polynomials
    let alpha: Vec<Polynomial<T>> = (0..h11)
        .into_par_iter()
        .map(|t| {
            let mut a = fp.c0_inv.mul(&fp.c1[t], poly_props);
            a.clean_up(poly_props);
            a
        })
        .collect();

    // Compute beta polynomials
    let beta: HashMap<(usize, usize), Polynomial<T>> = intnum_idxpairs
        .par_iter()
        .map(|&(t0, t1)| {
            let mut a = fp.c0_inv.mul(&fp.c2[&(t0, t1)], poly_props);
            a.clean_up(poly_props);
            ((t0, t1), a)
        })
        .collect();

    // Compute F polynomials
    let f_poly: HashMap<(usize, usize), Polynomial<T>> = intnum_idxpairs
        .par_iter()
        .map(|&(t0, t1)| {
            let mut p = alpha[t0].mul(&alpha[t1], poly_props);
            p.sub_assign(&beta[&(t0, t1)], &poly_props.zero);
            p.mul_scalar_assign(-1);
            p.clean_up(poly_props);
            ((t0, t1), p)
        })
        .collect();

    // Compute instanton corrections
    let inst: Vec<Polynomial<T>> = (0..n_indices)
        .into_par_iter()
        .map(|t| compute_inst(t, &f_poly, poly_props, intnum_dict, cy_kind))
        .collect();

    // Compute expalpha polynomials. Collecting into a `Result` stops the
    // remaining curves as soon as one of them fails.
    let expalpha: Vec<(Polynomial<T>, Polynomial<T>)> = (0..h11)
        .into_par_iter()
        .map(|t| {
            let mut p = alpha[t].exp_pos_neg(poly_props)?;
            p.0.clean_up(poly_props);
            p.1.clean_up(poly_props);
            Ok(p)
        })
        .collect::<Result<_, PolynomialError>>()?;

    Ok(InstantonData { inst, expalpha })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{fundamental_period::compute_omega, misc::process_int_nums, Semigroup};
    use nalgebra::{DMatrix, RowDVector};
    use rug::Rational;

    #[test]
    fn test_instanton() {
        let generators = DMatrix::from_column_slice(2, 2, &[0, -1, 1, 2]);
        let grading_vector = RowDVector::from_row_slice(&[3, -1]);
        let sg = Semigroup::with_min_elements(generators, grading_vector, 15).unwrap();
        let zero_rat = Rational::new();
        let poly_props = PolynomialProperties::new(&sg, &zero_rat);

        let q = DMatrix::from_column_slice(6, 2, &[1, 1, 1, 0, 1, 2, 0, 0, -1, 1, 1, -1]);
        let nefpart = Vec::new();
        let intnum_idxpairs = [(0, 0), (0, 1), (1, 1)].iter().cloned().collect();

        let fp = compute_omega(&poly_props, &sg, &q, &nefpart, &intnum_idxpairs);
        assert!(fp.is_ok());
        let fp = fp.unwrap();

        let intnums = HashMap::from([
            ((0, 0, 0), 2),
            ((0, 0, 1), 1),
            ((0, 1, 1), -1),
            ((1, 1, 1), -5),
        ]);
        let result = process_int_nums(intnums.clone(), CYKind::Threefold);
        assert!(result.is_ok());
        let (intnum_dict, intnum_idxpairs, n_indices) = result.unwrap();

        let inst_data = compute_instanton_data(
            fp,
            &poly_props,
            &intnum_idxpairs,
            n_indices,
            &intnum_dict,
            CYKind::Threefold,
        );
        assert!(inst_data.is_ok());
        let inst_data = inst_data.unwrap();

        let inst0_size = 10;
        let inst0_coeffs = vec![
            (252, 1),
            (7524, 1),
            (-33561, 1),
            (624024, 1),
            (6958792, 1),
            (-168359184, 1),
            (-7042450329_i64, 4_i64),
            (33379857, 1),
        ];
        let inst0_coeffs: HashSet<_> = inst0_coeffs.into_iter().map(Rational::from).collect();
        assert_eq!(inst_data.inst[0].nonzero.len(), inst0_size);
        assert_eq!(
            inst_data.inst[0]
                .coeffs
                .clone()
                .into_values()
                .collect::<HashSet<_>>(),
            inst0_coeffs
        );

        let inst1_size = 14;
        let inst1_coeffs = vec![
            (-6, 1),
            (504, 1),
            (7164, 1),
            (-67122, 1),
            (-3, 2),
            (-2, 3),
            (-3420, 1),
            (13917584, 1),
            (1248048, 1),
            (1248, 1),
            (-7042450329_i64, 2_i64),
            (-336718368, 1),
            (-3, 8),
            (32846391, 1),
        ];
        let inst1_coeffs: HashSet<_> = inst1_coeffs.into_iter().map(Rational::from).collect();
        assert_eq!(inst_data.inst[1].nonzero.len(), inst1_size);
        assert_eq!(
            inst_data.inst[1]
                .coeffs
                .clone()
                .into_values()
                .collect::<HashSet<_>>(),
            inst1_coeffs
        );

        let expalpha_pos0_size = 10;
        let expalpha_pos0_coeffs = vec![
            1, 60, 3312, -5130, 2772, 343440, 981560, -42569280, 17579592, -240879255,
        ];
        let expalpha_pos0_coeffs: HashSet<_> = expalpha_pos0_coeffs
            .into_iter()
            .map(Rational::from)
            .collect();
        assert_eq!(inst_data.expalpha[0].0.nonzero.len(), expalpha_pos0_size);
        assert_eq!(
            inst_data.expalpha[0]
                .0
                .coeffs
                .clone()
                .into_values()
                .collect::<HashSet<_>>(),
            expalpha_pos0_coeffs
        );

        let expalpha_pos1_size = 11;
        let expalpha_pos1_coeffs = vec![
            1, -60, -540, 8730, 540, -112320, -1813160, 38230920, -2269350, 60, 453347355,
        ];
        let expalpha_pos1_coeffs: HashSet<_> = expalpha_pos1_coeffs
            .into_iter()
            .map(Rational::from)
            .collect();
        assert_eq!(inst_data.expalpha[1].0.nonzero.len(), expalpha_pos1_size);
        assert_eq!(
            inst_data.expalpha[1]
                .0
                .coeffs
                .clone()
                .into_values()
                .collect::<HashSet<_>>(),
            expalpha_pos1_coeffs
        );

        let expalpha_neg0_size = 10;
        let expalpha_neg0_coeffs = vec![
            1, -60, 8730, -3312, -2772, -1813160, 54000, 14031360, -6277608, 453347355,
        ];
        let expalpha_neg0_coeffs: HashSet<_> = expalpha_neg0_coeffs
            .into_iter()
            .map(Rational::from)
            .collect();
        assert_eq!(inst_data.expalpha[0].1.nonzero.len(), expalpha_neg0_size);
        assert_eq!(
            inst_data.expalpha[0]
                .1
                .coeffs
                .clone()
                .into_values()
                .collect::<HashSet<_>>(),
            expalpha_neg0_coeffs
        );

        let expalpha_neg1_size = 11;
        let expalpha_neg1_coeffs = vec![
            1, 60, -5130, 540, -540, 981560, 177120, -28348920, 2496150, -240879255, -60,
        ];
        let expalpha_neg1_coeffs: HashSet<_> = expalpha_neg1_coeffs
            .into_iter()
            .map(Rational::from)
            .collect();
        assert_eq!(inst_data.expalpha[1].1.nonzero.len(), expalpha_neg1_size);
        assert_eq!(
            inst_data.expalpha[1]
                .1
                .coeffs
                .clone()
                .into_values()
                .collect::<HashSet<_>>(),
            expalpha_neg1_coeffs
        );
    }
}
