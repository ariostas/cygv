//! Fundamental period and its derivatives.

pub mod error;

use crate::factorial::{factorial_prod, harmonic};
use crate::polynomial::{coefficient::PolynomialCoeff, Polynomial};
use crate::semigroup::Semigroup;
use crate::PolynomialProperties;
use error::FundamentalPeriodError;
use nalgebra::{DMatrix, DMatrixView, DVector};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

/// Group curves by the number of negative intersections with the GLSM
/// basis.
fn group_by_neg_int(curves_dot_q: DMatrixView<i32>) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    let mut neg0 = Vec::new();
    let mut neg1 = Vec::new();
    let mut neg2 = Vec::new();

    for (i, col) in curves_dot_q.column_iter().enumerate() {
        let neg_ints = col.iter().filter(|s| s.is_negative()).count();
        match neg_ints {
            0 => neg0.push(i),
            1 => neg1.push(i),
            2 => neg2.push(i),
            _ => (),
        }
    }

    (neg0, neg1, neg2)
}

/// One coefficient of the fundamental period or of one of its derivatives: the
/// index of the curve it belongs to, the derivative indices (none for `c0`, one
/// for `c1`, two for `c2`), and the value.
type CCoeff<T> = (usize, Option<usize>, Option<usize>, T);

/// Scratch space that the `compute_c_*neg` routines reuse across curves, so
/// that the bignums backing it are allocated once per rayon job rather than
/// once per curve.
struct CScratch<T> {
    /// The A vector, with one entry per element of the GLSM basis.
    a: Vec<T>,
    fact: T,
    tmp0: T,
    tmp1: T,
    res: T,
}

impl<T: PolynomialCoeff<T>> CScratch<T> {
    fn new(template_var: &T, h11: usize) -> Self {
        Self {
            a: (0..h11).map(|_| template_var.clone()).collect(),
            fact: template_var.clone(),
            tmp0: template_var.clone(),
            tmp1: template_var.clone(),
            res: template_var.clone(),
        }
    }
}

/// Computes the c coefficients for a curve that has zero negative
/// intersections with the GLSM basis, pushing them onto `out`.
#[allow(clippy::too_many_arguments)]
fn compute_c_0neg<T>(
    t: usize,
    scratch: &mut CScratch<T>,
    out: &mut Vec<CCoeff<T>>,
    q: DMatrixView<i32>,
    q0: DMatrixView<i32>,
    curves_dot_q: DMatrixView<i32>,
    curves_dot_q0: DMatrixView<i32>,
    beta_pairs: &[(usize, usize)],
) where
    T: PolynomialCoeff<T>,
{
    let CScratch {
        a,
        fact: c0fact,
        tmp0: tmp_num0,
        tmp1: tmp_num1,
        res: tmp_final,
    } = scratch;
    // compute the common c0 factor
    let n: Vec<_> = curves_dot_q0.column(t).iter().map(|&c| c as u32).collect();
    let d: Vec<_> = curves_dot_q.column(t).iter().map(|&c| c as u32).collect();
    factorial_prod(&n, &d, c0fact);
    out.push((t, None, None, c0fact.clone()));
    // Compute A vector
    for (i, aa) in a.iter_mut().enumerate() {
        aa.assign(0);
        for (&qq0, &cdq0) in q0.column(i).iter().zip(curves_dot_q0.column(t).iter()) {
            harmonic(cdq0 as u32, 1, tmp_num0, tmp_num1);
            *tmp_num0 *= qq0;
            *aa += &*tmp_num0;
        }
        for (&qq, &cdq) in q.column(i).iter().zip(curves_dot_q.column(t).iter()) {
            harmonic(cdq as u32, 1, tmp_num0, tmp_num1);
            *tmp_num0 *= qq;
            *aa -= &*tmp_num0;
        }
        tmp_final.assign(&*c0fact);
        *tmp_final *= &*aa;
        out.push((t, Some(i), None, tmp_final.clone()));
    }
    // Finally, compute B elements
    for &(aa, bb) in beta_pairs.iter() {
        tmp_final.assign(0);
        for ((&q0a, &q0b), &cdq0) in q0
            .column(aa)
            .iter()
            .zip(q0.column(bb).iter())
            .zip(curves_dot_q0.column(t).iter())
        {
            harmonic(cdq0 as u32, 2, tmp_num0, tmp_num1);
            *tmp_num0 *= q0a * q0b;
            *tmp_final -= &*tmp_num0;
        }
        for ((&qa, &qb), &cdq) in q
            .column(aa)
            .iter()
            .zip(q.column(bb).iter())
            .zip(curves_dot_q.column(t).iter())
        {
            harmonic(cdq as u32, 2, tmp_num0, tmp_num1);
            *tmp_num0 *= qa * qb;
            *tmp_final += &*tmp_num0;
        }
        tmp_num0.assign(&a[aa]);
        *tmp_num0 *= &a[bb];
        *tmp_final += &*tmp_num0;
        *tmp_final *= &*c0fact;
        out.push((t, Some(aa), Some(bb), tmp_final.clone()));
    }
}

/// Computes the c coefficients for a curve that has one negative
/// intersection with the GLSM basis, pushing them onto `out`.
#[allow(clippy::too_many_arguments)]
fn compute_c_1neg<T>(
    t: usize,
    scratch: &mut CScratch<T>,
    out: &mut Vec<CCoeff<T>>,
    q: DMatrixView<i32>,
    q0: DMatrixView<i32>,
    curves_dot_q: DMatrixView<i32>,
    curves_dot_q0: DMatrixView<i32>,
    beta_pairs: &[(usize, usize)],
) where
    T: PolynomialCoeff<T>,
{
    let CScratch {
        a,
        fact: tmp_fact,
        tmp0: tmp_num0,
        tmp1: tmp_num1,
        res: tmp_final,
    } = scratch;
    let neg_ints: Vec<_> = curves_dot_q
        .column(t)
        .iter()
        .enumerate()
        .filter(|(_, c)| c.is_negative())
        .map(|(i, c)| (i, *c))
        .collect();
    let mut neg_ints_iter = neg_ints.into_iter();
    let (neg_idx, neg_int) = neg_ints_iter
        .next()
        .expect("the curve doesn't have negative intersections");
    assert!(
        neg_ints_iter.next().is_none(),
        "the curve has more than one negative intersection"
    );
    let sn = if neg_int % 2 == 0 { -1 } else { 1 };

    let mut n: Vec<_> = curves_dot_q0.column(t).iter().map(|c| *c as u32).collect();
    n.push((-neg_int - 1) as u32);
    let d: Vec<_> = curves_dot_q
        .column(t)
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != neg_idx)
        .map(|(_, c)| *c as u32)
        .collect();
    factorial_prod(&n, &d, tmp_fact);
    // Compute A vector
    for (i, aa) in a.iter_mut().enumerate() {
        aa.assign(0);
        for (&qq0, &cdq0) in q0.column(i).iter().zip(curves_dot_q0.column(t).iter()) {
            harmonic(cdq0 as u32, 1, tmp_num0, tmp_num1);
            *tmp_num0 *= qq0;
            *aa += &*tmp_num0;
        }
        for (&qq, &cdq) in q.column(i).iter().zip(curves_dot_q.column(t).iter()) {
            harmonic(
                if cdq.is_negative() {
                    (-cdq - 1) as u32
                } else {
                    cdq as u32
                },
                1,
                tmp_num0,
                tmp_num1,
            );
            *tmp_num0 *= qq;
            *aa -= &*tmp_num0;
        }
        tmp_final.assign(&*tmp_fact);
        *tmp_final *= sn;
        *tmp_final *= q[(neg_idx, i)];
        out.push((t, Some(i), None, tmp_final.clone()));
    }
    for &(aa, bb) in beta_pairs.iter() {
        tmp_final.assign(&*tmp_fact);
        tmp_num0.assign(&a[aa]);
        tmp_num1.assign(&a[bb]);
        *tmp_num0 *= q[(neg_idx, bb)];
        *tmp_num1 *= q[(neg_idx, aa)];
        *tmp_num0 += &*tmp_num1;
        *tmp_num0 *= sn;
        *tmp_final *= &*tmp_num0;
        out.push((t, Some(aa), Some(bb), tmp_final.clone()));
    }
}

/// Computes the c coefficients for a curve that has two negative
/// intersections with the GLSM basis, pushing them onto `out`.
fn compute_c_2neg<T>(
    t: usize,
    scratch: &mut CScratch<T>,
    out: &mut Vec<CCoeff<T>>,
    q: DMatrixView<i32>,
    curves_dot_q: DMatrixView<i32>,
    curves_dot_q0: DMatrixView<i32>,
    beta_pairs: &[(usize, usize)],
) where
    T: PolynomialCoeff<T>,
{
    let CScratch {
        fact: tmp_fact,
        res: tmp_final,
        ..
    } = scratch;
    let neg_ints: Vec<_> = curves_dot_q
        .column(t)
        .iter()
        .enumerate()
        .filter(|(_, c)| c.is_negative())
        .map(|(i, c)| (i, *c))
        .collect();
    let mut neg_ints_iter = neg_ints.into_iter();
    let (neg_idx1, neg_int1) = neg_ints_iter
        .next()
        .expect("the curve doesn't have negative intersections");
    let (neg_idx2, neg_int2) = neg_ints_iter
        .next()
        .expect("the curve only has one negative intersection");
    assert!(
        neg_ints_iter.next().is_none(),
        "the curve has more than two negative intersections"
    );
    let sn = if (neg_int1 + neg_int2) % 2 == 0 {
        1
    } else {
        -1
    };

    let mut n: Vec<_> = curves_dot_q0.column(t).iter().map(|c| *c as u32).collect();
    n.push((-neg_int1 - 1) as u32);
    n.push((-neg_int2 - 1) as u32);
    let d: Vec<_> = curves_dot_q
        .column(t)
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != neg_idx1 && *i != neg_idx2)
        .map(|(_, c)| *c as u32)
        .collect();
    factorial_prod(&n, &d, tmp_fact);
    *tmp_fact *= sn;
    for &(aa, bb) in beta_pairs.iter() {
        tmp_final.assign(&*tmp_fact);
        *tmp_final *= q[(neg_idx1, aa)] * q[(neg_idx2, bb)] + q[(neg_idx1, bb)] * q[(neg_idx2, aa)];
        out.push((t, Some(aa), Some(bb), tmp_final.clone()));
    }
}

/// A struct to keep information about the fundamental period.
/// c0 is the fundamental period, while c1 and c2 are the first and second derivatives of the coefficients.
/// c0_inv is the inverse of the fundamental period.
pub struct FundamentalPeriod<T> {
    pub c0: Polynomial<T>,
    pub c1: Vec<Polynomial<T>>,
    pub c2: HashMap<(usize, usize), Polynomial<T>>,
    pub c0_inv: Polynomial<T>,
}

/// Runs one of the `compute_c_*neg` routines over a set of curves in parallel.
///
/// The coefficients come back as one flat list per rayon job. The scratch space
/// and the list are both built once per job rather than once per curve, and the
/// caller files the coefficients into their polynomials afterwards.
fn compute_c<T, F>(
    curves: &[usize],
    poly_props: &PolynomialProperties<T>,
    h11: usize,
    compute: F,
) -> Vec<Vec<CCoeff<T>>>
where
    T: PolynomialCoeff<T>,
    F: Fn(usize, &mut CScratch<T>, &mut Vec<CCoeff<T>>) + Sync,
{
    curves
        .par_iter()
        .fold(
            || (CScratch::new(&poly_props.zero_cutoff, h11), Vec::new()),
            |(mut scratch, mut out), &t| {
                compute(t, &mut scratch, &mut out);
                (scratch, out)
            },
        )
        .map(|(_, out)| out)
        .collect()
}

/// Files the coefficients that `compute_c` produced into the polynomials they
/// belong to.
fn file_away<T>(
    coeffs: Vec<Vec<CCoeff<T>>>,
    c0: &mut Polynomial<T>,
    c1: &mut [Polynomial<T>],
    c2: &mut HashMap<(usize, usize), Polynomial<T>>,
) where
    T: PolynomialCoeff<T>,
{
    for (i, a, b, c) in coeffs.into_iter().flatten() {
        match (a, b) {
            (None, _) => {
                c0.coeffs.insert(i, c);
            }
            (Some(a), None) => {
                c1[a].coeffs.insert(i, c);
            }
            (Some(a), Some(b)) => {
                c2.get_mut(&(a, b)).unwrap().coeffs.insert(i, c);
            }
        }
    }
}

// -> Result<(Polynomial<T>,Vec<Polynomial<T>>,HashMap<(u32,u32),Polynomial<T>>),FundamentalPeriodError>
/// Computes the fundamental period $\omega$, its inverse $\omega^{-1}$, and polynomials with first
/// and second derivatives of its coefficients.
pub fn compute_omega<T>(
    poly_props: &PolynomialProperties<T>,
    sg: &Semigroup,
    q: &DMatrix<i32>,
    nefpart: &[DVector<i32>],
    intnum_idxpairs: &HashSet<(usize, usize)>,
) -> Result<FundamentalPeriod<T>, FundamentalPeriodError>
where
    T: PolynomialCoeff<T>,
{
    let curves = &sg.elements;
    let h11 = q.ncols();
    let h11pd = q.nrows();
    let ambient_dim = (h11pd as i32) - (h11 as i32);
    let cy_codim = if nefpart.is_empty() { 1 } else { nefpart.len() };
    let cy_dim = ambient_dim - (cy_codim as i32);

    // Run some basic checks on the input data
    if cy_dim < 3 {
        return Err(FundamentalPeriodError::CYDimLessThanThree);
    } else if nefpart
        .iter()
        .map(|v| v.iter().max().unwrap_or(&0))
        .any(|c| *c >= (h11pd as i32) || *c < 0)
    {
        return Err(FundamentalPeriodError::InconsistentNefPartition);
    }

    let mut q0 = DMatrix::<i32>::zeros(cy_codim, h11);
    if nefpart.is_empty() {
        for (qq0, q_col) in q0.iter_mut().zip(q.column_iter()) {
            *qq0 = q_col.iter().sum();
        }
    } else {
        for (mut qq0_row, part) in q0.row_iter_mut().zip(nefpart.iter()) {
            for (qq0, q_col) in qq0_row.iter_mut().zip(q.column_iter()) {
                *qq0 = part.iter().map(|c| q_col[*c as usize]).sum();
            }
        }
    }

    let curves_dot_q = q.clone() * curves;
    let curves_dot_q0 = q0.clone() * curves;
    let beta_pairs: Vec<_> = intnum_idxpairs.iter().cloned().collect();
    let (neg0, neg1, neg2) = group_by_neg_int(curves_dot_q.as_view());

    let mut c0 = Polynomial::new();
    let mut c1: Vec<Polynomial<T>> = (0..h11).map(|_| Polynomial::new()).collect();
    let mut c2: HashMap<_, _> = beta_pairs
        .iter()
        .map(|&(a, b)| ((a, b), Polynomial::new()))
        .collect();

    // Start by using curves with zero negative intersections.
    let coeffs_c0 = compute_c(&neg0, poly_props, h11, |t, scratch, out| {
        compute_c_0neg(
            t,
            scratch,
            out,
            q.as_view(),
            q0.as_view(),
            curves_dot_q.as_view(),
            curves_dot_q0.as_view(),
            &beta_pairs,
        )
    });
    file_away(coeffs_c0, &mut c0, &mut c1, &mut c2);

    c0.nonzero = c0.coeffs.keys().cloned().collect();
    c0.nonzero.sort_unstable();
    c0.clean_up(poly_props);

    // Now compute the inverse and the derivatives in parallel.
    let (mut c0_inv, (coeffs_c1, coeffs_c2)) = rayon::join(
        || c0.recipr(poly_props).unwrap(),
        || {
            rayon::join(
                || {
                    compute_c(&neg1, poly_props, h11, |t, scratch, out| {
                        compute_c_1neg(
                            t,
                            scratch,
                            out,
                            q.as_view(),
                            q0.as_view(),
                            curves_dot_q.as_view(),
                            curves_dot_q0.as_view(),
                            &beta_pairs,
                        )
                    })
                },
                || {
                    compute_c(&neg2, poly_props, h11, |t, scratch, out| {
                        compute_c_2neg(
                            t,
                            scratch,
                            out,
                            q.as_view(),
                            curves_dot_q.as_view(),
                            curves_dot_q0.as_view(),
                            &beta_pairs,
                        )
                    })
                },
            )
        },
    );
    file_away(coeffs_c1, &mut c0, &mut c1, &mut c2);
    file_away(coeffs_c2, &mut c0, &mut c1, &mut c2);

    c0_inv.clean_up(poly_props);
    for p in c1.iter_mut() {
        p.nonzero = p.coeffs.keys().cloned().collect();
        p.nonzero.sort_unstable();
        p.clean_up(poly_props);
    }
    for p in c2.values_mut() {
        p.nonzero = p.coeffs.keys().cloned().collect();
        p.nonzero.sort_unstable();
        p.clean_up(poly_props);
    }

    Ok(FundamentalPeriod { c0, c1, c2, c0_inv })
}

// Notes regarding shapes of matrices. Recall that nalgebra is column-major.
// Q should be h11pd x h11
// Q0 should be codim x h11
// curves_dot_q should be h11 x ncurves
// curves_dot_q0 should be codim x ncurves

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::RowDVector;
    use rug::Rational;
    use std::collections::HashSet;

    #[test]
    fn test_omega() {
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

        let c0_size = 4;
        let c0_coeffs = vec![1, 360, 1247400];
        let c0_coeffs: HashSet<_> = c0_coeffs.into_iter().map(Rational::from).collect();
        assert_eq!(fp.c0.nonzero.len(), c0_size);
        assert_eq!(
            fp.c0.coeffs.clone().into_values().collect::<HashSet<_>>(),
            c0_coeffs
        );

        let c1_0_size = 9;
        let c1_0_coeffs = vec![
            60, -6930, 3312, 166320, 2772, 1361360, -36756720, -334639305, 13142520,
        ];
        let c1_0_coeffs: HashSet<_> = c1_0_coeffs.into_iter().map(Rational::from).collect();
        assert_eq!(fp.c1[0].nonzero.len(), c1_0_size);
        assert_eq!(
            fp.c1[0]
                .coeffs
                .clone()
                .into_values()
                .collect::<HashSet<_>>(),
            c1_0_coeffs
        );
        let c1_1_size = 10;
        let c1_1_coeffs = vec![
            -60, 6930, -540, -1361360, -166320, 540, -2598750, 60, 36756720, 334639305,
        ];
        let c1_1_coeffs: HashSet<_> = c1_1_coeffs.into_iter().map(Rational::from).collect();
        assert_eq!(fp.c1[1].nonzero.len(), c1_1_size);
        assert_eq!(
            fp.c1[1]
                .coeffs
                .clone()
                .into_values()
                .collect::<HashSet<_>>(),
            c1_1_coeffs
        );

        let c2_00_size = 9;
        let c2_00_coeffs = vec![
            (1304, 1),
            (-168666, 1),
            (13752, 1),
            (5256, 1),
            (3770784, 1),
            (317945960, 9),
            (-851735880, 1),
            (-18140848109_i64, 2_i64),
            (79322526, 1),
        ];
        let c2_00_coeffs: HashSet<_> = c2_00_coeffs.into_iter().map(Rational::from).collect();
        assert_eq!(fp.c2[&(0, 0)].nonzero.len(), c2_00_size);
        assert_eq!(
            fp.c2[&(0, 0)]
                .coeffs
                .clone()
                .into_values()
                .collect::<HashSet<_>>(),
            c2_00_coeffs
        );
        let c2_01_size = 14;
        let c2_01_coeffs = vec![
            (1, 1),
            (-852, 1),
            (1, 4),
            (-5238, 1),
            (108819, 1),
            (-2403756, 1),
            (1, 9),
            (3258, 1),
            (-68424602, 3),
            (536181858, 1),
            (1, 16),
            (452, 1),
            (-57445875, 2),
            (23478658899_i64, 4_i64),
        ];
        let c2_01_coeffs: HashSet<_> = c2_01_coeffs.into_iter().map(Rational::from).collect();
        assert_eq!(fp.c2[&(0, 1)].nonzero.len(), c2_01_size);
        assert_eq!(
            fp.c2[&(0, 1)]
                .coeffs
                .clone()
                .into_values()
                .collect::<HashSet<_>>(),
            c2_01_coeffs
        );
        let c2_11_size = 14;
        let c2_11_coeffs = vec![
            (2, 1),
            (400, 1),
            (1, 2),
            (1980, 1),
            (-48972, 1),
            (1036728, 1),
            (2, 9),
            (1980, 1),
            (92601652, 9),
            (-220627836, 1),
            (1, 8),
            (400, 1),
            (10308375, 1),
            (-2668905395_i64, 1_i64),
        ];
        let c2_11_coeffs: HashSet<_> = c2_11_coeffs.into_iter().map(Rational::from).collect();
        assert_eq!(fp.c2[&(1, 1)].nonzero.len(), c2_11_size);
        assert_eq!(
            fp.c2[&(1, 1)]
                .coeffs
                .clone()
                .into_values()
                .collect::<HashSet<_>>(),
            c2_11_coeffs
        );
    }
}
