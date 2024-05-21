//! Fundamental period and its derivatives.

pub mod error;

use crate::factorial::{factorial_prod, harmonic};
use crate::polynomial::{coefficient::PolynomialCoeff, Polynomial};
use crate::pool::NumberPool;
use crate::semigroup::Semigroup;
use crate::PolynomialProperties;
use core::slice::Iter;
use error::FundamentalPeriodError;
use nalgebra::{DMatrix, DMatrixView, DVector};
use std::collections::{HashMap, HashSet};
use std::sync::mpsc::{channel, Sender};
use std::sync::{Arc, Mutex};
use std::thread;

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

/// Computes the c coefficients for curves that have zero negative
/// intersections with the GLSM basis. The results are sent so
/// that the main thread assembles the polynomials.
#[allow(clippy::too_many_arguments)]
fn compute_c_0neg<T>(
    tasks: Arc<Mutex<Iter<usize>>>,
    tx: Sender<(u32, Option<u32>, Option<u32>, T)>,
    template_var: &T,
    q: DMatrixView<i32>,
    q0: DMatrixView<i32>,
    curves_dot_q: DMatrixView<i32>,
    curves_dot_q0: DMatrixView<i32>,
    beta_pairs: &[(usize, usize)],
) where
    T: PolynomialCoeff<T>,
{
    let mut a: Vec<_> = (0..q.ncols()).map(|_| template_var.clone()).collect();
    let mut tmp_num0 = template_var.clone();
    let mut tmp_num1 = template_var.clone();
    let mut c0fact = template_var.clone();
    let mut tmp_final = template_var.clone();
    loop {
        let t;
        {
            let Some(i) = tasks.lock().unwrap().next() else {
                break;
            };
            t = *i;
        }
        // compute the common c0 factor
        let n: Vec<_> = curves_dot_q0.column(t).iter().map(|&c| c as u32).collect();
        let d: Vec<_> = curves_dot_q.column(t).iter().map(|&c| c as u32).collect();
        factorial_prod(&n, &d, &mut c0fact);
        tx.send((t as u32, None, None, c0fact.clone())).unwrap();
        // Compute A vector
        for (i, aa) in a.iter_mut().enumerate() {
            aa.assign(0);
            for (&qq0, &cdq0) in q0.column(i).iter().zip(curves_dot_q0.column(t).iter()) {
                harmonic(cdq0 as u32, 1, &mut tmp_num0, &mut tmp_num1);
                tmp_num0 *= qq0;
                *aa += &tmp_num0;
            }
            for (&qq, &cdq) in q.column(i).iter().zip(curves_dot_q.column(t).iter()) {
                harmonic(cdq as u32, 1, &mut tmp_num0, &mut tmp_num1);
                tmp_num0 *= qq;
                *aa -= &tmp_num0;
            }
            tmp_final.assign(&c0fact);
            tmp_final *= &*aa;
            tx.send((t as u32, Some(i as u32), None, tmp_final.clone()))
                .unwrap();
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
                harmonic(cdq0 as u32, 2, &mut tmp_num0, &mut tmp_num1);
                tmp_num0 *= q0a * q0b;
                tmp_final -= &tmp_num0;
            }
            for ((&qa, &qb), &cdq) in q
                .column(aa)
                .iter()
                .zip(q.column(bb).iter())
                .zip(curves_dot_q.column(t).iter())
            {
                harmonic(cdq as u32, 2, &mut tmp_num0, &mut tmp_num1);
                tmp_num0 *= qa * qb;
                tmp_final += &tmp_num0;
            }
            tmp_num0.assign(&a[aa]);
            tmp_num0 *= &a[bb];
            tmp_final += &tmp_num0;
            tmp_final *= &c0fact;
            tx.send((
                t as u32,
                Some(aa as u32),
                Some(bb as u32),
                tmp_final.clone(),
            ))
            .unwrap();
        }
    }
}

/// Computes the c coefficients for curves that have one negative
/// intersection with the GLSM basis. The results are sent so
/// that the main thread assembles the polynomials.
#[allow(clippy::too_many_arguments)]
fn compute_c_1neg<T>(
    tasks: Arc<Mutex<Iter<usize>>>,
    tx: Sender<(u32, Option<u32>, Option<u32>, T)>,
    template_var: &T,
    q: DMatrixView<i32>,
    q0: DMatrixView<i32>,
    curves_dot_q: DMatrixView<i32>,
    curves_dot_q0: DMatrixView<i32>,
    beta_pairs: &[(usize, usize)],
) where
    T: PolynomialCoeff<T>,
{
    let mut a: Vec<_> = (0..q.ncols()).map(|_| template_var.clone()).collect();
    let mut tmp_fact = template_var.clone();
    let mut tmp_num0 = template_var.clone();
    let mut tmp_num1 = template_var.clone();
    let mut tmp_final = template_var.clone();
    loop {
        let t;
        {
            let Some(i) = tasks.lock().unwrap().next() else {
                break;
            };
            t = *i;
        }
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
        factorial_prod(&n, &d, &mut tmp_fact);
        // Compute A vector
        for (i, aa) in a.iter_mut().enumerate() {
            aa.assign(0);
            for (&qq0, &cdq0) in q0.column(i).iter().zip(curves_dot_q0.column(t).iter()) {
                harmonic(cdq0 as u32, 1, &mut tmp_num0, &mut tmp_num1);
                tmp_num0 *= qq0;
                *aa += &tmp_num0;
            }
            for (&qq, &cdq) in q.column(i).iter().zip(curves_dot_q.column(t).iter()) {
                harmonic(
                    if cdq.is_negative() {
                        (-cdq - 1) as u32
                    } else {
                        cdq as u32
                    },
                    1,
                    &mut tmp_num0,
                    &mut tmp_num1,
                );
                tmp_num0 *= qq;
                *aa -= &tmp_num0;
            }
            tmp_final.assign(&tmp_fact);
            tmp_final *= sn;
            tmp_final *= q[(neg_idx, i)];
            tx.send((t as u32, Some(i as u32), None, tmp_final.clone()))
                .unwrap();
        }
        for &(aa, bb) in beta_pairs.iter() {
            tmp_final.assign(&tmp_fact);
            tmp_num0.assign(&a[aa]);
            tmp_num1.assign(&a[bb]);
            tmp_num0 *= q[(neg_idx, bb)];
            tmp_num1 *= q[(neg_idx, aa)];
            tmp_num0 += &tmp_num1;
            tmp_num0 *= sn;
            tmp_final *= &tmp_num0;
            tx.send((
                t as u32,
                Some(aa as u32),
                Some(bb as u32),
                tmp_final.clone(),
            ))
            .unwrap();
        }
    }
}

/// Computes the c coefficients for curves that have two negative
/// intersection with the GLSM basis. The results are sent so
/// that the main thread assembles the polynomials.
fn compute_c_2neg<T>(
    tasks: Arc<Mutex<Iter<usize>>>,
    tx: Sender<(u32, Option<u32>, Option<u32>, T)>,
    template_var: &T,
    q: DMatrixView<i32>,
    curves_dot_q: DMatrixView<i32>,
    curves_dot_q0: DMatrixView<i32>,
    beta_pairs: &[(usize, usize)],
) where
    T: PolynomialCoeff<T>,
{
    let mut tmp_fact = template_var.clone();
    let mut tmp_final = template_var.clone();
    loop {
        let t;
        {
            let Some(i) = tasks.lock().unwrap().next() else {
                break;
            };
            t = *i;
        }
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
        factorial_prod(&n, &d, &mut tmp_fact);
        tmp_fact *= sn;
        for &(aa, bb) in beta_pairs.iter() {
            tmp_final.assign(&tmp_fact);
            tmp_final *=
                q[(neg_idx1, aa)] * q[(neg_idx2, bb)] + q[(neg_idx1, bb)] * q[(neg_idx2, aa)];
            tx.send((
                t as u32,
                Some(aa as u32),
                Some(bb as u32),
                tmp_final.clone(),
            ))
            .unwrap();
        }
    }
}

/// A struct to keep information about the fundamental period.
/// c0 is the fundamental period, while c1 and c2 are the first and second derivatives of the coefficients.
/// c0_inv is the inverse of the fundamental period.
pub struct FundamentalPeriod<T> {
    pub c0: Polynomial<T>,
    pub c1: Vec<Polynomial<T>>,
    pub c2: HashMap<(u32, u32), Polynomial<T>>,
    pub c0_inv: Polynomial<T>,
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
    let mut coeff_pool = NumberPool::new(poly_props.zero_cutoff.clone(), 1000);

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

    let n_threads = thread::available_parallelism()
        .unwrap_or(std::num::NonZeroUsize::new(1).unwrap())
        .get();
    let mut c0 = Polynomial::new();
    let mut c1 = Vec::new();
    for _ in 0..h11 {
        c1.push(Polynomial::new());
    }
    let mut c2 = HashMap::new();
    for &(a, b) in beta_pairs.iter() {
        c2.insert((a as u32, b as u32), Polynomial::new());
    }
    let mut c0_inv = Polynomial::new();

    // Start by using curves with zero negative intersections.
    let tasks_c0 = Arc::new(Mutex::new(neg0.iter()));

    thread::scope(|s| {
        let (tx, rx) = channel();

        for _ in 0..n_threads {
            let tx = tx.clone();
            let tasks = Arc::clone(&tasks_c0);
            s.spawn(|| {
                compute_c_0neg(
                    tasks,
                    tx,
                    &poly_props.zero_cutoff,
                    q.as_view(),
                    q0.as_view(),
                    curves_dot_q.as_view(),
                    curves_dot_q0.as_view(),
                    &beta_pairs,
                )
            });
        }
        drop(tx);

        while let Ok(msg) = rx.recv() {
            match msg {
                (i, None, None, c) => {
                    c0.coeffs.insert(i, c);
                }
                (i, Some(a), None, c) => {
                    c1[a as usize].coeffs.insert(i, c);
                }
                (i, Some(a), Some(b), c) => {
                    c2.get_mut(&(a, b)).unwrap().coeffs.insert(i, c);
                }
                _ => {}
            }
        }
    });
    c0.nonzero = c0.coeffs.keys().cloned().collect();
    c0.nonzero.sort_unstable();
    c0.clean_up(poly_props, &mut coeff_pool);

    // Now compute the inverse and derivatives in parallel
    let tasks_c1 = Arc::new(Mutex::new(neg1.iter()));
    let tasks_c2 = Arc::new(Mutex::new(neg2.iter()));

    thread::scope(|s| {
        // Compute inverse of fundamental period
        s.spawn(|| {
            let tmp_poly = c0.recipr(poly_props, &mut coeff_pool).unwrap();
            tmp_poly.move_into(&mut c0_inv, &mut coeff_pool);
        });

        let (tx, rx) = channel();

        // Compute c1
        for _ in 0..n_threads {
            let tx = tx.clone();
            let tasks = Arc::clone(&tasks_c1);
            s.spawn(|| {
                compute_c_1neg(
                    tasks,
                    tx,
                    &poly_props.zero_cutoff,
                    q.as_view(),
                    q0.as_view(),
                    curves_dot_q.as_view(),
                    curves_dot_q0.as_view(),
                    &beta_pairs,
                )
            });
        }

        // Compute c2
        for _ in 0..n_threads {
            let tx = tx.clone();
            let tasks = Arc::clone(&tasks_c2);
            s.spawn(|| {
                compute_c_2neg(
                    tasks,
                    tx,
                    &poly_props.zero_cutoff,
                    q.as_view(),
                    curves_dot_q.as_view(),
                    curves_dot_q0.as_view(),
                    &beta_pairs,
                )
            });
        }

        drop(tx);

        while let Ok(msg) = rx.recv() {
            match msg {
                (i, Some(a), None, c) => {
                    c1[a as usize].coeffs.insert(i, c);
                }
                (i, Some(a), Some(b), c) => {
                    c2.get_mut(&(a, b)).unwrap().coeffs.insert(i, c);
                }
                _ => {}
            }
        }
    });
    for p in c1.iter_mut() {
        p.nonzero = p.coeffs.keys().cloned().collect();
        p.nonzero.sort_unstable();
        p.clean_up(poly_props, &mut coeff_pool);
    }
    for p in c2.values_mut() {
        p.nonzero = p.coeffs.keys().cloned().collect();
        p.nonzero.sort_unstable();
        p.clean_up(poly_props, &mut coeff_pool);
    }

    for p in c1.iter_mut() {
        p.nonzero = p.coeffs.keys().cloned().collect();
        p.nonzero.sort_unstable();
    }
    for p in c2.values_mut() {
        p.nonzero = p.coeffs.keys().cloned().collect();
        p.nonzero.sort_unstable();
    }

    c0_inv.clean_up(poly_props, &mut coeff_pool);
    for p in c1.iter_mut() {
        p.clean_up(poly_props, &mut coeff_pool);
    }
    for p in c2.values_mut() {
        p.clean_up(poly_props, &mut coeff_pool);
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
        let intnum_idxpairs = vec![(0, 0), (0, 1), (1, 1)].iter().cloned().collect();

        let fp = compute_omega(&poly_props, &sg, &q, &nefpart, &intnum_idxpairs);
        assert!(fp.is_ok());
        let fp = fp.unwrap();

        let c0_size = 4;
        let c0_coeffs = vec![1, 360, 1247400];
        let c0_coeffs: HashSet<_> = c0_coeffs.into_iter().map(|c| Rational::from(c)).collect();
        assert_eq!(fp.c0.nonzero.len(), c0_size);
        assert_eq!(
            fp.c0.coeffs.clone().into_values().collect::<HashSet<_>>(),
            c0_coeffs
        );

        let c1_0_size = 9;
        let c1_0_coeffs = vec![
            60, -6930, 3312, 166320, 2772, 1361360, -36756720, -334639305, 13142520,
        ];
        let c1_0_coeffs: HashSet<_> = c1_0_coeffs.into_iter().map(|c| Rational::from(c)).collect();
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
        let c1_1_coeffs: HashSet<_> = c1_1_coeffs.into_iter().map(|c| Rational::from(c)).collect();
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
        let c2_00_coeffs: HashSet<_> = c2_00_coeffs
            .into_iter()
            .map(|c| Rational::from(c))
            .collect();
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
        let c2_01_coeffs: HashSet<_> = c2_01_coeffs
            .into_iter()
            .map(|c| Rational::from(c))
            .collect();
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
        let c2_11_coeffs: HashSet<_> = c2_11_coeffs
            .into_iter()
            .map(|c| Rational::from(c))
            .collect();
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
