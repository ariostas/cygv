//! Compute the fundamental period and its derivatives.

use crate::factorial::{factorial_prod, harmonic, RecipMut};
use core::ops::{AddAssign, DivAssign, MulAssign, SubAssign};
use core::slice::Iter;
use nalgebra::{DMatrix, DMatrixView, DVector};
use rug::Assign;
use std::sync::mpsc::{channel, Receiver, Sender};
use std::sync::{Arc, Mutex};

/// Organize curves by the number of negative intersections with the GLSM
/// basis.
pub fn negative_intersections(
    curves_dot_q: DMatrixView<i32>,
) -> (Vec<usize>, Vec<(usize, usize)>, Vec<(usize, usize, usize)>) {
    let mut neg0 = Vec::new();
    let mut neg1 = Vec::new();
    let mut neg2 = Vec::new();

    for (i, col) in curves_dot_q.column_iter().enumerate() {
        let neg_ints: Vec<_> = col
            .iter()
            .enumerate()
            .filter_map(|s| if s.1.is_negative() { Some(s.0) } else { None })
            .collect();
        match neg_ints.len() {
            0 => neg0.push(i),
            1 => neg1.push((i, neg_ints[0])),
            2 => neg2.push((i, neg_ints[0], neg_ints[1])),
            _ => (),
        }
    }

    (neg0, neg1, neg2)
}

/// Computes the c coefficients for curves that have zero negative
/// intersections with the GLSM basis. The results are sent so
/// that the main thread assembles the polynomials.
pub fn compute_c_0neg<T>(
    tasks: Arc<Mutex<Iter<usize>>>,
    tx: Sender<(u32, Option<u32>, Option<u32>, T)>,
    template_var: &T,
    q: DMatrixView<i32>,
    q0: DMatrixView<i32>,
    curves_dot_q: DMatrixView<i32>,
    curves_dot_q0: DMatrixView<i32>,
    beta_pairs: &[(usize, usize)],
) where
    for<'a> T: Clone
        + RecipMut
        + Assign<u64>
        + MulAssign<i32>
        + MulAssign<u64>
        + DivAssign<u64>
        + Assign<&'a T>
        + AddAssign<&'a T>
        + SubAssign<&'a T>
        + MulAssign<&'a T>,
{
    let mut a: Vec<T> = (0..q.ncols()).map(|_| template_var.clone()).collect();
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
        let n: Vec<_> = curves_dot_q0.column(t).iter().map(|c| *c as u32).collect();
        let d: Vec<_> = curves_dot_q.column(t).iter().map(|c| *c as u32).collect();
        factorial_prod(&n, &d, &mut c0fact);
        tx.send((t as u32, None, None, c0fact.clone())).unwrap();
        // Compute A vector
        for (i, aa) in a.iter_mut().enumerate() {
            aa.assign(0);
            for (qq0, cdq0) in q0.column(i).iter().zip(curves_dot_q0.column(t).iter()) {
                harmonic((*cdq0 as u32) + 1, 2, &mut tmp_num0, &mut tmp_num1);
                tmp_num0 *= *qq0;
                *aa += &tmp_num0;
            }
            for (qq, cdq) in q.column(i).iter().zip(curves_dot_q.column(i).iter()) {
                harmonic((*cdq as u32) + 1, 2, &mut tmp_num0, &mut tmp_num1);
                tmp_num0 *= *qq;
                *aa -= &tmp_num0;
            }
            tmp_final.assign(&c0fact);
            tmp_final *= &*aa;
            tx.send((t as u32, Some(i as u32), None, tmp_final.clone()))
                .unwrap();
        }
        // Finally, compute B elements
        for (aa, bb) in beta_pairs.iter() {
            tmp_final.assign(0);
            for ((q0a, q0b), cdq0) in q0
                .column(*aa)
                .iter()
                .zip(q0.column(*bb).iter())
                .zip(curves_dot_q0.column(t).iter())
            {
                harmonic((*cdq0 as u32) + 1, 3, &mut tmp_num0, &mut tmp_num1);
                tmp_num0 *= q0a * q0b;
                tmp_final += &tmp_num0;
            }
            for ((qa, qb), cdq) in q
                .column(*aa)
                .iter()
                .zip(q.column(*bb).iter())
                .zip(curves_dot_q.column(t).iter())
            {
                harmonic((*cdq as u32) + 1, 3, &mut tmp_num0, &mut tmp_num1);
                tmp_num0 *= qa * qb;
                tmp_final -= &tmp_num0;
            }
            tmp_num0.assign(&a[*aa]);
            tmp_num0 *= &a[*bb];
            tmp_final += &tmp_num0;
            tmp_final *= &c0fact;
            tx.send((
                t as u32,
                Some(*aa as u32),
                Some(*bb as u32),
                tmp_final.clone(),
            ))
            .unwrap();
        }
    }
}

/// Computes the c coefficients for curves that have one negative
/// intersection with the GLSM basis. The results are sent so
/// that the main thread assembles the polynomials.
pub fn compute_c_1neg<T>(
    tasks: Arc<Mutex<Iter<(usize, usize)>>>,
    tx: Sender<(u32, Option<u32>, Option<u32>, T)>,
    template_var: &T,
    q: DMatrixView<i32>,
    q0: DMatrixView<i32>,
    curves_dot_q: DMatrixView<i32>,
    curves_dot_q0: DMatrixView<i32>,
    beta_pairs: &[(usize, usize)],
) where
    for<'a> T: Clone
        + RecipMut
        + Assign<u64>
        + MulAssign<i32>
        + MulAssign<u64>
        + DivAssign<u64>
        + Assign<&'a T>
        + AddAssign<&'a T>
        + SubAssign<&'a T>
        + MulAssign<&'a T>,
{
    let mut a: Vec<T> = (0..q.ncols()).map(|_| template_var.clone()).collect();
    let mut tmp_num0 = template_var.clone();
    let mut tmp_num1 = template_var.clone();
    let mut tmp_final = template_var.clone();
    loop {
        let t;
        let k;
        {
            let Some((i, ii)) = tasks.lock().unwrap().next() else {
                break;
            };
            t = *i;
            k = *ii;
        }
        let n: Vec<_> = curves_dot_q0.column(t).iter().map(|c| *c as u32).collect();
        let d: Vec<_> = curves_dot_q.column(t).iter().map(|c| *c as u32).collect();
        factorial_prod(&n, &d, &mut c0fact);
        tx.send((t as u32, None, None, c0fact.clone())).unwrap();
        // Compute A vector
        for (i, aa) in a.iter_mut().enumerate() {
            aa.assign(0);
            for (qq0, cdq0) in q0.column(i).iter().zip(curves_dot_q0.column(t).iter()) {
                harmonic((*cdq0 as u32) + 1, 2, &mut tmp_num0, &mut tmp_num1);
                tmp_num0 *= *qq0;
                *aa += &tmp_num0;
            }
            for (qq, cdq) in q.column(i).iter().zip(curves_dot_q.column(i).iter()) {
                harmonic((*cdq as u32) + 1, 2, &mut tmp_num0, &mut tmp_num1);
                tmp_num0 *= *qq;
                *aa -= &tmp_num0;
            }
            tmp_final.assign(&c0fact);
            tmp_final *= &*aa;
            tx.send((t as u32, Some(i as u32), None, tmp_final.clone()))
                .unwrap();
        }
        // Finally, compute B elements
        for (aa, bb) in beta_pairs.iter() {
            tmp_final.assign(0);
            for ((q0a, q0b), cdq0) in q0
                .column(*aa)
                .iter()
                .zip(q0.column(*bb).iter())
                .zip(curves_dot_q0.column(t).iter())
            {
                harmonic((*cdq0 as u32) + 1, 3, &mut tmp_num0, &mut tmp_num1);
                tmp_num0 *= q0a * q0b;
                tmp_final += &tmp_num0;
            }
            for ((qa, qb), cdq) in q
                .column(*aa)
                .iter()
                .zip(q.column(*bb).iter())
                .zip(curves_dot_q.column(t).iter())
            {
                harmonic((*cdq as u32) + 1, 3, &mut tmp_num0, &mut tmp_num1);
                tmp_num0 *= qa * qb;
                tmp_final -= &tmp_num0;
            }
            tmp_num0.assign(&a[*aa]);
            tmp_num0 *= &a[*bb];
            tmp_final += &tmp_num0;
            tmp_final *= &c0fact;
            tx.send((
                t as u32,
                Some(*aa as u32),
                Some(*bb as u32),
                tmp_final.clone(),
            ))
            .unwrap();
        }
    }
}

// Notes so I can keep track of shapes of matrices. It is important that nalgebra is column-major.
// Q should be h11pd x h11
// Q0 should be codim x h11
// curves_dot_q should be h11 x ncurves
// curves_dot_q0 should be codim x ncurves
