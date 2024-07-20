//! Series inversion algorithm.

pub mod error;

use crate::polynomial::{coefficient::PolynomialCoeff, error::PolynomialError};
use crate::{instanton::InstantonData, NumberPool, Polynomial, PolynomialProperties};
use core::cmp::Ordering;
use core::slice::Iter;
use error::SeriesInversionError;
use nalgebra::DVector;
use std::collections::{HashMap, HashSet, VecDeque};
use std::sync::mpsc::{channel, Sender};
use std::sync::{Arc, Mutex};
use std::thread;

/// Computes qN for generator curves
fn compute_qn<T>(
    closest_curve: &Polynomial<T>,
    closest_curve_diff: &DVector<i32>,
    expalpha: &[(Polynomial<T>, Polynomial<T>)],
    poly_props: &PolynomialProperties<T>,
    np: &mut NumberPool<T>,
) -> Polynomial<T>
where
    T: PolynomialCoeff<T>,
{
    let mut res = closest_curve.clone(np);
    for (i, diff) in closest_curve_diff.iter().enumerate() {
        let tmp_poly = match diff.cmp(&0) {
            Ordering::Greater => expalpha[i].0.pow(*diff, poly_props, np).unwrap(),
            Ordering::Less => expalpha[i].1.pow(-*diff, poly_props, np).unwrap(),
            _ => {
                continue;
            }
        };
        let tmp_poly2 = res.mul(&tmp_poly, poly_props, np);
        res.clear(np);
        tmp_poly2.move_into(&mut res, np);
        tmp_poly.drop(np);
    }
    res
}

/// Computes qN and Li2(qN)
#[allow(clippy::type_complexity)]
fn compute_li2qn_thread<T, const FIND_GV: bool>(
    tasks: Arc<Mutex<Iter<usize>>>,
    tx: Sender<(usize, Polynomial<T>, Result<Polynomial<T>, PolynomialError>)>,
    previous_qn: &VecDeque<HashMap<usize, Polynomial<T>>>,
    previous_qn_ind: &VecDeque<Vec<usize>>,
    expalpha: &[(Polynomial<T>, Polynomial<T>)],
    poly_props: &PolynomialProperties<T>,
    np: &mut NumberPool<T>,
) where
    T: PolynomialCoeff<T>,
{
    let h11 = poly_props.semigroup.elements.nrows();
    let mut closest_curve = Polynomial::new();
    let mut closest_curve_diff = DVector::zeros(h11);
    let mut tmp_curve_diff = DVector::zeros(h11);
    let mut closest_dist: f32;
    let mut tmp_dist: f32;
    loop {
        let t;
        {
            let Some(i) = tasks.lock().unwrap().next() else {
                break;
            };
            t = *i;
        }
        closest_curve.clear(np);
        let mut tmp_num = np.pop();
        tmp_num.assign(1);
        closest_curve.coeffs.insert(t, tmp_num);
        closest_curve.nonzero.push(t);
        closest_curve_diff
            .iter_mut()
            .zip(poly_props.semigroup.elements.column(t).iter())
            .for_each(|(d, s)| *d = *s);
        closest_dist = closest_curve_diff
            .iter()
            .map(|d| {
                if *d == 0 {
                    0_f32
                } else {
                    (*d as f32).abs().log2() + 1_f32
                }
            })
            .sum();
        // Now check to see if there is a better starting curve
        for (prev_inds, prev_qns) in previous_qn_ind.iter().zip(previous_qn.iter()) {
            for i in prev_inds {
                poly_props
                    .semigroup
                    .elements
                    .column(t)
                    .iter()
                    .zip(poly_props.semigroup.elements.column(*i).iter())
                    .zip(tmp_curve_diff.iter_mut())
                    .for_each(|((s1, s2), d)| *d = s1 - s2);
                tmp_dist = tmp_curve_diff
                    .iter()
                    .map(|d| {
                        if *d == 0 {
                            0_f32
                        } else {
                            (*d as f32).abs().log2() + 1_f32
                        }
                    })
                    .sum();
                if tmp_dist < closest_dist {
                    let Some(ind) = poly_props.monomial_map.get(&tmp_curve_diff.as_view()) else {
                        continue;
                    };
                    let mut tmp_poly = Polynomial::new();
                    let mut tmp_num = np.pop();
                    tmp_num.assign(1);
                    tmp_poly.coeffs.insert(*ind, tmp_num);
                    tmp_poly.nonzero.push(*ind);
                    closest_curve.clear(np);
                    closest_curve = prev_qns[i].mul(&tmp_poly, poly_props, np);
                    tmp_poly.drop(np);
                    closest_dist = tmp_dist;
                    closest_curve_diff
                        .iter_mut()
                        .zip(tmp_curve_diff.iter())
                        .for_each(|(d, s)| *d = *s);
                }
            }
        }
        // Now we compute qN and Li2(qN)
        let tmp_qn = compute_qn(
            &closest_curve,
            &closest_curve_diff,
            expalpha,
            poly_props,
            np,
        );
        let tmp_li2qn = if FIND_GV {
            tmp_qn.li_2(poly_props, np)
        } else {
            Ok(tmp_qn.clone(np))
        };
        tx.send((t, tmp_qn, tmp_li2qn)).unwrap();
    }
}

/// Find the coefficients of the inverse series, i.e. the GV or GW invariants.
pub fn invert_series<T, const FIND_GV: bool, const IS_THREEFOLD: bool>(
    inst_data: InstantonData<T>,
    poly_props: &PolynomialProperties<T>,
    all_pools: &mut (NumberPool<T>, Vec<NumberPool<T>>),
) -> Result<HashMap<(usize, usize), T>, SeriesInversionError>
where
    T: PolynomialCoeff<T>,
{
    let mut final_gv = HashMap::new();

    let h11 = poly_props.semigroup.elements.nrows();
    let n_previous_levels = if h11 < 4 {
        2
    } else if h11 < 10 {
        5
    } else {
        10
    };
    let mut tmp_gv = poly_props.zero_cutoff.clone();
    let mut tmp_gv_rounded = poly_props.zero_cutoff.clone();
    let mut previous_qn: VecDeque<_> = (1..=n_previous_levels).map(|_| HashMap::new()).collect();
    let mut previous_qn_ind: VecDeque<_> = (1..=n_previous_levels).map(|_| Vec::new()).collect();

    let (main_pool, pools) = all_pools;

    let InstantonData { mut inst, expalpha } = inst_data;

    let all_degs: HashSet<_> = poly_props
        .semigroup
        .degrees
        .iter()
        .cloned()
        .filter(|c| *c != 0)
        .collect();
    let mut distinct_degs: Vec<_> = all_degs.into_iter().collect();
    distinct_degs.sort_unstable();

    for d in distinct_degs.iter() {
        let mut vec_deg = Vec::new();
        let mut qn_to_compute = Vec::new();
        let mut gv_qn_to_compute = HashMap::new();
        let mut h22gv_qn_to_compute = HashMap::new();
        // First find the points of interest
        for (i, dd) in poly_props.semigroup.degrees.iter().enumerate() {
            match dd.cmp(d) {
                Ordering::Equal => {
                    vec_deg.push(i);
                }
                Ordering::Greater => {
                    break;
                }
                _ => {}
            }
        }
        if IS_THREEFOLD {
            for j in vec_deg {
                let kk = poly_props
                    .semigroup
                    .elements
                    .column(j)
                    .iter()
                    .cloned()
                    .enumerate()
                    .find(|(_, c)| *c != 0)
                    .unwrap();
                let Some(gv_ref) = inst[kk.0].coeffs.get(&(j)) else {
                    continue;
                };
                tmp_gv.assign(gv_ref);
                tmp_gv /= kk.1;
                if FIND_GV {
                    tmp_gv_rounded.assign(&tmp_gv);
                    tmp_gv_rounded.round_mut();
                    tmp_gv -= &tmp_gv_rounded;
                    tmp_gv.abs_mut();
                    if tmp_gv > 1e-3 {
                        return Err(SeriesInversionError::NonIntegerGVError);
                    }
                    tmp_gv.assign(&tmp_gv_rounded);
                    tmp_gv.abs_mut();
                    if tmp_gv < 0.5 {
                        continue;
                    }
                    final_gv.insert((j, 0), tmp_gv_rounded.clone());
                    qn_to_compute.push(j);
                    gv_qn_to_compute.insert(j, tmp_gv_rounded.clone());
                } else {
                    tmp_gv_rounded.assign(&tmp_gv);
                    tmp_gv_rounded.abs_mut();
                    if tmp_gv_rounded <= poly_props.zero_cutoff {
                        continue;
                    }
                    final_gv.insert((j, 0), tmp_gv.clone());
                    qn_to_compute.push(j);
                    gv_qn_to_compute.insert(j, tmp_gv.clone());
                }
            }
        } else {
            for j in vec_deg {
                for (k, inst_k) in inst.iter().enumerate() {
                    let Some(gv_ref) = inst_k.coeffs.get(&(j)) else {
                        continue;
                    };
                    tmp_gv.assign(gv_ref);
                    if FIND_GV {
                        tmp_gv_rounded.assign(&tmp_gv);
                        tmp_gv_rounded.round_mut();
                        tmp_gv -= &tmp_gv_rounded;
                        tmp_gv.abs_mut();
                        if tmp_gv > 1e-3 {
                            return Err(SeriesInversionError::NonIntegerGVError);
                        }
                        tmp_gv.assign(&tmp_gv_rounded);
                        tmp_gv.abs_mut();
                        if tmp_gv < 0.5 {
                            continue;
                        }
                        final_gv.insert((j, k), tmp_gv_rounded.clone());
                        let h22list = h22gv_qn_to_compute.entry(j).or_insert_with(|| {
                            qn_to_compute.push(j);
                            Vec::new()
                        });
                        h22list.push((k, tmp_gv_rounded.clone()));
                    } else {
                        tmp_gv_rounded.assign(&tmp_gv);
                        tmp_gv_rounded.abs_mut();
                        if tmp_gv_rounded <= poly_props.zero_cutoff {
                            continue;
                        }
                        final_gv.insert((j, k), tmp_gv.clone());
                        let h22list = h22gv_qn_to_compute.entry(j).or_insert_with(|| {
                            qn_to_compute.push(j);
                            Vec::new()
                        });
                        h22list.push((k, tmp_gv.clone()));
                    }
                }
            }
        }
        let mut computed_qn = HashMap::new();
        // compute qN and Li2(qN) in parallel and subtract from the instanton corrections
        let tasks_iter = Arc::new(Mutex::new(qn_to_compute.iter()));
        let mut error = None;
        thread::scope(|s| {
            let (tx, rx) = channel();
            for np in pools.iter_mut() {
                let tx = tx.clone();
                let tasks = Arc::clone(&tasks_iter);
                s.spawn(|| {
                    compute_li2qn_thread::<T, FIND_GV>(
                        tasks,
                        tx,
                        &previous_qn,
                        &previous_qn_ind,
                        &expalpha,
                        poly_props,
                        np,
                    );
                });
            }
            drop(tx);
            while let Ok((j, qn, li2qn_r)) = rx.recv() {
                let Ok(li2qn) = li2qn_r else {
                    error = li2qn_r.err();
                    break;
                };
                computed_qn.insert(j, qn.clone(main_pool));
                if IS_THREEFOLD {
                    for (k, inst_k) in inst.iter_mut().enumerate() {
                        if poly_props.semigroup.elements[(k, j)] == 0 {
                            continue;
                        }
                        let mut tmp_poly = li2qn.clone(main_pool);
                        tmp_gv.assign(&gv_qn_to_compute[&j]);
                        tmp_gv *= poly_props.semigroup.elements[(k, j)];
                        tmp_poly.mul_scalar_assign(&tmp_gv);
                        inst_k.sub_assign(&tmp_poly, main_pool);
                        tmp_poly.drop(main_pool);
                    }
                } else {
                    for kk in h22gv_qn_to_compute[&j].iter() {
                        let mut tmp_poly = li2qn.clone(main_pool);
                        tmp_poly.mul_scalar_assign(&kk.1);
                        inst[kk.0].sub_assign(&tmp_poly, main_pool);
                        tmp_poly.drop(main_pool);
                    }
                }
            }
        });
        // Now we update the cache of previous qN
        previous_qn.pop_front();
        previous_qn_ind.pop_front();
        previous_qn.push_back(computed_qn);
        previous_qn_ind.push_back(qn_to_compute);
    }

    Ok(final_gv)
}
