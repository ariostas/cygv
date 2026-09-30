//! Series inversion algorithm.

pub mod error;

use crate::polynomial::{coefficient::PolynomialCoeff, error::PolynomialError};
use crate::{instanton::InstantonData, CYKind, InvariantKind, Polynomial, PolynomialProperties};
use core::cmp::Ordering;
use error::SeriesInversionError;
use nalgebra::DVector;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet, VecDeque};

/// Computes qN for generator curves
fn compute_qn<T>(
    closest_curve: &Polynomial<T>,
    closest_curve_diff: &DVector<i32>,
    expalpha: &[(Polynomial<T>, Polynomial<T>)],
    poly_props: &PolynomialProperties<T>,
) -> Polynomial<T>
where
    T: PolynomialCoeff<T>,
{
    let mut res = closest_curve.clone(&poly_props.zero);
    for (i, diff) in closest_curve_diff.iter().enumerate() {
        let tmp_poly = match diff.cmp(&0) {
            Ordering::Greater => expalpha[i].0.pow(*diff, poly_props).unwrap(),
            Ordering::Less => expalpha[i].1.pow(-*diff, poly_props).unwrap(),
            _ => {
                continue;
            }
        };
        let tmp_poly2 = res.mul(&tmp_poly, poly_props);
        res.clear();
        tmp_poly2.move_into(&mut res);
    }
    res
}

/// Computes qN for a single curve, starting from whichever of the qN of the
/// previous levels is closest to it.
fn compute_qn_from_previous<T>(
    t: usize,
    previous_qn: &VecDeque<HashMap<usize, Polynomial<T>>>,
    previous_qn_ind: &VecDeque<Vec<usize>>,
    expalpha: &[(Polynomial<T>, Polynomial<T>)],
    poly_props: &PolynomialProperties<T>,
) -> Polynomial<T>
where
    T: PolynomialCoeff<T>,
{
    let h11 = poly_props.semigroup.elements.nrows();
    let mut closest_curve = Polynomial::new();
    let mut closest_curve_diff = DVector::zeros(h11);
    let mut tmp_curve_diff = DVector::zeros(h11);
    let mut tmp_num = poly_props.zero.clone();
    tmp_num.assign(1);
    closest_curve.coeffs.insert(t, tmp_num);
    closest_curve.nonzero.push(t);
    closest_curve_diff
        .iter_mut()
        .zip(poly_props.semigroup.elements.column(t).iter())
        .for_each(|(d, s)| *d = *s);
    let mut closest_dist: f32 = closest_curve_diff.iter().map(curve_diff_cost).sum();
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
            let tmp_dist: f32 = tmp_curve_diff.iter().map(curve_diff_cost).sum();
            if tmp_dist < closest_dist {
                let Some(ind) = poly_props.monomial_map.get(&tmp_curve_diff.as_view()) else {
                    continue;
                };
                let mut tmp_poly = Polynomial::new();
                let mut tmp_num = poly_props.zero.clone();
                tmp_num.assign(1);
                tmp_poly.coeffs.insert(*ind, tmp_num);
                tmp_poly.nonzero.push(*ind);
                closest_curve = prev_qns[i].mul(&tmp_poly, poly_props);
                closest_dist = tmp_dist;
                closest_curve_diff.copy_from(&tmp_curve_diff);
            }
        }
    }
    compute_qn(&closest_curve, &closest_curve_diff, expalpha, poly_props)
}

/// How expensive it is to walk one step of a curve-class difference, used to
/// pick the cheapest already-computed qN to start from.
fn curve_diff_cost(d: &i32) -> f32 {
    if *d == 0 {
        0_f32
    } else {
        (*d as f32).abs().log2() + 1_f32
    }
}

/// Find the coefficients of the inverse series, i.e. the GV or GW invariants.
pub fn invert_series<T>(
    inst_data: InstantonData<T>,
    poly_props: &PolynomialProperties<T>,
    invariant_kind: InvariantKind,
    cy_kind: CYKind,
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
        match cy_kind {
            CYKind::Threefold => {
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
                    match invariant_kind {
                        InvariantKind::GV => {
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
                        }
                        InvariantKind::GW => {
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
                }
            }
            CYKind::Nfold => {
                for j in vec_deg {
                    for (k, inst_k) in inst.iter().enumerate() {
                        let Some(gv_ref) = inst_k.coeffs.get(&(j)) else {
                            continue;
                        };
                        tmp_gv.assign(gv_ref);
                        match invariant_kind {
                            InvariantKind::GV => {
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
                            }
                            InvariantKind::GW => {
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
            }
        }
        // Compute qN and Li2(qN) in parallel, then subtract them from the
        // instanton corrections. GW invariants subtract qN itself, so there is
        // no second polynomial to compute or hold on to for them.
        let computed: Vec<(usize, Polynomial<T>, Option<Polynomial<T>>)> = qn_to_compute
            .par_iter()
            .map(|&j| {
                let qn = compute_qn_from_previous(
                    j,
                    &previous_qn,
                    &previous_qn_ind,
                    &expalpha,
                    poly_props,
                );
                let li2qn = match invariant_kind {
                    InvariantKind::GV => Some(qn.li_2(poly_props)?),
                    InvariantKind::GW => None,
                };
                Ok((j, qn, li2qn))
            })
            .collect::<Result<_, PolynomialError>>()?;

        // Every instanton correction is updated independently of the others, so
        // the subtraction is parallelized over them. Each one walks the curves
        // in the same order, which keeps the result independent of scheduling.
        inst.par_iter_mut().enumerate().for_each(|(k, inst_k)| {
            let mut tmp_gv = poly_props.zero.clone();
            for (j, qn, li2qn) in computed.iter() {
                let li2qn = li2qn.as_ref().unwrap_or(qn);
                if cy_kind.is_threefold() {
                    let e = poly_props.semigroup.elements[(k, *j)];
                    if e == 0 {
                        continue;
                    }
                    let mut tmp_poly = li2qn.clone(&poly_props.zero);
                    tmp_gv.assign(&gv_qn_to_compute[j]);
                    tmp_gv *= e;
                    tmp_poly.mul_scalar_assign(&tmp_gv);
                    inst_k.sub_assign(&tmp_poly, &poly_props.zero);
                } else {
                    for kk in h22gv_qn_to_compute[j].iter().filter(|kk| kk.0 == k) {
                        let mut tmp_poly = li2qn.clone(&poly_props.zero);
                        tmp_poly.mul_scalar_assign(&kk.1);
                        inst_k.sub_assign(&tmp_poly, &poly_props.zero);
                    }
                }
            }
        });
        let computed_qn: HashMap<_, _> = computed.into_iter().map(|(j, qn, _)| (j, qn)).collect();
        // Now we update the cache of previous qN
        previous_qn.pop_front();
        previous_qn_ind.pop_front();
        previous_qn.push_back(computed_qn);
        previous_qn_ind.push_back(qn_to_compute);
    }

    Ok(final_gv)
}
