//! Fundamental period and its derivatives.

pub mod error;

use crate::reference::factorial::{factorial_prod, harmonic};
use crate::reference::polynomial::{coefficient::PolynomialCoeff, Polynomial};
use crate::reference::semigroup::Semigroup;
use crate::{check_strata, CYKind, PolynomialProperties, Stratum};
use error::FundamentalPeriodError;
use nalgebra::{DMatrix, DMatrixView, DVector};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::sync::mpsc::{sync_channel, SyncSender};
use std::thread;

/// The curves that carry period data, grouped by how they get it.
#[derive(Debug, Default, PartialEq)]
struct CurveGroups {
    /// The curves with zero, one and two negative intersections with the toric divisors and
    /// non-negative anticanonical degrees.
    neg: [Vec<usize>; 3],
    /// The curves of negative anticanonical degree whose three negative intersections are
    /// with the divisors of a stratum, with the index of that stratum.
    vex: Vec<(usize, usize)>,
}

/// Groups the curves by the number of negative intersections with the toric divisors.
///
/// A curve with a negative anticanonical degree only contributes through a stratum, when it
/// meets exactly the divisors of one negatively. Every other such curve contributes nothing,
/// and it is an error for one to have fewer than two negative intersections.
fn group_by_neg_int(
    curves_dot_q: DMatrixView<i32>,
    curves_dot_q0: DMatrixView<i32>,
    strata: &[Stratum],
) -> Result<CurveGroups, FundamentalPeriodError> {
    let strata: Vec<[usize; 3]> = strata
        .iter()
        .map(|s| {
            let mut d = s.divisors;
            d.sort_unstable();
            d
        })
        .collect();
    let mut groups = CurveGroups::default();
    let mut neg = Vec::new();
    for (i, (col, col0)) in curves_dot_q
        .column_iter()
        .zip(curves_dot_q0.column_iter())
        .enumerate()
    {
        neg.clear();
        neg.extend((0..col.len()).filter(|&r| col[r] < 0));
        if col0.iter().all(|&m| m >= 0) {
            if let Some(group) = groups.neg.get_mut(neg.len()) {
                group.push(i);
            }
        } else if neg.len() < 2 {
            return Err(FundamentalPeriodError::NegativeAnticanonicalDegree);
        } else if let Some(k) = strata.iter().position(|d| d[..] == neg[..]) {
            groups.vex.push((i, k));
        }
    }
    Ok(groups)
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

/// Computes the coefficients of the stratum terms at a curve of negative anticanonical degree
/// $m$ whose negative intersections $-k_i$ are with the divisors of a stratum: $U_0\, d_{S,t}$
/// for each $t$, with
/// $U_0 = (-1)^{e} \prod_{i \in S} (k_i - 1)! \,/\, \big(\prod_{i \notin S} b_i!\,
/// (\lvert m \rvert - 1)!\big)$, $e = \lvert m \rvert - 1 + \sum_{i \in S} (k_i - 1)$ and $b_i$ the
/// other intersections. The formula is the one of the C program of PR #94.
fn compute_c_vex<T>(
    t: usize,
    stratum: &Stratum,
    template_var: &T,
    curves_dot_q: DMatrixView<i32>,
    curves_dot_q0: DMatrixView<i32>,
) -> Vec<(usize, T)>
where
    T: PolynomialCoeff<T>,
{
    let m = -curves_dot_q0[(0, t)];
    let mut n = Vec::new();
    let mut d = vec![(m - 1) as u32];
    let mut e = m - 1;
    for &c in curves_dot_q.column(t).iter() {
        if c < 0 {
            n.push((-c - 1) as u32);
            e += -c - 1;
        } else {
            d.push(c as u32);
        }
    }
    let mut u0 = template_var.clone();
    factorial_prod(&n, &d, &mut u0);
    if e % 2 != 0 {
        u0 *= -1;
    }
    stratum
        .class
        .iter()
        .enumerate()
        .filter(|(_, &(num, _))| num != 0)
        .map(|(a, &(num, den))| {
            let mut x = u0.clone();
            x *= num;
            x /= den;
            (a, x)
        })
        .collect()
}

/// A struct to keep information about the fundamental period.
/// c0 is the fundamental period, while c1 and c2 are the first and second derivatives of the coefficients.
/// c0_inv is the inverse of the fundamental period.
///
/// In a vex phase, `vex[t]` holds $U_0\, d_{S,t}$ at each curve that meets the divisors of a
/// stratum $S$ negatively (see [`Stratum`]). These terms add to the second derivatives
/// contracted with the intersection numbers, so they enter instanton correction $t$ as
/// $c_0^{-1}$ times `vex[t]`. Without strata the polynomials are empty.
pub struct FundamentalPeriod<T> {
    pub c0: Polynomial<T>,
    pub c1: Vec<Polynomial<T>>,
    pub c2: HashMap<(usize, usize), Polynomial<T>>,
    pub c0_inv: Polynomial<T>,
    pub vex: Vec<Polynomial<T>>,
}

/// Runs one of the `compute_c_*neg` routines over a set of curves in parallel,
/// sending the coefficients of each curve down `tx` as soon as they are ready.
///
/// The scratch space is built once per rayon job rather than once per curve.
fn send_c<T, F>(
    curves: &[usize],
    poly_props: &PolynomialProperties<T>,
    h11: usize,
    tx: &SyncSender<Vec<CCoeff<T>>>,
    compute: F,
) where
    T: PolynomialCoeff<T>,
    F: Fn(usize, &mut CScratch<T>, &mut Vec<CCoeff<T>>) + Sync,
{
    curves.par_iter().for_each_init(
        || CScratch::new(&poly_props.zero, h11),
        |scratch, &t| {
            let mut out = Vec::new();
            compute(t, scratch, &mut out);
            // The receiver only goes away early if filing panicked, in which
            // case that panic is the one worth reporting.
            let _ = tx.send(out);
        },
    );
}

/// Runs `produce`, which sends coefficients down the channel it is given, while
/// a separate thread files them into their polynomials with `file` as they
/// arrive.
///
/// Filing the coefficients as they arrive, rather than collecting them all and
/// filing them afterwards, is what keeps the memory of this stage down: the
/// collected lists would hold every coefficient at once while the polynomials
/// fill up, which puts the peak up by half on the realistic models, and the
/// filing would become a serial tail. The filing thread sits outside the rayon
/// pool, so that waiting on the channel can never take a worker away from the
/// producers, even when the pool only has one.
///
/// The channel is bounded, so that the producers wait for the filing thread
/// whenever it falls behind, e.g. because the machine is busy, rather than
/// piling up coefficients faster than they can be filed. It never waits on the
/// producers, so this cannot deadlock.
fn file_while<T, R>(
    mut file: impl FnMut(CCoeff<T>) + Send,
    produce: impl FnOnce(&SyncSender<Vec<CCoeff<T>>>) -> R,
) -> R
where
    T: PolynomialCoeff<T>,
{
    // A few curves' worth per worker is plenty to smooth out the jitter.
    let (tx, rx) = sync_channel::<Vec<CCoeff<T>>>(8 * rayon::current_num_threads());
    thread::scope(|s| {
        s.spawn(move || {
            for coeffs in rx {
                coeffs.into_iter().for_each(&mut file);
            }
        });
        let res = produce(&tx);
        // Hang up, so that the filing thread stops once the channel is drained.
        drop(tx);
        res
    })
}

/// Files a coefficient that a `compute_c_*neg` routine produced into the
/// polynomial it belongs to. `c0` is only needed while it is being computed.
fn file_coeff<T>(
    (i, a, b, c): CCoeff<T>,
    c0: Option<&mut Polynomial<T>>,
    c1: &mut [Polynomial<T>],
    c2: &mut HashMap<(usize, usize), Polynomial<T>>,
) where
    T: PolynomialCoeff<T>,
{
    match (a, b) {
        (None, _) => {
            c0.expect("a coefficient of c0 arrived after c0 was complete")
                .coeffs
                .insert(i, c);
        }
        (Some(a), None) => {
            c1[a].coeffs.insert(i, c);
        }
        (Some(a), Some(b)) => {
            c2.get_mut(&(a, b)).unwrap().coeffs.insert(i, c);
        }
    }
}

/// Computes the fundamental period $\omega$, its inverse $\omega^{-1}$, and polynomials with first
/// and second derivatives of its coefficients, plus the terms of the strata of a vex phase.
pub fn compute_omega<T>(
    poly_props: &PolynomialProperties<T>,
    sg: &Semigroup,
    q: &DMatrix<i32>,
    nefpart: &[DVector<i32>],
    intnum_idxpairs: &HashSet<(usize, usize)>,
    cy_kind: CYKind,
    strata: &[Stratum],
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
    check_strata(strata, h11pd, h11, cy_kind, nefpart.is_empty())
        .map_err(FundamentalPeriodError::InvalidStrata)?;

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
    let CurveGroups {
        neg: [neg0, neg1, neg2],
        vex: vex_curves,
    } = group_by_neg_int(curves_dot_q.as_view(), curves_dot_q0.as_view(), strata)?;

    let mut c0 = Polynomial::new();
    let mut c1: Vec<Polynomial<T>> = (0..h11).map(|_| Polynomial::new()).collect();
    let mut c2: HashMap<_, _> = beta_pairs
        .iter()
        .map(|&(a, b)| ((a, b), Polynomial::new()))
        .collect();

    // Start by using curves with zero negative intersections.
    file_while(
        |coeff| file_coeff(coeff, Some(&mut c0), &mut c1, &mut c2),
        |tx| {
            send_c(&neg0, poly_props, h11, tx, |t, scratch, out| {
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
            })
        },
    );

    c0.nonzero = c0.coeffs.keys().cloned().collect();
    c0.nonzero.sort_unstable();
    c0.clean_up();

    // Now compute the inverse and the derivatives in parallel.
    let (mut c0_inv, _) = file_while(
        |coeff| file_coeff(coeff, None, &mut c1, &mut c2),
        |tx| {
            rayon::join(
                || c0.recipr(poly_props).unwrap(),
                || {
                    rayon::join(
                        || {
                            send_c(&neg1, poly_props, h11, tx, |t, scratch, out| {
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
                            send_c(&neg2, poly_props, h11, tx, |t, scratch, out| {
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
            )
        },
    );

    c0_inv.clean_up();
    for p in c1.iter_mut() {
        p.nonzero = p.coeffs.keys().cloned().collect();
        p.nonzero.sort_unstable();
        p.clean_up();
    }
    for p in c2.values_mut() {
        p.nonzero = p.coeffs.keys().cloned().collect();
        p.nonzero.sort_unstable();
        p.clean_up();
    }

    let mut vex: Vec<Polynomial<T>> = (0..h11).map(|_| Polynomial::new()).collect();
    let vex_coeffs: Vec<_> = vex_curves
        .par_iter()
        .map(|&(t, k)| {
            let c = compute_c_vex(
                t,
                &strata[k],
                &poly_props.zero,
                curves_dot_q.as_view(),
                curves_dot_q0.as_view(),
            );
            (t, c)
        })
        .collect();
    for (t, c) in vex_coeffs {
        for (a, x) in c {
            vex[a].coeffs.insert(t, x);
        }
    }
    for p in vex.iter_mut() {
        p.nonzero = p.coeffs.keys().cloned().collect();
        p.nonzero.sort_unstable();
    }

    Ok(FundamentalPeriod {
        c0,
        c1,
        c2,
        c0_inv,
        vex,
    })
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
    fn test_group_by_neg_int() {
        // The intersections of six curves with five divisors, and the anticanonical degrees.
        let dq = DMatrix::from_column_slice(
            5,
            6,
            &[
                1, 0, 2, 0, 1, // no negatives
                -1, 2, 0, 0, 1, // one
                -1, -1, 1, 1, 2, // two
                -1, -1, -1, 1, 1, // three, m < 0, on a stratum
                -1, -1, 0, -1, 2, // three, m < 0, on no stratum
                -2, -1, 0, 0, 1, // two, m < 0
            ],
        );
        let m: Vec<i32> = dq.column_iter().map(|c| c.sum()).collect();
        let dq0 = DMatrix::from_row_slice(1, 6, &m);
        let strata = [Stratum {
            divisors: [2, 0, 1],
            class: vec![(1, 1)],
        }];
        let groups = group_by_neg_int(dq.as_view(), dq0.as_view(), &strata).unwrap();
        assert_eq!(
            groups,
            CurveGroups {
                neg: [vec![0], vec![1], vec![2]],
                vex: vec![(3, 0)],
            }
        );
        // Without strata, curves of negative anticanonical degree contribute nothing.
        let groups = group_by_neg_int(dq.as_view(), dq0.as_view(), &[]).unwrap();
        assert!(groups.vex.is_empty());

        // A negative anticanonical degree with a single negative intersection is an error.
        let dq = DMatrix::from_column_slice(3, 1, &[-3, 1, 1]);
        let dq0 = DMatrix::from_row_slice(1, 1, &[-1]);
        assert_eq!(
            group_by_neg_int(dq.as_view(), dq0.as_view(), &[]),
            Err(FundamentalPeriodError::NegativeAnticanonicalDegree)
        );
    }

    #[test]
    fn test_c_vex() {
        // k = (1, 2, 1) on the stratum, b = (0, 3) off it, and m = -1, so that
        // U_0 = (-1)^1 0! 1! 0! / (0! 3! 0!) = -1/6.
        let dq = DMatrix::from_column_slice(5, 1, &[-1, -2, -1, 0, 3]);
        let dq0 = DMatrix::from_row_slice(1, 1, &[-1]);
        let stratum = Stratum {
            divisors: [0, 1, 2],
            class: vec![(1, 2), (0, 1), (-3, 1)],
        };
        let res = compute_c_vex(0, &stratum, &Rational::new(), dq.as_view(), dq0.as_view());
        assert_eq!(
            res,
            vec![(0, Rational::from((-1, 12))), (2, Rational::from((1, 2)))]
        );
    }

    #[test]
    fn test_omega() {
        let generators = DMatrix::from_column_slice(2, 2, &[0, -1, 1, 2]);
        let grading_vector = RowDVector::from_row_slice(&[3, -1]);
        let sg = Semigroup::with_min_elements(generators, grading_vector, 15).unwrap();
        let poly_props = PolynomialProperties::<Rational>::new(&sg);

        let q = DMatrix::from_column_slice(6, 2, &[1, 1, 1, 0, 1, 2, 0, 0, -1, 1, 1, -1]);
        let nefpart = Vec::new();
        let intnum_idxpairs = [(0, 0), (0, 1), (1, 1)].iter().cloned().collect();

        let fp = compute_omega(
            &poly_props,
            &sg,
            &q,
            &nefpart,
            &intnum_idxpairs,
            CYKind::Threefold,
            &[],
        );
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
