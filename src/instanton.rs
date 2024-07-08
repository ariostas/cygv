//! Instanton corrections.

use crate::fundamental_period::FundamentalPeriod;
use crate::polynomial::{
    coefficient::PolynomialCoeff, error::PolynomialError, properties::PolynomialProperties,
    Polynomial,
};
use std::collections::{HashMap, HashSet};
use std::sync::mpsc::{channel, Sender};
use std::sync::{Arc, Mutex};
use std::thread;

pub struct InstantonData<T> {
    pub inst: Vec<Polynomial<T>>,
    pub expalpha: Vec<(Polynomial<T>, Polynomial<T>)>,
}

fn compute_alpha_thread<T>(
    tasks: Arc<Mutex<core::slice::Iter<usize>>>,
    tx: Sender<(usize, Polynomial<T>)>,
    fp: &FundamentalPeriod<T>,
    poly_props: &PolynomialProperties<T>,
) where
    T: PolynomialCoeff<T>,
{
    loop {
        let t;
        {
            let Some(i) = tasks.lock().unwrap().next() else {
                break;
            };
            t = *i;
        }
        let mut a = fp.c0_inv.mul(&fp.c1[t], poly_props);
        a.clean_up(poly_props);
        tx.send((t, a)).unwrap();
    }
}

fn compute_beta_thread<T>(
    tasks: Arc<Mutex<std::collections::hash_set::Iter<(usize, usize)>>>,
    tx: Sender<((usize, usize), Polynomial<T>)>,
    fp: &FundamentalPeriod<T>,
    poly_props: &PolynomialProperties<T>,
) where
    T: PolynomialCoeff<T>,
{
    loop {
        let (t0, t1);
        {
            let Some(i) = tasks.lock().unwrap().next() else {
                break;
            };
            t0 = i.0;
            t1 = i.1;
        }
        let mut a = fp.c0_inv.mul(&fp.c2[&(t0, t1)], poly_props);
        a.clean_up(poly_props);
        tx.send(((t0, t1), a)).unwrap();
    }
}

fn compute_f_thread<T>(
    tasks: Arc<Mutex<std::collections::hash_set::Iter<(usize, usize)>>>,
    tx: Sender<((usize, usize), Polynomial<T>)>,
    alpha: &[Polynomial<T>],
    beta: &HashMap<(usize, usize), Polynomial<T>>,
    poly_props: &PolynomialProperties<T>,
) where
    T: PolynomialCoeff<T>,
{
    loop {
        let (t0, t1);
        {
            let Some(i) = tasks.lock().unwrap().next() else {
                break;
            };
            t0 = i.0;
            t1 = i.1;
        }
        let mut p = alpha[t0].mul(&alpha[t1], poly_props);
        p.sub_assign(&beta[&(t0, t1)], poly_props);
        p.mul_scalar_assign(-1);
        p.clean_up(poly_props);
        tx.send(((t0, t1), p)).unwrap();
    }
}

fn compute_inst_thread<T>(
    tasks: Arc<Mutex<core::slice::Iter<usize>>>,
    tx: Sender<(usize, Polynomial<T>)>,
    f_poly: &HashMap<(usize, usize), Polynomial<T>>,
    poly_props: &PolynomialProperties<T>,
    intnum_dict: &HashMap<(usize, usize, usize), i32>,
    is_threefold: bool,
) where
    T: PolynomialCoeff<T>,
{
    let h11 = poly_props.semigroup.elements.nrows();
    let mut intnum_ind = [0_usize; 3];
    let mut tmp_num = poly_props.new_variable();
    loop {
        let t;
        {
            let Some(i) = tasks.lock().unwrap().next() else {
                break;
            };
            t = *i;
        }
        let mut p = Polynomial::new();
        for a in 0..h11 {
            for b in a..h11 {
                intnum_ind[0] = t;
                intnum_ind[1] = a;
                intnum_ind[2] = b;
                if is_threefold {
                    intnum_ind.sort_unstable();
                }
                let Some(x) = intnum_dict.get(&(intnum_ind[0], intnum_ind[1], intnum_ind[2]))
                else {
                    continue;
                };
                let mut tmp_poly = f_poly[&(a, b)].clone();
                if a != b {
                    tmp_poly.mul_scalar_assign(*x);
                } else {
                    tmp_num.assign(*x);
                    tmp_num /= 2;
                    tmp_poly.mul_scalar_assign(&tmp_num);
                }
                p.add_assign(&tmp_poly, poly_props);
            }
        }
        p.clean_up(poly_props);
        tx.send((t, p)).unwrap();
    }
}

#[allow(clippy::type_complexity)]
fn compute_expalpha_thread<T>(
    tasks: Arc<Mutex<core::slice::Iter<usize>>>,
    tx: Sender<(
        usize,
        Result<(Polynomial<T>, Polynomial<T>), PolynomialError>,
    )>,
    alpha: &[Polynomial<T>],
    poly_props: &PolynomialProperties<T>,
) where
    T: PolynomialCoeff<T>,
{
    loop {
        let t;
        {
            let Some(i) = tasks.lock().unwrap().next() else {
                break;
            };
            t = *i;
        }
        let p = alpha[t].exp_pos_neg(poly_props);
        let p = match p {
            Ok(mut pp) => {
                pp.0.clean_up(poly_props);
                pp.1.clean_up(poly_props);
                Ok(pp)
            }
            Err(e) => Err(e),
        };
        tx.send((t, p)).unwrap();
    }
}

/// Compute instanton corrections, as well as other objects needed for the series inversion.
pub fn compute_instanton_data<T>(
    fp: FundamentalPeriod<T>,
    poly_props: &PolynomialProperties<T>,
    intnum_idxpairs: &HashSet<(usize, usize)>,
    n_indices: usize,
    intnum_dict: &HashMap<(usize, usize, usize), i32>,
    is_threefold: bool,
) -> Result<InstantonData<T>, PolynomialError>
where
    T: PolynomialCoeff<T>,
{
    let h11 = poly_props.semigroup.elements.nrows();
    let n_threads = thread::available_parallelism()
        .unwrap_or(std::num::NonZeroUsize::new(1).unwrap())
        .get();

    // Compute alpha polynomials
    let mut alpha: Vec<_> = (0..h11).map(|_| Polynomial::<T>::new()).collect();
    let tasks_alpha: Vec<_> = (0..h11).collect();
    let tasks_alpha_iter = Arc::new(Mutex::new(tasks_alpha.iter()));
    thread::scope(|s| {
        let (tx, rx) = channel();
        for _ in 0..n_threads {
            let tx = tx.clone();
            let tasks = Arc::clone(&tasks_alpha_iter);
            s.spawn(|| {
                compute_alpha_thread(tasks, tx, &fp, poly_props);
            });
        }
        drop(tx);
        while let Ok((t, p)) = rx.recv() {
            alpha[t] = p;
        }
    });

    // Compute beta polynomials
    let mut beta = HashMap::new();
    let tasks_beta_iter = Arc::new(Mutex::new(intnum_idxpairs.iter()));
    thread::scope(|s| {
        let (tx, rx) = channel();
        for _ in 0..n_threads {
            let tx = tx.clone();
            let tasks = Arc::clone(&tasks_beta_iter);
            s.spawn(|| {
                compute_beta_thread(tasks, tx, &fp, poly_props);
            });
        }
        drop(tx);
        while let Ok((t, p)) = rx.recv() {
            beta.insert(t, p);
        }
    });

    // Compute F polynomials
    let mut f_poly = HashMap::new();
    let tasks_f_iter = Arc::new(Mutex::new(intnum_idxpairs.iter()));
    thread::scope(|s| {
        let (tx, rx) = channel();
        for _ in 0..n_threads {
            let tx = tx.clone();
            let tasks = Arc::clone(&tasks_f_iter);
            s.spawn(|| {
                compute_f_thread(tasks, tx, &alpha, &beta, poly_props);
            });
        }
        drop(tx);
        while let Ok((t, p)) = rx.recv() {
            f_poly.insert(t, p);
        }
    });

    // Compute instanton corrections
    let mut inst: Vec<_> = (0..n_indices).map(|_| Polynomial::<T>::new()).collect();
    let tasks_inst: Vec<_> = (0..n_indices).collect();
    let tasks_inst_iter = Arc::new(Mutex::new(tasks_inst.iter()));
    thread::scope(|s| {
        let (tx, rx) = channel();
        for _ in 0..n_threads {
            let tx = tx.clone();
            let tasks = Arc::clone(&tasks_inst_iter);
            s.spawn(|| {
                compute_inst_thread(tasks, tx, &f_poly, poly_props, intnum_dict, is_threefold);
            });
        }
        drop(tx);
        while let Ok((t, p)) = rx.recv() {
            inst[t] = p;
        }
    });

    // Compute expalpha polynomials
    let mut expalpha: Vec<_> = (0..h11)
        .map(|_| (Polynomial::<T>::new(), Polynomial::<T>::new()))
        .collect();
    let tasks_expalpha: Vec<_> = (0..h11).collect();
    let tasks_expalpha_iter = Arc::new(Mutex::new(tasks_expalpha.iter()));
    let mut error = None;
    thread::scope(|s| {
        let (tx, rx) = channel();
        for _ in 0..n_threads {
            let tx = tx.clone();
            let tasks = Arc::clone(&tasks_expalpha_iter);
            s.spawn(|| {
                compute_expalpha_thread(tasks, tx, &alpha, poly_props);
            });
        }
        drop(tx);
        while let Ok((t, p)) = rx.recv() {
            match p {
                Ok(pp) => expalpha[t] = pp,
                Err(e) => {
                    error = Some(e);
                    break;
                }
            }
        }
    });

    if let Some(e) = error {
        return Err(e);
    }

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
        let result = process_int_nums(intnums.clone(), true);
        assert!(result.is_ok());
        let (intnum_dict, intnum_idxpairs, n_indices) = result.unwrap();

        let inst_data = compute_instanton_data(
            fp,
            &poly_props,
            &intnum_idxpairs,
            n_indices,
            &intnum_dict,
            true,
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
