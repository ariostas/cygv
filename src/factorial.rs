//! Factorials and related numbers.

use core::cmp::Ordering;
use core::ops::{AddAssign, DivAssign, MulAssign};
use rug::{Assign, Float, Rational};
use std::collections::HashMap;

/// Sets the value of a number to its reciprocal.
pub trait RecipMut {
    fn recip_mut(&mut self);
}

impl RecipMut for Rational {
    #[inline]
    fn recip_mut(&mut self) {
        Rational::recip_mut(self);
    }
}

impl RecipMut for Float {
    #[inline]
    fn recip_mut(&mut self) {
        Float::recip_mut(self);
    }
}

/// Computes the generalized harmonic number $H_{n,m}$.
///
/// # Examples
///
/// We compute $H_{10,2}=\frac{1968329}{1270080}$.
/// ```
/// use rug::Rational;
/// use cygv::factorial;
///
/// let mut h_10_2 = Rational::new();
/// let mut tmp_var = Rational::new();
///
/// factorial::harmonic(10, 2, &mut h_10_2, &mut tmp_var);
///
/// assert_eq!(h_10_2, Rational::from((1_968_329, 1_270_080)));
/// ```
pub fn harmonic<T>(n: u32, m: u32, res: &mut T, tmp_var: &mut T)
where
    T: RecipMut + Assign<u64> + for<'a> AddAssign<&'a T>,
{
    res.assign(0);
    for i in 1..=n {
        tmp_var.assign((i as u64).pow(m));
        tmp_var.recip_mut();
        *res += &*tmp_var;
    }
}

/// A cache of generalized harmonic numbers $H_{n,m}$ for a fixed $m$.
///
/// The values are computed incrementally, so requesting $H_{n,m}$ only costs
/// one reciprocal and one addition for each $n$ beyond the largest one
/// requested so far. This is much cheaper than [`harmonic`] when the same
/// values are needed repeatedly.
///
/// # Examples
///
/// ```
/// use rug::Rational;
/// use cygv::factorial::HarmonicCache;
///
/// let mut cache = HarmonicCache::new(2, &Rational::new());
///
/// assert_eq!(*cache.get(10), Rational::from((1_968_329, 1_270_080)));
/// assert_eq!(*cache.get(1), Rational::from(1));
/// ```
pub struct HarmonicCache<T> {
    m: u32,
    values: Vec<T>,
    tmp_var: T,
}

impl<T> HarmonicCache<T>
where
    T: Clone + RecipMut + Assign<u64> + for<'a> AddAssign<&'a T>,
{
    /// Create a new cache for the power $m$, using a clonable template number.
    pub fn new(m: u32, template_num: &T) -> Self {
        let mut zero = template_num.clone();
        zero.assign(0u64);
        Self {
            m,
            values: vec![zero],
            tmp_var: template_num.clone(),
        }
    }

    /// Return $H_{n,m}$, extending the cache when needed.
    pub fn get(&mut self, n: u32) -> &T {
        let n = n as usize;
        while self.values.len() <= n {
            let i = self.values.len() as u64;
            self.tmp_var.assign(i.pow(self.m));
            self.tmp_var.recip_mut();
            let mut next = self.values[self.values.len() - 1].clone();
            next += &self.tmp_var;
            self.values.push(next);
        }
        &self.values[n]
    }
}

/// Efficienty computes products and fractions of factorials.
///
/// # Examples
///
/// We perform the following computation.
/// $$\frac{0!2!4!6!8!10!}{1!3!5!7!9!}=3840$$
/// ```
/// use rug::Rational;
/// use cygv::factorial;
///
/// let mut res = Rational::new();
///
/// let num = vec![0,2,4,6,8,10];
/// let den = vec![1,3,5,7,9];
///
/// factorial::factorial_prod(&num, &den, &mut res);
///
/// assert_eq!(res, Rational::from((3_840, 1)));
/// ```
pub fn factorial_prod<T>(n: &[u32], d: &[u32], res: &mut T)
where
    T: Assign<u64> + MulAssign<u64> + DivAssign<u64>,
{
    let mut power_diffs: HashMap<u32, i32> = HashMap::new();
    let mut power = 0i32;

    for c in n.iter() {
        if *c == 0 {
            continue;
        }
        power += 1;
        *power_diffs.entry(*c).or_insert(0) -= 1;
    }
    for c in d.iter() {
        if *c == 0 {
            continue;
        }
        power -= 1;
        *power_diffs.entry(*c).or_insert(0) += 1;
    }

    res.assign(1);

    let Some(max_power) = power_diffs.keys().max() else {
        // res is 1 when given empty slices
        return;
    };

    for i in 1..=*max_power {
        match power.cmp(&0) {
            Ordering::Greater => {
                *res *= (i as u64).pow(power as u32);
            }
            Ordering::Less => {
                *res /= (i as u64).pow(-power as u32);
            }
            _ => {}
        }
        if let Some(c) = power_diffs.get(&i) {
            power += *c;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rug::ops::CompleteRound;

    #[test]
    fn test_harmonic() {
        let mut res_rat = Rational::new();
        let mut tmp_rat = Rational::new();
        let mut true_rat = Rational::new();

        let prec = 512;
        let mut res_float = Float::new(prec);
        let mut tmp_float = Float::new(prec);
        let mut true_float = Float::new(prec);

        let n_m = [(0, 1), (1, 1), (10, 1), (10, 2)];
        let h_n_m = [(0, 1), (1, 1), (7_381, 2_520), (1_968_329, 1_270_080)];

        for (h, (n, m)) in h_n_m.iter().zip(n_m.iter()) {
            harmonic(*n, *m, &mut res_rat, &mut tmp_rat);
            harmonic(*n, *m, &mut res_float, &mut tmp_float);

            true_rat.assign(*h);
            true_float.assign(&true_rat);

            assert_eq!(res_rat, true_rat);
            assert!((&res_float - &true_float).complete(512).abs() < 1e-100);
        }
    }

    #[test]
    fn test_harmonic_cache() {
        let mut res_rat = Rational::new();
        let mut tmp_rat = Rational::new();

        let prec = 512;
        let mut res_float = Float::new(prec);
        let mut tmp_float = Float::new(prec);

        for m in [1, 2] {
            let mut cache_rat = HarmonicCache::new(m, &Rational::new());
            let mut cache_float = HarmonicCache::new(m, &Float::new(prec));

            // Request the values out of order to exercise the incremental fill.
            for n in [10, 0, 3, 10, 25] {
                harmonic(n, m, &mut res_rat, &mut tmp_rat);
                harmonic(n, m, &mut res_float, &mut tmp_float);

                assert_eq!(*cache_rat.get(n), res_rat);
                assert!((cache_float.get(n) - &res_float).complete(512).abs() < 1e-100);
            }
        }
    }

    #[test]
    fn test_factorial_prod() {
        let mut res_rat = Rational::new();
        let mut true_rat = Rational::new();

        let prec = 512;
        let mut res_float = Float::new(prec);
        let mut true_float = Float::new(prec);

        let mut n_d = Vec::new();
        let mut f_n_d = Vec::new();

        n_d.push((vec![], vec![]));
        f_n_d.push((1, 1));

        n_d.push((vec![0, 1, 2, 3, 4, 5], vec![5, 4, 3, 2, 1, 0]));
        f_n_d.push((1, 1));

        n_d.push((vec![0, 2, 4, 6, 8, 10], vec![1, 3, 5, 7, 9]));
        f_n_d.push((3_840, 1));

        n_d.push((vec![90, 25], vec![100]));
        f_n_d.push((34_898_688_000_i64, 141_329));

        for (f, (n, d)) in f_n_d.iter().zip(n_d.iter()) {
            factorial_prod(n, d, &mut res_rat);
            factorial_prod(n, d, &mut res_float);

            true_rat.assign(*f);
            true_float.assign(&true_rat);

            assert_eq!(res_rat, true_rat);
            assert!((&res_float - &true_float).complete(512).abs() < 1e-100);
        }
    }
}
