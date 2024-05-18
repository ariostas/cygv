//! Affine semigroups that generate SCRP cones.
//!
//! This module contains functions to construct (truncated) affine semigroups
//! that generate scrongly-convex rational polyhedral cones. In other words, the
//! semigroups are of the form $S_\sigma=\sigma\cap\mathbb{Z}^n$ for a
//! strongly-convex rational polyhedal cone $\sigma$ and some
//! $n\in\mathbb{Z}_{>0}$.

pub mod error;

use core::cmp::Ordering;
use error::SemigroupError;
use itertools::Itertools;
use nalgebra::{DMatrix, DVector, RowDVector};
use std::collections::HashSet;

/// A structure for an affine truncated semigroup.
///
/// The affine semigroup needs a grading vector that results in a
/// positive-definite grading, which equivalently means that the semigroup
/// generates a SCRP cone and the grading vector is in the dual cone.
///
/// Functions that use this structure assume that the elements are sorted by
/// degree and that the data is consistent.
#[derive(Clone, Debug, PartialEq)]
pub struct Semigroup {
    pub elements: DMatrix<i32>,
    pub grading_vector: RowDVector<i32>,
    pub degrees: RowDVector<u32>,
    pub max_degree: u32,
}

impl Semigroup {
    /// Constructs a semigroup from given data while only performing essential
    /// checks.
    ///
    /// The matrix of elements must be in column-major format.
    pub fn from_data(
        mut elements: DMatrix<i32>,
        grading_vector: RowDVector<i32>,
    ) -> Result<Self, SemigroupError> {
        let degrees = sort_elements(&mut elements, &grading_vector);
        let degrees = check_final_degrees(degrees)?;
        if elements.column(0).iter().any(|x| *x != 0) {
            return Err(SemigroupError::MissingIdentityError);
        }

        let max_degree = degrees[degrees.len() - 1];

        Ok(Self {
            elements,
            grading_vector,
            degrees,
            max_degree,
        })
    }

    /// Constructs a semigroup given a list of elements containing the
    /// generators, a grading vector, and a maximum degree.
    pub fn with_max_degree(
        elements: DMatrix<i32>,
        grading_vector: RowDVector<i32>,
        max_degree: u32,
    ) -> Result<Self, SemigroupError> {
        // Make sure that the input elements are valid.
        check_degrees(&elements, &grading_vector)?;

        let generators = find_generators(&elements);

        let mut elements_set = HashSet::new();
        for c in elements.column_iter() {
            let tmp_vec = DVector::from_column_slice(c.as_slice());
            elements_set.insert(tmp_vec);
        }

        let mut starting_elements = HashSet::new();
        for c in elements.column_iter() {
            let tmp_vec = DVector::from_column_slice(c.as_slice());
            starting_elements.insert(tmp_vec);
        }
        drop(elements);

        // TODO: This part might be easy to parallelize, so it's worth checking out.
        loop {
            let new_elements = find_new_elements_until_max_deg(
                &generators,
                &starting_elements,
                &elements_set,
                &grading_vector,
                max_degree,
            );
            if new_elements.is_empty() {
                break;
            }
            for c in new_elements.iter() {
                elements_set.insert(c.clone());
            }
            starting_elements = new_elements;
        }
        // Make sure that the zero vector is in the set of elements.
        elements_set.insert(DVector::zeros(grading_vector.len()));

        let mut elements = DMatrix::zeros(grading_vector.len(), elements_set.len());
        elements
            .column_iter_mut()
            .zip(elements_set)
            .for_each(|(mut d, s)| d.copy_from(&s));

        Self::from_data(elements, grading_vector)
    }

    /// Constructs a semigroup by increasing the maximum degree until the minimum number of elements is achieved.
    pub fn with_min_elements(
        elements: DMatrix<i32>,
        grading_vector: RowDVector<i32>,
        min_elements: usize,
    ) -> Result<Self, SemigroupError> {
        // Make sure that the input elements are valid.
        check_degrees(&elements, &grading_vector)?;

        let generators = find_generators(&elements);

        let mut elements_set = HashSet::new();
        for c in elements.column_iter() {
            let tmp_vec = DVector::from_column_slice(c.as_slice());
            elements_set.insert(tmp_vec);
        }
        // Make sure that the zero vector is in the set of elements.
        elements_set.insert(DVector::zeros(grading_vector.len()));

        drop(elements);

        let mut max_degree = 0;
        while elements_set.len() < min_elements {
            max_degree += 1;
            loop {
                let new_elements = find_new_elements_until_max_deg(
                    &generators,
                    &elements_set,
                    &elements_set,
                    &grading_vector,
                    max_degree,
                );
                if new_elements.is_empty() {
                    break;
                }
                for c in new_elements.iter() {
                    elements_set.insert(c.clone());
                }
            }
        }

        let mut elements = DMatrix::zeros(grading_vector.len(), elements_set.len());
        elements
            .column_iter_mut()
            .zip(elements_set)
            .for_each(|(mut d, s)| d.copy_from(&s));

        Self::from_data(elements, grading_vector)
    }
}

/// Sort the elements by degree and make sure that only the identity has degree
/// zero.
fn sort_elements(elements: &mut DMatrix<i32>, grading_vector: &RowDVector<i32>) -> RowDVector<i32> {
    let signed_degrees = grading_vector * &*elements;

    // TODO: figure out if there's a way to sort in place.

    let mut degs_vecs: Vec<(i32, DVector<i32>)> = signed_degrees
        .iter()
        .cloned()
        .zip(elements.column_iter())
        .map(|(d, v)| (d, v.clone_owned()))
        .collect();

    degs_vecs.sort_unstable_by_key(|k| k.0);

    let mut degrees = RowDVector::<i32>::zeros(elements.ncols());
    degrees
        .iter_mut()
        .zip(elements.column_iter_mut())
        .zip(degs_vecs)
        .for_each(|((d, mut v), d_v)| {
            *d = d_v.0;
            v.copy_from(&d_v.1);
        });

    degrees
}

/// Check that the degrees are positive, except for the first one, which must be
/// zero.
fn check_final_degrees(degrees: RowDVector<i32>) -> Result<RowDVector<u32>, SemigroupError> {
    let mut final_degrees = RowDVector::<u32>::zeros(degrees.len());

    let mut degs_iter = degrees.iter();
    let Some(zero_deg) = degs_iter.next() else {
        return Err(SemigroupError::MissingIdentityError);
    };
    match zero_deg.cmp(&0) {
        Ordering::Less => return Err(SemigroupError::NonPositiveDegreeError),
        Ordering::Greater => return Err(SemigroupError::MissingIdentityError),
        _ => (),
    }

    for (d, s) in final_degrees.iter_mut().skip(1).zip(degs_iter) {
        if *s <= 0 {
            return Err(SemigroupError::NonPositiveDegreeError);
        }
        *d = *s as u32
    }

    Ok(final_degrees)
}

fn check_degrees(
    elements: &DMatrix<i32>,
    grading_vector: &RowDVector<i32>,
) -> Result<(), SemigroupError> {
    let signed_degrees = grading_vector * elements;
    for (d, c) in signed_degrees.into_iter().zip(elements.column_iter()) {
        if *d < 0 || (*d == 0 && c.iter().any(|x| *x != 0)) {
            return Err(SemigroupError::NonPositiveDegreeError);
        }
    }
    Ok(())
}

/// Tries to find a smaller subset of elements that generates the semigroup. It
/// does not necessarily return the minimal set of generators since doing so is
/// very difficult.
fn find_generators(elements: &DMatrix<i32>) -> DMatrix<i32> {
    // TODO: Need to check if it is worth to do this in parallel.

    // TODO: Need to check if this is reasonable. For the original code it was
    // 5, but that was probably too high.
    let max_sum_elements = 3;

    let dim = elements.nrows();
    let zero_vec = DVector::<i32>::zeros(dim);
    let mut tmp_vec = zero_vec.clone();

    let mut generators: HashSet<_> = elements.column_iter().collect();
    generators.remove(&zero_vec.as_view());

    let mut to_remove = HashSet::new();

    for n in 2..max_sum_elements {
        for v in generators.iter().combinations_with_replacement(n) {
            tmp_vec.copy_from(&zero_vec);
            for c in v.into_iter() {
                tmp_vec += *c;
            }
            let view = tmp_vec.column(0);
            if generators.contains(&view) {
                to_remove.insert(tmp_vec.clone());
            }
        }
    }

    for c in to_remove.iter() {
        generators.remove(&c.as_view());
    }

    let mut generators_mat = DMatrix::zeros(dim, generators.len());
    generators
        .into_iter()
        .zip(generators_mat.column_iter_mut())
        .for_each(|(s, mut d)| d.copy_from(&s));
    generators_mat
}

/// Find new elements up to a maximum degree using the generators and the starting elements.
fn find_new_elements_until_max_deg(
    generators: &DMatrix<i32>,
    starting_elements: &HashSet<DVector<i32>>,
    elements_set: &HashSet<DVector<i32>>,
    grading_vector: &RowDVector<i32>,
    max_degree: u32,
) -> HashSet<DVector<i32>> {
    let mut new_elements = HashSet::new();
    let mut tmp_vec = DVector::zeros(generators.nrows());
    for c1 in generators.column_iter() {
        for c2 in starting_elements.iter() {
            tmp_vec.copy_from(c2);
            tmp_vec += c1;
            let deg = grading_vector.tr_dot(&tmp_vec) as u32;
            if deg <= max_degree && !elements_set.contains(&tmp_vec) {
                new_elements.insert(tmp_vec.clone());
            }
        }
    }
    new_elements
}

#[cfg(test)]
mod tests {
    use super::*;

    fn example_elements_and_grading_vector() -> (DMatrix<i32>, RowDVector<i32>) {
        #[rustfmt::skip]
        let elements = DMatrix::from_column_slice(2, 6,
            &[
                0, 0,
                1, 0,
                0, 1,
                2, 0,
                1, 1,
                0, 2
            ]
        );
        let grading_vector = RowDVector::from_row_slice(&[1, 1]);
        (elements, grading_vector)
    }

    #[test]
    fn test_semigroup_from_data() {
        let (elements, grading_vector) = example_elements_and_grading_vector();

        let sg_result = Semigroup::from_data(elements.clone(), grading_vector.clone());
        assert!(sg_result.is_ok());
        let sg = sg_result.unwrap();
        assert_eq!(sg.degrees, RowDVector::from_row_slice(&[0, 1, 1, 2, 2, 2]));

        let grading_vector = RowDVector::from_row_slice(&[1, -1]);

        let sg_result = Semigroup::from_data(elements, grading_vector);
        assert!(sg_result.is_err());
        let e = sg_result.err().unwrap();
        assert_eq!(e, SemigroupError::NonPositiveDegreeError);
    }

    #[test]
    fn test_semigroup_with_max_degree() {
        let (elements, grading_vector) = example_elements_and_grading_vector();

        let sg_result = Semigroup::with_max_degree(elements, grading_vector, 3);
        assert!(sg_result.is_ok());
        let sg = sg_result.unwrap();
        assert_eq!(sg.degrees.len(), 10);
        assert_eq!(
            sg.degrees,
            RowDVector::from_row_slice(&[0, 1, 1, 2, 2, 2, 3, 3, 3, 3])
        );
    }

    #[test]
    fn test_semigroup_with_min_elements() {
        let (elements, grading_vector) = example_elements_and_grading_vector();

        let sg_result = Semigroup::with_min_elements(elements, grading_vector, 11);
        assert!(sg_result.is_ok());
        let sg = sg_result.unwrap();
        assert_eq!(sg.degrees.len(), 15);
        assert_eq!(
            sg.degrees,
            RowDVector::from_row_slice(&[0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4])
        );
    }

    #[test]
    fn test_find_generators() {
        let (elements, _) = example_elements_and_grading_vector();

        let generators = find_generators(&elements);
        assert_eq!(generators.ncols(), 2);
    }
}
