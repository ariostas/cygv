//! Affine semigroups that generate SCRP cones
//!
//! This module contains functions to construct (truncated) affine semigroups
//! that generate scrongly-convex rational polyhedral cones. In other words, the
//! semigroups are of the form $S_\sigma=\sigma\cap\mathbb{Z}^n$ for a
//! strongly-convex rational polyhedal cone $\sigma$ and some
//! $n\in\mathbb{Z}_{>0}$.

pub mod error;

use error::SemigroupError;
use nalgebra::{DMatrix, DVector, RowDVector};

/// A structure for an affine truncated semigroup
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
        let degrees = check_degrees(degrees)?;
        if elements.column(0).iter().any(|x| *x != 0) {
            return Err(SemigroupError::MissingIdentityError);
        }

        let max_degree = degrees[(0, degrees.len() - 1)];

        Ok(Self {
            elements,
            grading_vector,
            degrees,
            max_degree,
        })
    }

    // TODO: implement other constructors
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
        .zip(degs_vecs.into_iter())
        .for_each(|((d, mut v), d_v)| {
            *d = d_v.0;
            v.copy_from(&d_v.1);
        });

    degrees
}

/// Check that the degrees are positive, except for the first one, which must be
/// zero.
fn check_degrees(degrees: RowDVector<i32>) -> Result<RowDVector<u32>, SemigroupError> {
    let mut final_degrees = RowDVector::<u32>::zeros(degrees.len());

    let mut degs_iter = degrees.iter();
    let Some(zero_deg) = degs_iter.next() else {
        return Err(SemigroupError::MissingIdentityError);
    };
    if *zero_deg < 0 {
        return Err(SemigroupError::NonPositiveDegreeError);
    } else if *zero_deg > 0 {
        return Err(SemigroupError::MissingIdentityError);
    }

    for (d, s) in final_degrees.iter_mut().skip(1).zip(degs_iter) {
        if *s <= 0 {
            return Err(SemigroupError::NonPositiveDegreeError);
        }
        *d = *s as u32
    }

    Ok(final_degrees)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_semigroup_from_data() {
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
}
