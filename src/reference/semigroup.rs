//! Affine semigroups that generate SCRP cones.
//!
//! This module contains functions to construct (truncated) affine semigroups
//! that generate strongly-convex rational polyhedral cones. In other words, the
//! semigroups are of the form $S_\sigma=\sigma\cap\mathbb{Z}^n$ for a
//! strongly-convex rational polyhedral cone $\sigma$ and some
//! $n\in\mathbb{Z}_{>0}$.

pub mod error;

use core::cmp::Ordering;
use error::SemigroupError;
use nalgebra::{DMatrix, DVector, RowDVector};
use rayon::prelude::*;
use rustc_hash::FxHashSet as HashSet;

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
    /// The matrix of elements must be in column-major format, with each element once.
    pub fn from_data(
        mut elements: DMatrix<i32>,
        grading_vector: RowDVector<i32>,
    ) -> Result<Self, SemigroupError> {
        let degrees = sort_elements(&mut elements, &grading_vector);
        // Sorted, an element given twice is next to itself.
        let cols: Vec<&[i32]> = crate::cone::columns(&elements).collect();
        if cols.windows(2).any(|w| w[0] == w[1]) {
            return Err(SemigroupError::DuplicateElementError);
        }
        drop(cols);
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

        if let Some(elements) = enumerate(&elements, &grading_vector, max_degree, usize::MAX) {
            return Self::from_data(elements?, grading_vector);
        }

        // Remove elements that exceed the maximum degree from the start.
        let elements = trim_by_max_deg(&elements, &grading_vector, max_degree);

        let generators = find_generators(&elements);

        let mut elements_set = matrix_columns_to_set(&elements);
        let mut starting_elements: Vec<_> = elements_set.iter().cloned().collect();
        drop(elements);

        loop {
            crate::interrupt::check(SemigroupError::Interrupted)?;
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
            elements_set.extend(new_elements.iter().cloned());
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

        // Up to a degree that is doubled until it holds enough elements, or the keys no longer
        // hold them.
        let top_generator = elements
            .column_iter()
            .map(|g| degree(&grading_vector, g))
            .max()
            .unwrap_or(0);
        let mut bound = u32::try_from(top_generator.max(1)).unwrap_or(u32::MAX);
        loop {
            match enumerate(&elements, &grading_vector, bound, min_elements) {
                Some(Ok(found)) if found.ncols() >= min_elements || top_generator == 0 => {
                    return Self::from_data(found, grading_vector);
                }
                Some(Err(e)) => return Err(e),
                Some(Ok(_)) if bound < u32::MAX / 2 => bound *= 2,
                _ => break,
            }
        }

        let generators = find_generators(&elements);

        let mut elements_set = HashSet::default();
        elements_set.insert(DVector::zeros(grading_vector.len()));

        drop(elements);

        let top_generator = generators
            .column_iter()
            .map(|g| degree(&grading_vector, g))
            .max();
        // Without generators the semigroup is the identity alone.
        let mut max_degree: u32 = 0;
        while elements_set.len() < min_elements && generators.ncols() > 0 {
            max_degree += 1;
            // Raising the degree bound lets through the sums of the elements within the
            // largest degree of a generator below it, so each round starts from those,
            // and afterwards only the frontier of newly found elements is expanded.
            let low = max_degree as i64 - top_generator.unwrap_or(0);
            let mut starting_elements: Vec<DVector<i32>> = elements_set
                .iter()
                .filter(|&c| degree(&grading_vector, c) >= low)
                .cloned()
                .collect();
            loop {
                crate::interrupt::check(SemigroupError::Interrupted)?;
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
                elements_set.extend(new_elements.iter().cloned());
                starting_elements = new_elements;
            }
        }

        let mut elements = DMatrix::zeros(grading_vector.len(), elements_set.len());
        elements
            .column_iter_mut()
            .zip(elements_set)
            .for_each(|(mut d, s)| d.copy_from(&s));

        Self::from_data(elements, grading_vector)
    }

    /// Constructs the part of the semigroup needed to compute invariants at
    /// the requested target points.
    pub fn with_target_points(
        elements: DMatrix<i32>,
        grading_vector: RowDVector<i32>,
        target_points: DMatrix<i32>,
    ) -> Result<Self, SemigroupError> {
        if target_points.ncols() == 0 {
            return Self::from_data(elements, grading_vector);
        }
        if elements.nrows() != target_points.nrows()
            || grading_vector.len() != target_points.nrows()
        {
            return Err(SemigroupError::DimensionMismatchError);
        }

        check_degrees(&elements, &grading_vector)?;
        check_degrees(&target_points, &grading_vector)?;

        let target_max_degree = target_points
            .column_iter()
            .map(|t| degree(&grading_vector, t))
            .max()
            .unwrap_or(0);
        let target_max_degree =
            u32::try_from(target_max_degree).map_err(|_| SemigroupError::DegreeTooLargeError)?;
        let generation_bound = target_max_degree;
        let base_bound = target_max_degree.div_ceil(2);

        let generators = find_generators(&elements);
        let base_semigroup = Self::with_max_degree(elements, grading_vector.clone(), base_bound)?;

        let mut elements_set = matrix_columns_to_set(&base_semigroup.elements);
        let mut reflected_set = HashSet::default();
        for target in target_points.column_iter() {
            for element in elements_set.iter() {
                reflected_set.insert(target - element);
            }
        }

        let mut starting_elements: Vec<_> = elements_set.iter().cloned().collect();
        loop {
            crate::interrupt::check(SemigroupError::Interrupted)?;
            let new_elements = find_new_elements_in_reflected_targets(
                &generators,
                &starting_elements,
                &elements_set,
                &reflected_set,
                &grading_vector,
                generation_bound,
            );
            if new_elements.is_empty() {
                break;
            }
            elements_set.extend(new_elements.iter().cloned());
            starting_elements = new_elements;
        }

        let mut final_elements = HashSet::default();
        for target in target_points.column_iter() {
            for element in elements_set.iter() {
                let reflected = target - element;
                if elements_set.contains(&reflected) {
                    final_elements.insert(element.clone());
                }
            }
        }
        final_elements.insert(DVector::zeros(grading_vector.len()));

        let mut elements = DMatrix::zeros(grading_vector.len(), final_elements.len());
        elements
            .column_iter_mut()
            .zip(final_elements)
            .for_each(|(mut d, s)| d.copy_from(&s));

        Self::from_data(elements, grading_vector)
    }
}

/// The elements of the semigroup generated by the columns of `generators` of degree at most
/// `max_degree`, or of the lowest degrees that hold at least `count` of them, with the
/// modular pipeline's enumeration, which keeps them in packed keys (many times faster than
/// the sets of vectors here); `None` if those cannot hold them.
fn enumerate(
    generators: &DMatrix<i32>,
    grading_vector: &RowDVector<i32>,
    max_degree: u32,
    count: usize,
) -> Option<Result<DMatrix<i32>, SemigroupError>> {
    use crate::modular::candidates::lowest_elements;
    use crate::modular::key::{coordinate_bounds, KeyLayout, PackedLayout};
    let g = grading_vector.as_slice();
    let gens: Vec<&[i32]> = crate::cone::columns(generators).collect();
    let bounds = coordinate_bounds(gens.iter().copied(), g, max_degree).ok()?;
    let layout = PackedLayout::<u128>::new(&bounds).ok()?;
    let points = match lowest_elements(&layout, gens, g, max_degree, count) {
        Ok(points) => points,
        Err(crate::modular::Error::Interrupted) => return Some(Err(SemigroupError::Interrupted)),
        Err(_) => return None,
    };
    let dim = g.len();
    let mut elements = DMatrix::zeros(dim, points.len());
    for (mut c, &(key, _)) in elements.column_iter_mut().zip(&points) {
        layout.unpack(key, c.as_mut_slice());
    }
    Some(Ok(elements))
}

/// Sort the elements by degree and make sure that only the identity has degree
/// zero.
fn sort_elements(elements: &mut DMatrix<i32>, grading_vector: &RowDVector<i32>) -> Vec<i64> {
    // By degree, then lexicographically, so that the order does not depend on how the
    // elements were found: the indices are sorted, in parallel, and the columns copied once.
    let cols: Vec<&[i32]> = crate::cone::columns(elements).collect();
    let degrees: Vec<i64> = cols.iter().map(|&c| degree(grading_vector, c)).collect();
    let mut order: Vec<usize> = (0..cols.len()).collect();
    order.par_sort_unstable_by(|&a, &b| (degrees[a], cols[a]).cmp(&(degrees[b], cols[b])));
    let mut sorted = DMatrix::zeros(elements.nrows(), elements.ncols());
    for (mut c, &i) in sorted.column_iter_mut().zip(&order) {
        c.copy_from_slice(cols[i]);
    }
    let degrees = order.iter().map(|&i| degrees[i]).collect();
    *elements = sorted;
    degrees
}

/// Check that the degrees are positive, except for the first one, which must be
/// zero.
fn check_final_degrees(degrees: Vec<i64>) -> Result<RowDVector<u32>, SemigroupError> {
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
        *d = u32::try_from(*s).map_err(|_| SemigroupError::DegreeTooLargeError)?;
    }

    Ok(final_degrees)
}

fn check_degrees(
    elements: &DMatrix<i32>,
    grading_vector: &RowDVector<i32>,
) -> Result<(), SemigroupError> {
    for c in elements.column_iter() {
        let d = degree(grading_vector, c);
        if d < 0 || (d == 0 && c.iter().any(|x| *x != 0)) {
            return Err(SemigroupError::NonPositiveDegreeError);
        }
    }
    Ok(())
}

/// The degree of `c` under `grading_vector`, in 64 bits, which the products of 32-bit entries
/// and their sums fit.
fn degree<'a>(grading_vector: &RowDVector<i32>, c: impl IntoIterator<Item = &'a i32>) -> i64 {
    grading_vector
        .iter()
        .zip(c)
        .map(|(&g, &x)| g as i64 * x as i64)
        .sum()
}

fn matrix_columns_to_set(elements: &DMatrix<i32>) -> HashSet<DVector<i32>> {
    elements.column_iter().map(|c| c.clone_owned()).collect()
}

/// Returns only the elements with degree up to the specified maximum degree.
fn trim_by_max_deg(
    elements: &DMatrix<i32>,
    grading_vector: &RowDVector<i32>,
    max_degree: u32,
) -> DMatrix<i32> {
    let selected: Vec<_> = elements
        .column_iter()
        .enumerate()
        .filter(|(_, c)| degree(grading_vector, c) <= max_degree as i64)
        .map(|(i, _)| i)
        .collect();
    let mut trimmed = DMatrix::<i32>::zeros(elements.nrows(), selected.len());
    trimmed
        .column_iter_mut()
        .zip(selected)
        .for_each(|(mut c, i)| c.copy_from(&elements.column(i)));
    trimmed
}

/// Tries to find a smaller subset of elements that generates the semigroup. It
/// does not necessarily return the minimal set of generators since doing so is
/// very difficult.
fn find_generators(elements: &DMatrix<i32>) -> DMatrix<i32> {
    // TODO: Need to check if it is worth to do this in parallel.

    // Only elements that are the sum of two others are discarded. Checking sums
    // of three as well makes this cubic in the number of elements, which takes
    // hours on inputs with thousands of lattice points, like the ones CYTools
    // produces, while on such inputs it discards at most a couple more elements
    // out of the thousands that pairs already do.
    let dim = elements.nrows();
    let zero_vec = DVector::<i32>::zeros(dim);
    let mut tmp_vec = zero_vec.clone();

    let mut generators: HashSet<_> = elements.column_iter().collect();
    generators.remove(&zero_vec.as_view());

    let mut to_remove = HashSet::default();

    let candidates: Vec<_> = generators.iter().collect();
    for (i, c1) in candidates.iter().enumerate() {
        for c2 in candidates[i..].iter() {
            tmp_vec.copy_from(*c1);
            tmp_vec += **c2;
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

/// The sums of a starting element and a generator that `keep` takes and that are not among
/// the elements, once each.
fn find_new_elements(
    generators: &DMatrix<i32>,
    starting_elements: &[DVector<i32>],
    elements_set: &HashSet<DVector<i32>>,
    keep: impl Fn(&DVector<i32>) -> bool,
) -> Vec<DVector<i32>> {
    let mut new_elements = HashSet::default();
    let mut tmp_vec = DVector::zeros(generators.nrows());
    for c1 in generators.column_iter() {
        for c2 in starting_elements {
            tmp_vec.copy_from(c2);
            tmp_vec += c1;
            if keep(&tmp_vec) && !elements_set.contains(&tmp_vec) {
                new_elements.insert(tmp_vec.clone());
            }
        }
    }
    new_elements.into_iter().collect()
}

/// Find new elements up to a maximum degree using the generators and the starting elements.
fn find_new_elements_until_max_deg(
    generators: &DMatrix<i32>,
    starting_elements: &[DVector<i32>],
    elements_set: &HashSet<DVector<i32>>,
    grading_vector: &RowDVector<i32>,
    max_degree: u32,
) -> Vec<DVector<i32>> {
    find_new_elements(generators, starting_elements, elements_set, |v| {
        degree(grading_vector, v) <= max_degree as i64
    })
}

/// Find new elements that lie in the reflected target regions.
fn find_new_elements_in_reflected_targets(
    generators: &DMatrix<i32>,
    starting_elements: &[DVector<i32>],
    elements_set: &HashSet<DVector<i32>>,
    reflected_set: &HashSet<DVector<i32>>,
    grading_vector: &RowDVector<i32>,
    max_degree: u32,
) -> Vec<DVector<i32>> {
    find_new_elements(generators, starting_elements, elements_set, |v| {
        degree(grading_vector, v) <= max_degree as i64 && reflected_set.contains(v)
    })
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

        let sg_result = Semigroup::from_data(elements.clone(), grading_vector);
        assert!(sg_result.is_err());
        let e = sg_result.err().unwrap();
        assert_eq!(e, SemigroupError::NonPositiveDegreeError);

        // A second zero is a duplicate, not an element of degree zero.
        let doubled = elements.insert_column(0, 0);
        let e = Semigroup::from_data(doubled, RowDVector::from_row_slice(&[1, 1])).err();
        assert_eq!(e, Some(SemigroupError::DuplicateElementError));
    }

    #[test]
    fn test_semigroup_with_max_degree() {
        let (elements, grading_vector) = example_elements_and_grading_vector();

        let sg_result = Semigroup::with_max_degree(elements.clone(), grading_vector.clone(), 3);
        assert!(sg_result.is_ok());
        let sg = sg_result.unwrap();
        assert_eq!(sg.degrees.len(), 10);
        assert_eq!(
            sg.degrees,
            RowDVector::from_row_slice(&[0, 1, 1, 2, 2, 2, 3, 3, 3, 3])
        );

        let sg_result = Semigroup::with_max_degree(elements, grading_vector, 1);
        assert!(sg_result.is_ok());
        let sg = sg_result.unwrap();
        assert_eq!(sg.degrees.len(), 3);
        assert_eq!(sg.degrees, RowDVector::from_row_slice(&[0, 1, 1]));
    }

    #[test]
    fn test_semigroup_too_wide_for_keys() {
        // Five coordinates of up to 3 * 2^26 take more bits than the packed keys have, so
        // the sets of vectors find the elements, a g + b h with a + b <= 3.
        let big = 1 << 26;
        let elements = DMatrix::from_column_slice(
            6,
            3,
            &[
                0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, big, big, big, big, big,
            ],
        );
        let grading_vector = RowDVector::from_row_slice(&[1, 0, 0, 0, 0, 0]);
        assert!(enumerate(&elements, &grading_vector, 3, usize::MAX).is_none());
        let sg = Semigroup::with_max_degree(elements.clone(), grading_vector.clone(), 3).unwrap();
        assert_eq!(sg.elements.ncols(), 10);
        let b = 3 * big;
        assert_eq!(sg.elements.column(9).as_slice(), &[3, b, b, b, b, b]);
        let sg = Semigroup::with_min_elements(elements, grading_vector, 4).unwrap();
        assert_eq!(sg.elements.ncols(), 6);
    }

    #[test]
    fn test_semigroup_with_min_elements() {
        let (elements, grading_vector) = example_elements_and_grading_vector();

        let sg_result = Semigroup::with_min_elements(elements.clone(), grading_vector.clone(), 11);
        assert!(sg_result.is_ok());
        let sg = sg_result.unwrap();
        assert_eq!(sg.degrees.len(), 15);
        assert_eq!(
            sg.degrees,
            RowDVector::from_row_slice(&[0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4])
        );

        let sg_result = Semigroup::with_min_elements(elements, grading_vector, 3);
        assert!(sg_result.is_ok());
        let sg = sg_result.unwrap();
        assert_eq!(sg.degrees.len(), 3);
        assert_eq!(sg.degrees, RowDVector::from_row_slice(&[0, 1, 1]));
    }

    #[test]
    fn test_semigroup_with_target_points() {
        let generators = DMatrix::from_column_slice(2, 2, &[1, 0, 0, 1]);
        let grading_vector = RowDVector::from_row_slice(&[1, 1]);
        let target_points = DMatrix::from_column_slice(2, 1, &[2, 1]);

        let sg = Semigroup::with_target_points(
            generators.clone(),
            grading_vector.clone(),
            target_points,
        )
        .unwrap();
        let elements = matrix_columns_to_set(&sg.elements);
        let expected = HashSet::from_iter([
            DVector::from_column_slice(&[0, 0]),
            DVector::from_column_slice(&[1, 0]),
            DVector::from_column_slice(&[2, 0]),
            DVector::from_column_slice(&[0, 1]),
            DVector::from_column_slice(&[1, 1]),
            DVector::from_column_slice(&[2, 1]),
        ]);

        assert_eq!(elements, expected);

        let target_points = DMatrix::from_column_slice(2, 1, &[10, 10]);
        let sg = Semigroup::with_target_points(
            generators.clone(),
            grading_vector.clone(),
            target_points,
        )
        .unwrap();
        let elements = matrix_columns_to_set(&sg.elements);

        assert_eq!(elements.len(), 121);

        let target_points = DMatrix::from_column_slice(2, 1, &[0, 10]);
        let sg = Semigroup::with_target_points(
            generators.clone(),
            grading_vector.clone(),
            target_points,
        )
        .unwrap();
        let elements = matrix_columns_to_set(&sg.elements);

        assert_eq!(elements.len(), 11);
    }

    #[test]
    fn test_find_generators() {
        let (elements, _) = example_elements_and_grading_vector();

        let generators = find_generators(&elements);
        assert_eq!(generators.ncols(), 2);
    }
}
