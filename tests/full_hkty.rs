use cygv::dispatch::{
    compute_gv_rat_nfold, compute_gv_rat_threefold, compute_gw_rat_nfold, compute_gw_rat_threefold,
    ComputeError,
};
use cygv::malachite::Rational;
use cygv::Method;
use nalgebra::{dmatrix, dvector, DMatrix, DVector, RowDVector};
use std::collections::HashMap;

#[test]
fn test_threefold() {
    let generators = dmatrix![
         0, 1;
        -1, 2;
    ];
    let grading_vector = RowDVector::from_row_slice(&[3, -1]);

    let q = dmatrix![
        1,  0;
        1,  0;
        1, -1;
        0,  1;
        1,  1;
        2, -1;
    ];
    let nefpart = Vec::new();

    let intnums = HashMap::from([
        ((0, 0, 0), 2),
        ((0, 0, 1), 1),
        ((0, 1, 1), -1),
        ((1, 1, 1), -5),
    ]);

    // For now, these are just smoke tests
    compute_gv_rat_threefold(
        generators.clone(),
        grading_vector.clone(),
        None,
        Some(100),
        None,
        q.clone(),
        nefpart.clone(),
        intnums.clone(),
        Vec::new(),
        None,
        None,
        None,
        Method::Reference,
        false,
    )
    .unwrap();
    compute_gw_rat_threefold(
        generators.clone(),
        grading_vector.clone(),
        None,
        Some(100),
        None,
        q.clone(),
        nefpart.clone(),
        intnums.clone(),
        Vec::new(),
        None,
        None,
        None,
        Method::Reference,
        false,
    )
    .unwrap();
}

#[test]
fn test_fourfold() {
    let generators = dmatrix![
        0,1,0,0,0,0;
        0,0,0,1,-2,0;
        0,-4,1,0,0,0;
        -2,0,1,0,1,0;
        3,0,0,1,0,-2;
        0,0,0,-2,2,1;
    ];
    let grading_vector = RowDVector::from_row_slice(&[73, -30, 18, -17, -11, -21]);

    let q = dmatrix![
        1,0,0,0,0,0;
        1,0,0,0,0,0;
        1,0,0,0,0,0;
        1,0,0,0,0,0;
        0,1,0,0,0,0;
        0,0,1,0,0,0;
        8,-4,2,-2,-2,-3;
        12,-6,3,-3,-2,-4;
        4,-3,0,-2,-1,-2;
        0,0,0,1,0,0;
        0,0,0,0,1,0;
        0,0,0,0,0,1;
    ];
    let nefpart = vec![dvector![0, 1, 2, 3], dvector![4, 5, 6, 7, 8, 9, 10, 11]];

    let intnums = HashMap::from([
        ((0, 0, 2), 1),
        ((0, 1, 1), -8),
        ((0, 1, 3), 4),
        ((0, 1, 5), 8),
        ((0, 2, 2), -4),
        ((0, 3, 3), -8),
        ((0, 4, 4), -16),
        ((0, 4, 5), 8),
        ((0, 5, 5), -16),
        ((1, 0, 1), -8),
        ((1, 0, 3), 4),
        ((1, 0, 5), 8),
        ((1, 1, 3), 16),
        ((1, 1, 5), -32),
        ((1, 3, 3), -16),
        ((1, 5, 5), 64),
        ((2, 0, 0), 1),
        ((2, 0, 2), -4),
        ((2, 2, 2), 16),
        ((3, 0, 1), 4),
        ((3, 0, 3), -8),
        ((3, 1, 1), 16),
        ((3, 1, 3), -16),
        ((4, 0, 4), -16),
        ((4, 0, 5), 8),
        ((4, 4, 4), -128),
        ((4, 4, 5), 32),
        ((5, 0, 1), 8),
        ((5, 0, 4), 8),
        ((5, 0, 5), -16),
        ((5, 1, 1), -32),
        ((5, 1, 5), 64),
        ((5, 4, 4), 32),
        ((5, 5, 5), -128),
        ((6, 0, 0), -8),
        ((6, 0, 3), 16),
        ((6, 0, 5), -32),
        ((6, 1, 1), -128),
        ((6, 1, 3), 64),
        ((6, 1, 5), 128),
        ((6, 3, 3), -64),
        ((6, 5, 5), -256),
        ((7, 0, 0), 4),
        ((7, 0, 1), 16),
        ((7, 0, 3), -16),
        ((7, 1, 1), 64),
        ((7, 1, 3), -64),
        ((7, 3, 3), 64),
        ((8, 0, 0), 8),
        ((8, 0, 1), -32),
        ((8, 0, 5), 64),
        ((8, 1, 1), 128),
        ((8, 1, 5), -256),
        ((8, 5, 5), 512),
        ((9, 0, 0), -4),
        ((9, 0, 2), 16),
        ((9, 2, 2), -64),
        ((10, 0, 0), -8),
        ((10, 0, 1), -16),
        ((10, 1, 1), -64),
        ((10, 1, 3), 64),
        ((10, 3, 3), -128),
        ((11, 0, 0), -16),
        ((11, 0, 4), -128),
        ((11, 0, 5), 32),
        ((11, 4, 4), -768),
        ((11, 4, 5), 128),
        ((12, 0, 0), 8),
        ((12, 0, 4), 32),
        ((12, 4, 4), 128),
        ((13, 0, 0), -16),
        ((13, 0, 1), 64),
        ((13, 0, 5), -128),
        ((13, 1, 1), -256),
        ((13, 1, 5), 512),
        ((13, 5, 5), -1024),
    ]);

    // For now, these are just smoke tests
    compute_gv_rat_nfold(
        generators.clone(),
        grading_vector.clone(),
        None,
        Some(100),
        None,
        q.clone(),
        nefpart.clone(),
        intnums.clone(),
        Vec::new(),
        None,
        None,
        None,
        Method::Reference,
        false,
    )
    .unwrap();
    compute_gw_rat_nfold(
        generators.clone(),
        grading_vector.clone(),
        None,
        Some(100),
        None,
        q.clone(),
        nefpart.clone(),
        intnums.clone(),
        Vec::new(),
        None,
        None,
        None,
        Method::Reference,
        false,
    )
    .unwrap();

    // Both methods give the same invariants of this complete intersection.
    let gv = |method| {
        let mut res = compute_gv_rat_nfold(
            generators.clone(),
            grading_vector.clone(),
            Some(8),
            None,
            None,
            q.clone(),
            nefpart.clone(),
            intnums.clone(),
            Vec::new(),
            None,
            None,
            None,
            method,
            false,
        )
        .unwrap();
        res.sort_unstable_by(|a, b| (a.0 .0.as_slice(), a.0 .1).cmp(&(b.0 .0.as_slice(), b.0 .1)));
        res
    };
    let gw = |method| {
        let mut res = compute_gw_rat_nfold(
            generators.clone(),
            grading_vector.clone(),
            Some(8),
            None,
            None,
            q.clone(),
            nefpart.clone(),
            intnums.clone(),
            Vec::new(),
            None,
            None,
            None,
            method,
            false,
        )
        .unwrap();
        res.sort_unstable_by(|a, b| (a.0 .0.as_slice(), a.0 .1).cmp(&(b.0 .0.as_slice(), b.0 .1)));
        res
    };
    let expected = gv(Method::Reference);
    assert!(!expected.is_empty());
    assert_eq!(gv(Method::Modular), expected);
    assert_eq!(gv(Method::Auto), expected);
    let expected = gw(Method::Reference);
    assert!(expected.len() > gv(Method::Reference).len());
    assert_eq!(gw(Method::Modular), expected);
}

/// The quintic, through the entry points with each method.
#[test]
fn test_methods() {
    let generators = dmatrix![1];
    let grading_vector = RowDVector::from_row_slice(&[1]);
    let q = DMatrix::from_element(5, 1, 1);
    let intnums = HashMap::from([((0, 0, 0), 5)]);
    let gv = |method, max_deg, min_points| {
        compute_gv_rat_threefold(
            generators.clone(),
            grading_vector.clone(),
            max_deg,
            min_points,
            None,
            q.clone(),
            Vec::new(),
            intnums.clone(),
            Vec::new(),
            None,
            None,
            None,
            method,
            false,
        )
        .map(|mut res| {
            res.sort_unstable_by(|a, b| a.0.as_slice().cmp(b.0.as_slice()));
            res
        })
    };
    let gw = |method| {
        let mut res = compute_gw_rat_threefold(
            generators.clone(),
            grading_vector.clone(),
            Some(6),
            None,
            None,
            q.clone(),
            Vec::new(),
            intnums.clone(),
            Vec::new(),
            None,
            None,
            None,
            method,
            false,
        )
        .unwrap();
        res.sort_unstable_by(|a, b| a.0.as_slice().cmp(b.0.as_slice()));
        res
    };

    let expected = gv(Method::Reference, Some(6), None).unwrap();
    assert_eq!(expected.len(), 6);
    assert_eq!(expected[0], (dvector![1], 2875.into()));
    assert_eq!(gv(Method::Modular, Some(6), None).unwrap(), expected);
    assert_eq!(gv(Method::Auto, Some(6), None).unwrap(), expected);
    let expected = gw(Method::Reference);
    // GV_2 + GV_1 / 8
    assert_eq!(expected[1].1, Rational::from_signeds(609250 * 8 + 2875, 8));
    assert_eq!(gw(Method::Modular), expected);

    // The modular method only supports max_deg; auto falls back to the reference one.
    assert_eq!(
        gv(Method::Modular, None, Some(4)),
        Err(ComputeError::Unsupported("min_points"))
    );
    let expected = gv(Method::Reference, None, Some(4)).unwrap();
    assert!(!expected.is_empty());
    assert_eq!(gv(Method::Auto, None, Some(4)).unwrap(), expected);
}

/// Invalid nef partitions, intersection numbers and untruncated generators are errors with
/// each method, rather than panics of the reference, and untruncated generators get zero.
#[test]
fn test_invalid_input() {
    let q = DMatrix::from_element(5, 1, 1);
    let quintic = HashMap::from([((0, 0, 0), 5)]);
    let gv = |generators: DMatrix<i32>,
              max_deg,
              nefpart: Vec<DVector<i32>>,
              intnums: &HashMap<(usize, usize, usize), i32>,
              method| {
        compute_gv_rat_threefold(
            generators,
            RowDVector::from_row_slice(&[1]),
            max_deg,
            None,
            None,
            q.clone(),
            nefpart,
            intnums.clone(),
            Vec::new(),
            None,
            None,
            None,
            method,
            false,
        )
        .map(|mut res| {
            res.sort_unstable_by(|a, b| a.0.as_slice().cmp(b.0.as_slice()));
            res
        })
    };
    for method in [Method::Auto, Method::Modular, Method::Reference] {
        let run =
            |nefpart, intnums: &HashMap<_, _>| gv(dmatrix![1], Some(3), nefpart, intnums, method);
        for (nefpart, reason) in [
            (vec![dvector![0, 1, 2, 3]], "include every row"),
            (vec![dvector![0, 1, 2, 3, 4, 4]], "disjoint"),
            (vec![dvector![0, 1, 2, 3, 7]], "indices of rows"),
            (vec![dvector![0, 1, 2, 3, 4], dvector![]], "not be empty"),
        ] {
            match run(nefpart, &quintic) {
                Err(ComputeError::InvalidNefPartition(r)) => assert!(r.contains(reason), "{r}"),
                res => panic!("{method:?}, {reason}: {res:?}"),
            }
        }
        assert_eq!(
            run(vec![dvector![0, 1], dvector![2, 3, 4]], &quintic),
            Err(ComputeError::InvalidDimension(2))
        );
        // A partition into one part is the hypersurface.
        assert_eq!(
            run(vec![dvector![4, 0, 1, 2, 3]], &quintic),
            run(Vec::new(), &quintic)
        );
        for intnums in [
            HashMap::from([((0, 0, 0), 5), ((0, 0, 1), 1)]),
            HashMap::from([((0, 0, 0), 0)]),
        ] {
            assert!(
                matches!(
                    run(Vec::new(), &intnums),
                    Err(ComputeError::InvalidIntersectionNumbers(_))
                ),
                "{method:?}: {intnums:?}"
            );
        }
    }

    // Without truncation the generators are the classes, with zero added, each once.
    let expected = gv(
        dmatrix![1],
        Some(2),
        Vec::new(),
        &quintic,
        Method::Reference,
    )
    .unwrap();
    assert_eq!(expected.len(), 2);
    for generators in [dmatrix![1, 2], dmatrix![2, 0, 1], dmatrix![1, 2, 1]] {
        let res = gv(generators, None, Vec::new(), &quintic, Method::Auto);
        assert_eq!(res.unwrap(), expected);
    }
    assert!(matches!(
        gv(dmatrix![1, -1], None, Vec::new(), &quintic, Method::Auto),
        Err(ComputeError::InvalidGenerators(_))
    ));
}
