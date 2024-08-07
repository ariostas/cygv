use cygv::hkty::{
    compute_gv_float_nfold, compute_gv_float_threefold, compute_gv_rat_nfold,
    compute_gv_rat_threefold, compute_gw_float_nfold, compute_gw_float_threefold,
    compute_gw_rat_nfold, compute_gw_rat_threefold,
};
use nalgebra::{dmatrix, dvector, RowDVector};
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
        q.clone(),
        nefpart.clone(),
        intnums.clone(),
        None,
        1000,
    );
    compute_gw_rat_threefold(
        generators.clone(),
        grading_vector.clone(),
        None,
        Some(100),
        q.clone(),
        nefpart.clone(),
        intnums.clone(),
        None,
        1000,
    );
    compute_gv_float_threefold(
        generators.clone(),
        grading_vector.clone(),
        None,
        Some(100),
        q.clone(),
        nefpart.clone(),
        intnums.clone(),
        None,
        1000,
        200,
    );
    compute_gw_float_threefold(
        generators.clone(),
        grading_vector.clone(),
        None,
        Some(100),
        q.clone(),
        nefpart.clone(),
        intnums.clone(),
        None,
        1000,
        200,
    );
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
        q.clone(),
        nefpart.clone(),
        intnums.clone(),
        None,
        1000,
    );
    compute_gw_rat_nfold(
        generators.clone(),
        grading_vector.clone(),
        None,
        Some(100),
        q.clone(),
        nefpart.clone(),
        intnums.clone(),
        None,
        1000,
    );
    compute_gv_float_nfold(
        generators.clone(),
        grading_vector.clone(),
        None,
        Some(100),
        q.clone(),
        nefpart.clone(),
        intnums.clone(),
        None,
        1000,
        200,
    );
    compute_gw_float_nfold(
        generators.clone(),
        grading_vector.clone(),
        None,
        Some(100),
        q.clone(),
        nefpart.clone(),
        intnums.clone(),
        None,
        1000,
        200,
    );
}
