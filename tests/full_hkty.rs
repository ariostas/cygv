use cygv::hkty::{
    compute_gv_float_nfold, compute_gv_float_threefold, compute_gv_rat_nfold,
    compute_gv_rat_threefold, compute_gw_float_nfold, compute_gw_float_threefold,
    compute_gw_rat_nfold, compute_gw_rat_threefold,
};
use nalgebra::{dmatrix, dvector, RowDVector};

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

    let intnums = dmatrix![
        0, 0,  0,  1;
        0, 0,  1,  1;
        0, 1,  1,  1;
        2, 1, -1, -5;
    ];

    // For now, these are just smoke tests
    compute_gv_rat_threefold(
        generators.clone(),
        grading_vector.clone(),
        Some(20),
        None,
        q.clone(),
        nefpart.clone(),
        intnums.clone(),
    );
    // compute_gw_rat_threefold(
    //     generators.clone(),
    //     grading_vector.clone(),
    //     None,
    //     Some(100),
    //     q.clone(),
    //     nefpart.clone(),
    //     intnums.clone(),
    // );
    // compute_gv_float_threefold(
    //     generators.clone(),
    //     grading_vector.clone(),
    //     None,
    //     Some(100),
    //     q.clone(),
    //     nefpart.clone(),
    //     intnums.clone(),
    //     200,
    // );
    // compute_gw_float_threefold(
    //     generators.clone(),
    //     grading_vector.clone(),
    //     None,
    //     Some(100),
    //     q.clone(),
    //     nefpart.clone(),
    //     intnums.clone(),
    //     200,
    // );
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

    let intnums = dmatrix![
        0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,13,13,13,13,13,13;
        0,1,1,1,2,3,4,4,5,0,0,0,1,1,3,5,0,0,2,0,0,1,1,0,0,4,4,0,0,0,1,1,4,5,0,0,0,1,1,1,3,5,0,0,0,1,1,3,0,0,0,1,1,5,0,0,2,0,0,1,1,3,0,0,0,4,4,0,0,4,0,0,0,1,1,5;
        2,1,3,5,2,3,4,5,5,1,3,5,3,5,3,5,0,2,2,1,3,1,3,4,5,4,5,1,4,5,1,5,4,5,0,3,5,1,3,5,3,5,0,1,3,1,3,3,0,1,5,1,5,5,0,2,2,0,1,1,3,3,0,4,5,4,5,0,4,4,0,1,5,1,5,5;
        1,-8,4,8,-4,-8,-16,8,-16,-8,4,8,16,-32,-16,64,1,-4,16,4,-8,16,-16,-16,8,-128,32,8,8,-16,-32,64,32,-128,-8,16,-32,-128,64,128,-64,-256,4,16,-16,64,-64,64,8,-32,64,128,-256,512,-4,16,-64,-8,-16,-64,64,-128,-16,-128,32,-768,128,8,32,128,-16,64,-128,-256,512,-1024;
    ];

    // For now, these are just smoke tests
    compute_gv_rat_nfold(
        generators.clone(),
        grading_vector.clone(),
        None,
        Some(100),
        q.clone(),
        nefpart.clone(),
        intnums.clone(),
    );
    compute_gw_rat_nfold(
        generators.clone(),
        grading_vector.clone(),
        None,
        Some(100),
        q.clone(),
        nefpart.clone(),
        intnums.clone(),
    );
    compute_gv_float_nfold(
        generators.clone(),
        grading_vector.clone(),
        None,
        Some(100),
        q.clone(),
        nefpart.clone(),
        intnums.clone(),
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
        200,
    );
}
