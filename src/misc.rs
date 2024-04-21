//! Miscellaneous
//!
//! This module contains miscellaneous functions needed for some computations.

pub mod error;

use error::MiscError;
use nalgebra::DMatrix;
use std::collections::{HashMap, HashSet};

/// Process the input intersection numbers and find the relevant pairs of indices.
/// For threefolds, the relevant pair of indices are simply all distinct (sorted) pairs.
/// Othersize, the relevant pairs are the indices in the second and third columns,
/// since indices in the first column correspond to reference surfaces.
#[allow(clippy::type_complexity)]
pub fn process_int_nums(
    intnums: DMatrix<i32>,
    is_threefold: bool,
) -> Result<(HashMap<(usize, usize, usize), i32>, HashSet<(usize, usize)>), MiscError> {
    if intnums.ncols() == 0 {
        return Err(MiscError::EmptyIntNums);
    } else if intnums.nrows() != 4 {
        return Err(MiscError::WrongDimIntNums);
    }

    let mut intnum_dict = HashMap::new();
    let mut intnum_idxpairs = HashSet::new();

    for col in intnums.column_iter() {
        if col[3] == 0 {
            continue;
        }
        let mut tmp_vec: Vec<_> = col.iter().take(3).cloned().collect();
        if is_threefold {
            tmp_vec.sort_unstable();
            intnum_idxpairs.insert((tmp_vec[0] as usize, tmp_vec[1] as usize));
            intnum_idxpairs.insert((tmp_vec[0] as usize, tmp_vec[2] as usize));
        } else {
            tmp_vec[1..=2].sort_unstable();
        }
        intnum_idxpairs.insert((tmp_vec[1] as usize, tmp_vec[2] as usize));
        if *tmp_vec.iter().min().unwrap() < 0 {
            return Err(MiscError::NegativeIndex);
        }
        let tmp_tup = (
            tmp_vec[0] as usize,
            tmp_vec[1] as usize,
            tmp_vec[2] as usize,
        );
        if intnum_dict.insert(tmp_tup, col[3]).is_some() {
            return Err(MiscError::RepeatedIdxIntNums);
        }
    }

    Ok((intnum_dict, intnum_idxpairs))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::dmatrix;

    #[test]
    fn test_int_nums() {
        let intnums = dmatrix![ 0, 0, 0;
                                1, 1, 2;
                                1, 2, 3;
                               -1, 3, 2;];
        let result = process_int_nums(intnums.clone(), true);
        assert!(result.is_ok());
        let (intnum_dict, intnum_idxpairs) = result.unwrap();
        assert_eq!(intnum_dict.len(), 3);
        assert_eq!(intnum_idxpairs.len(), 6);

        let result = process_int_nums(intnums.clone(), false);
        assert!(result.is_ok());
        let (intnum_dict, intnum_idxpairs) = result.unwrap();
        assert_eq!(intnum_dict.len(), 3);
        assert_eq!(intnum_idxpairs.len(), 3);
    }
}
