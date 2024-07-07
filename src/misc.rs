//! Miscellaneous functions

pub mod error;

use error::MiscError;
use std::collections::{HashMap, HashSet};

/// Process the input intersection numbers and find the relevant pairs of indices.
/// For threefolds, the relevant pair of indices are simply all distinct (sorted) pairs.
/// Othersize, the relevant pairs are the indices in the second and third columns,
/// since indices in the first column correspond to reference surfaces.
#[allow(clippy::type_complexity)]
pub fn process_int_nums(
    intnums: HashMap<(usize, usize, usize), i32>,
    is_threefold: bool,
) -> Result<
    (
        HashMap<(usize, usize, usize), i32>,
        HashSet<(usize, usize)>,
        usize,
    ),
    MiscError,
> {
    if intnums.is_empty() {
        return Err(MiscError::EmptyIntNums);
    }

    let mut intnum_res = HashMap::new();
    let mut intnum_idxpairs = HashSet::new();

    for (idx, val) in intnums.iter() {
        if *val == 0 {
            continue;
        }
        let (i, j, k) = *idx;
        let mut tmp_vec = [i, j, k];
        if is_threefold {
            tmp_vec.sort_unstable();
            intnum_idxpairs.insert((tmp_vec[0], tmp_vec[1]));
            intnum_idxpairs.insert((tmp_vec[0], tmp_vec[2]));
        } else {
            tmp_vec[1..=2].sort_unstable();
        }
        intnum_idxpairs.insert((tmp_vec[1], tmp_vec[2]));
        let tmp_tup = (tmp_vec[0], tmp_vec[1], tmp_vec[2]);
        if intnum_res.insert(tmp_tup, *val).is_some() {
            return Err(MiscError::RepeatedIdxIntNums);
        }
    }
    let n_indices = if is_threefold {
        intnum_idxpairs.iter().map(|p| p.1).max().unwrap()
    } else {
        intnum_res.keys().map(|p| p.0).max().unwrap()
    } + 1;

    Ok((intnum_res, intnum_idxpairs, n_indices))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_int_nums() {
        let intnums = HashMap::from([((0, 1, 1), -1), ((0, 1, 2), 3), ((0, 2, 3), 2)]);
        let result = process_int_nums(intnums.clone(), true);
        assert!(result.is_ok());
        let (intnum_dict, intnum_idxpairs, n_indices) = result.unwrap();
        assert_eq!(intnum_dict.len(), 3);
        assert_eq!(intnum_idxpairs.len(), 6);
        assert_eq!(n_indices, 4);

        let result = process_int_nums(intnums.clone(), false);
        assert!(result.is_ok());
        let (intnum_dict, intnum_idxpairs, n_indices) = result.unwrap();
        assert_eq!(intnum_dict.len(), 3);
        assert_eq!(intnum_idxpairs.len(), 3);
        assert_eq!(n_indices, 1);
    }
}
