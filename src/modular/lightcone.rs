//! The backward lightcones of target points.
//!
//! The backward lightcone of a target $P$ is the set of curve classes $C$ with $P - C$ in the
//! Mori cone. A summand of such a class is again in it, so the invariants on the union of the
//! lightcones of a few targets only depend on the candidates in it, exactly as with the
//! truncation by degree. The pipeline filters at three points, as cgv
//! does: the candidates, the curves whose $q^C$ the extraction subtracts (a curve outside
//! only changes curves outside), and the output of each pass, since products of curves inside
//! reach curves outside, whose values miss the curves dropped.
//!
//! Membership is tested with whatever describes the cone the candidates come from: the
//! inequalities $h \cdot C \ge 0$ of the Mori cone, as $h \cdot C \le h \cdot P$ for every $h$,
//! when the candidates are the lattice points of the cones $K_T$; or the elements of the
//! semigroup of the generators, when the candidates are those, which agrees with the reference
//! pipeline.

use super::key::KeyLayout;
use super::lattice::Lattice;
use super::table::SlotMap;
use super::{dot, Error};
use crate::cone::{columns, ConeData};

/// The union of the backward lightcones of some targets, or every curve class.
pub enum Lightcone<K> {
    /// Every curve class.
    All,
    /// The classes $C$ with $h \cdot C \le h \cdot P$ for every inequality $h$ of the Mori cone,
    /// for some target $P$.
    Hyperplanes {
        /// The inequalities, each of `dim` entries, one after the other.
        hyperplanes: Vec<i32>,
        /// $h \cdot P$ for each inequality, for each target.
        bounds: Vec<Vec<i64>>,
    },
    /// The classes in a set, the elements $C$ of the semigroup of the generators with $P - C$
    /// in it as well, for some target $P$.
    Points(SlotMap<K, ()>),
}

impl<K: super::key::Key> Lightcone<K> {
    /// The lightcones of `targets` in the cone of `cone_data`. Targets outside the linear span
    /// of the cone have nothing in it below them, so they are left out.
    pub fn from_cone_data(cone_data: &ConeData, targets: &[&[i32]]) -> Result<Self, Error> {
        let dim = cone_data.hyperplanes.nrows();
        if targets.iter().any(|p| p.len() != dim)
            || (cone_data.equations.ncols() > 0 && cone_data.equations.nrows() != dim)
        {
            return Err(Error::DimensionMismatch);
        }
        let hyperplanes: Vec<&[i32]> = columns(&cone_data.hyperplanes).collect();
        let bounds = targets
            .iter()
            .filter(|p| columns(&cone_data.equations).all(|e| dot(e, p) == 0))
            .map(|p| hyperplanes.iter().map(|h| dot(h, p)).collect())
            .collect();
        Ok(Self::Hyperplanes {
            hyperplanes: hyperplanes.concat(),
            bounds,
        })
    }

    /// The lightcones of `targets` in `semigroup`, which must have every element of degree up
    /// to that of the targets, such as the output of
    /// [`generated_semigroup`](super::candidates::generated_semigroup).
    pub fn from_semigroup<L: KeyLayout<Key = K>>(
        layout: &L,
        semigroup: &[(K, u32)],
        targets: &[&[i32]],
    ) -> Result<Self, Error> {
        let dim = layout.dim();
        if targets.iter().any(|p| p.len() != dim) {
            return Err(Error::DimensionMismatch);
        }
        let mut all: SlotMap<K, ()> = SlotMap::new(1, 2 * semigroup.len());
        for &(key, deg) in semigroup {
            all.insert(key, deg);
        }
        let mut points = SlotMap::new(1, 1024);
        let (mut c, mut d) = (vec![0; dim], vec![0; dim]);
        for &(key, deg) in semigroup {
            layout.unpack(key, &mut c);
            let below = targets.iter().any(|p| {
                for ((y, &x), &t) in d.iter_mut().zip(&c).zip(*p) {
                    *y = t - x;
                }
                layout.fits(&d) && all.get(layout.pack(&d)).is_some()
            });
            if below {
                points.insert(key, deg);
            }
        }
        Ok(Self::Points(points))
    }

    /// Whether the curve class with key `key` and coordinates `curve` is in the lightcones.
    ///
    /// The equations of a Mori cone that is a face are not tested: every class asked about is
    /// a candidate or a sum of candidates, which lie in its span.
    pub fn contains(&self, key: K, curve: &[i32]) -> bool {
        match self {
            Self::All => true,
            Self::Hyperplanes {
                hyperplanes,
                bounds,
            } => {
                let dim = curve.len().max(1);
                bounds.iter().any(|hp| {
                    hyperplanes
                        .chunks_exact(dim)
                        .zip(hp)
                        .all(|(h, &b)| dot(h, curve) <= b)
                })
            }
            Self::Points(points) => points.get(key).is_some(),
        }
    }

    /// The same lightcones for curve classes in the coordinates of `lattice`: the inequalities
    /// of the Mori cone restricted to it, with the same bounds. The sets of classes of the
    /// generated semigroup are built in the coordinates of the run in the first place, and
    /// stay as they are.
    pub fn restrict(self, lattice: &Lattice) -> Result<Self, Error> {
        Ok(match self {
            Self::Hyperplanes {
                hyperplanes,
                bounds,
            } => Self::Hyperplanes {
                hyperplanes: hyperplanes
                    .chunks_exact(lattice.h11().max(1))
                    .map(|h| lattice.restrict(h))
                    .collect::<Result<Vec<_>, _>>()?
                    .concat(),
                bounds,
            },
            other => other,
        })
    }

    /// Whether every curve class is in the lightcones.
    pub fn is_all(&self) -> bool {
        matches!(self, Self::All)
    }
}

/// The largest degree of the targets, or zero without any.
pub(crate) fn max_degree(targets: &[&[i32]], grading_vector: &[i32]) -> Result<u32, Error> {
    let mut max = 0;
    for p in targets {
        if p.len() != grading_vector.len() {
            return Err(Error::DimensionMismatch);
        }
        max = max.max(dot(p, grading_vector).clamp(0, u32::MAX as i64) as u32);
    }
    Ok(max)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modular::candidates::generated_semigroup;
    use crate::modular::key::{PackedKey, PackedLayout};
    use nalgebra::{dmatrix, DMatrix};

    #[test]
    fn test_hyperplanes() {
        // The cone spanned by (1, 0) and (1, 2), and its face spanned by (1, 0).
        let cone_data = ConeData {
            hyperplanes: dmatrix![0, 2; 1, -1],
            equations: DMatrix::zeros(2, 0),
            cones: Vec::new(),
        };
        let lc = Lightcone::<PackedKey>::from_cone_data(&cone_data, &[&[2, 2]]).unwrap();
        let below = |c: &[i32]| lc.contains(PackedKey(0), c);
        for c in [[0, 0], [1, 0], [1, 1], [2, 2], [2, 1], [2, 0], [1, 2]] {
            assert_eq!(below(&c), c != [2, 1] && c != [2, 0], "{c:?}");
        }
        // Two targets give the union.
        let lc = Lightcone::<PackedKey>::from_cone_data(&cone_data, &[&[2, 2], &[3, 0]]).unwrap();
        assert!(lc.contains(PackedKey(0), &[2, 0]) && lc.contains(PackedKey(0), &[3, 0]));
        assert!(!lc.contains(PackedKey(0), &[2, 1]));

        let face = ConeData {
            hyperplanes: dmatrix![1; 0],
            equations: dmatrix![0; 1],
            ..cone_data
        };
        let lc = Lightcone::<PackedKey>::from_cone_data(&face, &[&[2, 2], &[3, 0]]).unwrap();
        assert!(lc.contains(PackedKey(0), &[3, 0]) && !lc.contains(PackedKey(0), &[4, 0]));
        // The target off the face has nothing below it.
        let lc = Lightcone::<PackedKey>::from_cone_data(&face, &[&[2, 2]]).unwrap();
        assert!(!lc.contains(PackedKey(0), &[0, 0]));
        assert!(Lightcone::<PackedKey>::from_cone_data(&face, &[&[2]]).is_err());
    }

    #[test]
    fn test_semigroup() {
        // The semigroup of (1, 0) and (1, 2), which misses (1, 1): below (2, 2) are only the
        // elements whose difference with it is an element.
        let layout = PackedLayout::<u128>::new(&[4, 4]).unwrap();
        let sg = generated_semigroup(&layout, [&[1, 0][..], &[1, 2]], &[1, 1], 4).unwrap();
        let lc = Lightcone::from_semigroup(&layout, &sg, &[&[2, 2]]).unwrap();
        let inside: Vec<[i32; 2]> = sg
            .iter()
            .map(|&(key, _)| {
                let mut c = [0; 2];
                layout.unpack(key, &mut c);
                (key, c)
            })
            .filter(|(key, c)| lc.contains(*key, c))
            .map(|(_, c)| c)
            .collect();
        assert_eq!(inside, [[0, 0], [1, 0], [1, 2], [2, 2]]);
        assert_eq!(max_degree(&[&[2, 2], &[3, 0]], &[1, 1]), Ok(4));
        assert_eq!(max_degree(&[&[-2, 0]], &[1, 1]), Ok(0));
        assert!(max_degree(&[&[2]], &[1, 1]).is_err());
    }
}
