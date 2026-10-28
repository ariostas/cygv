//! Phases of the ambient toric variety from its fan (a port of `cgv/tools/cgv_phase.py` of
//! cgv).
//!
//! A threefold hypersurface is given by the charges $Q$ and the intersection numbers alone, but
//! in a vex phase (MacFadden–Sheridan, arXiv:2512.14817) the invariants also need the strata
//! (see [`Stratum`]), which come from the fan: the 3-cones whose rays share no facet of the
//! convex hull of the rays. [`Fan`] derives them, together with the Mori cone (spanned by the
//! classes of the curves $V(\tau)$ of the walls $\tau$, the 3-cones) and a grading vector,
//! from the maximal cones and the charges, and checks intersection numbers against the fan's
//! intersection ring. The rays are the integer Gale dual of the charges, which fixes them up
//! to $GL(4, \mathbb{Z})$, under which everything here is invariant.
//!
//! The fan must be what the strata formula assumes: a complete simplicial fan in dimension
//! four whose rays are the boundary points of a reflexive polytope that are not interior to
//! facets. What can be checked from the fan is checked: the cones, completeness (every 3-cone
//! is a face of exactly two maximal cones), reflexivity, that no ray is interior to a facet,
//! that every vex face is a 3-cone (vex 2-cones are impossible for fine fans,
//! MacFadden–Sheridan Prop. 5), and regularity, through the Mori cone, which contains a line
//! otherwise. That the fan uses every boundary point not interior to a facet is not checked.

use crate::cone::{self, columns, det, kernel_basis, rank, Inequalities};
use crate::Stratum;
use malachite::base::num::arithmetic::traits::Abs;
use malachite::{Integer, Rational};
use nalgebra::{DMatrix, RowDVector};
use std::collections::{HashMap, HashSet};
use std::fmt;

/// The largest entry of a ray. The determinant of four rays by fraction-free elimination
/// multiplies two minors of three, each at most $3^{3/2} B^3$, which for $B = 2^{20}$ is
/// below $2^{125}$.
const MAX_RAY: i64 = 1 << 20;

/// Why the fan or the intersection numbers were rejected.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FanError {
    /// The charges do not describe a toric variety of dimension four, i.e. there are not
    /// $h^{1,1} + 4$ divisors with charges of rank $h^{1,1}$.
    NotFourfold,
    /// A maximal cone is not four distinct divisors whose rays are linearly independent, or
    /// appears twice.
    InvalidCone(Vec<usize>),
    /// The fan is not complete: a 3-cone is a face of other than two maximal cones, or a ray is
    /// in none.
    NotComplete,
    /// The convex hull of the rays is not a reflexive polytope.
    NotReflexive,
    /// A ray is interior to a facet of the convex hull of the rays.
    InteriorToFacet(usize),
    /// A ray is interior to the convex hull of the rays, on none of its facets.
    InteriorToPolytope(usize),
    /// A 2-cone of the fan is vex, which is impossible for a fine fan on a reflexive polytope.
    Vex2Cone([usize; 2]),
    /// The Mori cone contains a line, i.e. the fan is not regular.
    NotRegular,
    /// The grading vector does not have $h^{1,1}$ entries, or is not positive on the Mori cone:
    /// on this class of a wall, if it has the entries.
    NonpositiveGrading(Option<Vec<i32>>),
    /// The intersection numbers have an index above $h^{1,1}$, or two of them, of the same
    /// divisors, differ.
    InvalidIntnums,
    /// The intersection numbers do not match those of the fan, e.g. $D_i D_j D_k$ for the
    /// divisors (rows of the charges) $i, j, k$.
    IntnumsMismatch([usize; 3]),
    /// An entry of the result, or an intermediate one, does not fit its integer type.
    TooLarge,
    /// A fan was given for something other than a threefold hypersurface (returned by the
    /// entry points, such as [`crate::compute`]).
    NotHypersurface,
    /// The computation was interrupted ([`crate::Interrupt`]).
    Interrupted,
}

impl fmt::Display for FanError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::NotFourfold => write!(
                f,
                "the charges must describe a toric fourfold (h11 + 4 divisors, rank h11)"
            ),
            Self::InvalidCone(c) => write!(
                f,
                "the cone {c:?} is not four distinct divisors with independent rays, or repeats"
            ),
            Self::NotComplete => write!(f, "the fan is not complete"),
            Self::NotReflexive => write!(f, "the rays are not on a reflexive polytope"),
            Self::InteriorToFacet(i) => write!(
                f,
                "the ray of divisor {i} is interior to a facet of the polytope of the rays"
            ),
            Self::InteriorToPolytope(i) => write!(
                f,
                "the ray of divisor {i} is interior to the polytope of the rays"
            ),
            Self::Vex2Cone(s) => write!(
                f,
                "the 2-cone {s:?} is vex, which is impossible for a fine fan on a reflexive \
                 polytope (MacFadden-Sheridan, arXiv:2512.14817, Prop. 5); check that the fan \
                 is fine and uses every boundary point not interior to a facet"
            ),
            Self::NotRegular => write!(f, "the fan is not regular (its Mori cone has a line)"),
            Self::NonpositiveGrading(class) => {
                write!(
                    f,
                    "the grading vector must have h11 entries and be positive on the Mori cone, \
                     which has none if the fan is not regular"
                )?;
                match class {
                    Some(c) => write!(f, "; it is not positive on the class {c:?} of a wall"),
                    None => Ok(()),
                }
            }
            Self::InvalidIntnums => write!(
                f,
                "the intersection numbers must have indices below h11 and agree on permutations"
            ),
            Self::IntnumsMismatch([i, j, k]) => write!(
                f,
                "the intersection numbers do not match the fan's intersection ring (e.g. at the \
                 divisors {i}, {j}, {k}); check that the charges and the intersection numbers \
                 use the same basis"
            ),
            Self::TooLarge => write!(f, "the entries of the fan's data are too large"),
            Self::NotHypersurface => {
                write!(f, "a fan is only supported for threefold hypersurfaces")
            }
            Self::Interrupted => write!(f, "the computation was interrupted"),
        }
    }
}

impl std::error::Error for FanError {}

impl From<cone::TooLarge> for FanError {
    fn from(_: cone::TooLarge) -> Self {
        Self::TooLarge
    }
}

impl From<cone::HullError> for FanError {
    fn from(e: cone::HullError) -> Self {
        match e {
            cone::HullError::TooLarge => Self::TooLarge,
            cone::HullError::Interrupted => Self::Interrupted,
        }
    }
}

/// The Mori cone of the fan's toric variety, from [`Fan::mori_cone`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MoriCone {
    /// Its inequalities (no equations, since the fan is regular).
    pub inequalities: Inequalities,
    /// The sum of its facets, which is positive on it.
    pub grading_vector: RowDVector<i32>,
}

/// A complete simplicial fan in dimension four, given by its maximal cones and the charges of
/// its rays.
#[derive(Debug, Clone)]
pub struct Fan {
    /// The charges, one row per divisor.
    q: Vec<Vec<i64>>,
    rays: Vec<[i64; 4]>,
    /// The maximal cones, sorted.
    cones: Vec<[usize; 4]>,
    /// Every face of a maximal cone, the cones themselves and $\{0\}$ included, sorted.
    faces: HashSet<Vec<usize>>,
    /// A left inverse of the charges, from [`left_inverse`].
    inverse: Vec<(usize, Vec<Rational>)>,
    /// The divisors that share a 2-cone with each divisor, sorted.
    neighbors: Vec<Vec<usize>>,
}

impl Fan {
    /// The fan with maximal cones `cones`, each four indices of rows of the charge matrix `q`
    /// (one row per divisor, one column per Kähler parameter), checked as the module
    /// documentation describes.
    pub fn new(cones: &[Vec<usize>], q: &DMatrix<i32>) -> Result<Self, FanError> {
        let (n, h11) = q.shape();
        let charges: Vec<Vec<i64>> = columns(q)
            .map(|c| c.iter().map(|&x| x as i64).collect())
            .collect();
        if h11 == 0 || n != h11 + 4 || rank(&charges) != h11 {
            return Err(FanError::NotFourfold);
        }
        let kernel = kernel_basis(&charges, n)?;
        let rays: Vec<[i64; 4]> = (0..n)
            .map(|i| std::array::from_fn(|k| kernel[k][i]))
            .collect();
        // Bounded so that they fit an `i32` and the determinants of four of them an `i128`.
        if rays.iter().flatten().any(|x| x.abs() > MAX_RAY) {
            return Err(FanError::TooLarge);
        }
        let q: Vec<Vec<i64>> = (0..n)
            .map(|i| charges.iter().map(|r| r[i]).collect())
            .collect();
        let inverse = left_inverse(&q, h11);

        let mut sorted = Vec::with_capacity(cones.len());
        let mut faces = HashSet::new();
        for c in cones {
            let mut s: [usize; 4] = c
                .as_slice()
                .try_into()
                .map_err(|_| FanError::InvalidCone(c.clone()))?;
            s.sort_unstable();
            if s[3] >= n || s.windows(2).any(|w| w[0] == w[1]) || !faces.insert(s.to_vec()) {
                return Err(FanError::InvalidCone(c.clone()));
            }
            let rows: Vec<Vec<i64>> = s.iter().map(|&i| rays[i].to_vec()).collect();
            if det(&rows)? == 0 {
                return Err(FanError::InvalidCone(c.clone()));
            }
            sorted.push(s);
        }
        let mut walls: HashMap<Vec<usize>, usize> = HashMap::new();
        for s in &sorted {
            for mask in 0u32..15 {
                let face: Vec<usize> = (0..4)
                    .filter(|k| mask >> k & 1 == 1)
                    .map(|k| s[k])
                    .collect();
                if face.len() == 3 {
                    *walls.entry(face.clone()).or_default() += 1;
                }
                faces.insert(face);
            }
        }
        if walls.values().any(|&k| k != 2) || (0..n).any(|i| !faces.contains(&vec![i])) {
            return Err(FanError::NotComplete);
        }
        let mut neighbors = vec![Vec::new(); n];
        for f in faces.iter().filter(|f| f.len() == 2) {
            neighbors[f[0]].push(f[1]);
            neighbors[f[1]].push(f[0]);
        }
        neighbors.iter_mut().for_each(|v| v.sort_unstable());
        let fan = Self {
            q,
            rays,
            cones: sorted,
            faces,
            inverse,
            neighbors,
        };
        fan.polytope_facets()?;
        Ok(fan)
    }

    /// The number of divisors.
    fn n(&self) -> usize {
        self.rays.len()
    }

    fn h11(&self) -> usize {
        self.n() - 4
    }

    /// The rays, one per column, in the basis of the Gale dual of the charges.
    pub fn rays(&self) -> DMatrix<i32> {
        DMatrix::from_fn(4, self.n(), |k, i| self.rays[i][k] as i32)
    }

    /// The walls (3-cones) of the fan, sorted.
    pub fn walls(&self) -> Vec<[usize; 3]> {
        let mut walls: Vec<[usize; 3]> = self
            .faces
            .iter()
            .filter_map(|f| f.as_slice().try_into().ok())
            .collect();
        walls.sort_unstable();
        walls
    }

    /// The facets of the convex hull of the rays containing each ray, after checking that it
    /// is reflexive and that no ray is interior to a facet.
    fn polytope_facets(&self) -> Result<Vec<Vec<usize>>, FanError> {
        let points = DMatrix::from_fn(
            5,
            self.n(),
            |k, i| {
                if k < 4 {
                    self.rays[i][k] as i32
                } else {
                    1
                }
            },
        );
        let Inequalities {
            hyperplanes,
            equations,
        } = cone::facets(&points)?;
        // Each facet is h . v + h_4 >= 0 with h primitive, which for a polytope with the
        // origin in its interior is at lattice distance h_4 from it.
        if equations.ncols() > 0 || hyperplanes.row(4).iter().any(|&x| x != 1) {
            return Err(FanError::NotReflexive);
        }
        let on: Vec<Vec<usize>> = columns(&points)
            .map(|v| {
                columns(&hyperplanes)
                    .enumerate()
                    .filter(|(_, h)| {
                        h.iter()
                            .zip(v)
                            .map(|(&a, &b)| a as i64 * b as i64)
                            .sum::<i64>()
                            == 0
                    })
                    .map(|(k, _)| k)
                    .collect()
            })
            .collect();
        if let Some(i) = on.iter().position(Vec::is_empty) {
            return Err(FanError::InteriorToPolytope(i));
        }
        if let Some(i) = on.iter().position(|f| f.len() < 2) {
            return Err(FanError::InteriorToFacet(i));
        }
        Ok(on)
    }

    /// The faces of two or three rays that share no facet of the convex hull of the rays,
    /// sorted by size and then lexicographically.
    pub fn vex_cones(&self) -> Vec<Vec<usize>> {
        let on = self
            .polytope_facets()
            .expect("the polytope was checked when the fan was built");
        let mut vex: Vec<Vec<usize>> = self
            .faces
            .iter()
            .filter(|f| matches!(f.len(), 2 | 3))
            .filter(|f| {
                !on[f[0]]
                    .iter()
                    .any(|k| f[1..].iter().all(|&i| on[i].contains(k)))
            })
            .cloned()
            .collect();
        vex.sort_unstable_by(|a, b| (a.len(), a).cmp(&(b.len(), b)));
        vex
    }

    /// The class of the curve $V(\tau)$ of the wall $\tau$, as its intersections $J_a \cdot
    /// V(\tau)$ with the basis of the charges (Cox–Little–Schenck, Prop. 6.4.4).
    pub fn wall_class(&self, wall: [usize; 3]) -> Vec<Rational> {
        let mut apex = self
            .cones
            .iter()
            .filter(|c| wall.iter().all(|i| c.contains(i)))
            .map(|c| *c.iter().find(|i| !wall.contains(i)).unwrap());
        let (u1, u2) = (apex.next().unwrap(), apex.next().unwrap());
        // The relation v_u1 + l_u2 v_u2 + sum_s l_s v_s = 0, unique since u2 and the wall span
        // the space.
        let others = [u2, wall[0], wall[1], wall[2]];
        let a: Vec<Vec<Rational>> = (0..4)
            .map(|k| {
                others
                    .iter()
                    .map(|&i| Rational::from(self.rays[i][k]))
                    .collect()
            })
            .collect();
        let b: Vec<Rational> = (0..4).map(|k| Rational::from(-self.rays[u1][k])).collect();
        let lambda = solve(&a, &b, 4).expect("a maximal cone spans the space");
        // D_u1 . V(tau) = mult(tau) / mult(sigma_1), with mult(tau) the gcd of the 3x3 minors.
        let mut mult = 0i128;
        for skip in 0..4 {
            let minor: Vec<Vec<i64>> = wall
                .iter()
                .map(|&i| {
                    (0..4)
                        .filter(|&k| k != skip)
                        .map(|k| self.rays[i][k])
                        .collect()
                })
                .collect();
            mult = gcd(mult, det(&minor).expect("the rays are bounded"));
        }
        let sigma: Vec<Vec<i64>> = [u1, wall[0], wall[1], wall[2]]
            .iter()
            .map(|&i| self.rays[i].to_vec())
            .collect();
        let d_u1 = Rational::from_integers(
            Integer::from(mult),
            Integer::from(det(&sigma).expect("the rays are bounded").abs()),
        );
        let mut d = vec![Rational::from(0); self.n()];
        d[u1] = d_u1.clone();
        for (&i, l) in others.iter().zip(lambda) {
            d[i] = l * &d_u1;
        }
        self.curve_class(&d)
    }

    /// The curve class $n$, $n_a = J_a \cdot C$, with intersections `d` with the divisors,
    /// from $D_i \cdot C = \sum_a Q_{ia} n_a$.
    fn curve_class(&self, d: &[Rational]) -> Vec<Rational> {
        let mut n = vec![Rational::from(0); self.h11()];
        for (p, t) in &self.inverse {
            if d[*p] != 0u32 {
                n.iter_mut().zip(t).for_each(|(x, y)| *x += &d[*p] * y);
            }
        }
        n
    }

    /// The distinct primitive classes of the walls, nonzero ones only, one per column.
    fn wall_classes(&self) -> Result<DMatrix<i32>, FanError> {
        let mut classes = Vec::new();
        for w in self.walls() {
            let c = primitive(&self.wall_class(w))?;
            if c.iter().any(|&x| x != 0) {
                classes.push(c);
            }
        }
        classes.sort_unstable();
        classes.dedup();
        let n = classes.len();
        Ok(DMatrix::from_iterator(
            self.h11(),
            n,
            classes.into_iter().flatten(),
        ))
    }

    /// The Mori cone, spanned by the classes of the walls, by its inequalities, and its default
    /// grading vector, the sum of its facets. The fan is checked to be regular: the grading
    /// vector is then positive on every wall class, and zero on the line of the Mori cone of a
    /// fan that is not. This is the only expensive part of the fan, a double description in
    /// dimension $h^{1,1}$.
    pub fn mori_cone(&self) -> Result<MoriCone, FanError> {
        let classes = self.wall_classes()?;
        let inequalities = cone::facets(&classes)?;
        let mut grading = RowDVector::zeros(self.h11());
        for h in columns(&inequalities.hyperplanes) {
            for (g, &x) in grading.iter_mut().zip(h) {
                *g = i32::checked_add(*g, x).ok_or(FanError::TooLarge)?;
            }
        }
        if inequalities.equations.ncols() > 0 || !positive(&grading, &classes) {
            return Err(FanError::NotRegular);
        }
        Ok(MoriCone {
            inequalities,
            grading_vector: grading,
        })
    }

    /// Checks that the grading vector is positive on the Mori cone, i.e. on every wall class.
    /// That shows the Mori cone pointed, i.e. the fan regular, without computing it. A fan
    /// that is not regular has no such vector, which the error does not tell apart: that would
    /// take the Mori cone, which takes minutes for a large fan.
    pub fn check_grading(&self, grading_vector: &RowDVector<i32>) -> Result<(), FanError> {
        if grading_vector.len() != self.h11() {
            return Err(FanError::NonpositiveGrading(None));
        }
        match columns(&self.wall_classes()?).find(|c| degree(grading_vector, c) <= 0) {
            Some(c) => Err(FanError::NonpositiveGrading(Some(c.to_vec()))),
            None => Ok(()),
        }
    }

    /// The first class of a wall, primitive, that is outside the cone of the inequalities
    /// `hyperplanes` ($h \cdot C \ge 0$) and `equations` ($e \cdot C = 0$), one per column,
    /// or `None` if the cone contains the Mori cone. The check of a Mori cone given with the
    /// fan, without computing the fan's: a cone that misses a wall misses curves.
    pub fn wall_outside(
        &self,
        hyperplanes: &DMatrix<i32>,
        equations: &DMatrix<i32>,
    ) -> Result<Option<Vec<i32>>, FanError> {
        let dot = |a: &[i32], b: &[i32]| -> i64 {
            a.iter().zip(b).map(|(&x, &y)| x as i64 * y as i64).sum()
        };
        let classes = self.wall_classes()?;
        let outside = columns(&classes)
            .find(|c| {
                columns(hyperplanes).any(|h| dot(h, c) < 0)
                    || columns(equations).any(|e| dot(e, c) != 0)
            })
            .map(<[i32]>::to_vec);
        Ok(outside)
    }

    /// The strata of the vex 3-cones, sorted by their divisors; none in a non-vex phase. The
    /// vex faces are checked to be 3-cones.
    pub fn strata(&self) -> Result<Vec<Stratum>, FanError> {
        let mut strata = Vec::new();
        for s in self.vex_cones() {
            let divisors: [usize; 3] = match s[..] {
                [a, b] => return Err(FanError::Vex2Cone([a, b])),
                [a, b, c] => [a, b, c],
                _ => unreachable!(),
            };
            let class = self.wall_class(divisors);
            if class.iter().any(|x| crate::fraction_i64(x).is_none()) {
                return Err(FanError::TooLarge);
            }
            strata.push(Stratum { divisors, class });
        }
        Ok(strata)
    }

    /// Checks the intersection numbers of a threefold hypersurface, $\kappa_{abc} = J_a J_b
    /// J_c (\sum_i D_i)$, in the basis of the charges and with each unordered triple given at
    /// most once (the missing ones are zero), against the fan's intersection ring.
    ///
    /// Both are compared as trilinear forms on $h^{1,1}$ independent divisors $D_i$, on which
    /// the given numbers are integers. The cost grows as $(h^{1,1})^4$.
    pub fn check_intnums(
        &self,
        intnums: &HashMap<(usize, usize, usize), i32>,
    ) -> Result<(), FanError> {
        let h = self.h11();
        let mut seen = HashMap::new();
        for (&(a, b, c), &v) in intnums {
            let mut k = [a, b, c];
            k.sort_unstable();
            if k[2] >= h || *seen.entry(k).or_insert(v) != v {
                return Err(FanError::InvalidIntnums);
            }
        }
        // Independent divisors (those of the inverse), and the numbers on them,
        // sum_abc q_ia q_jb q_kc kappa_abc, on sorted triples of positions in `basis`. Both
        // sides are sparse: kappa has few entries, and the fan's numbers vanish off its faces.
        let mut basis: Vec<usize> = self.inverse.iter().map(|(p, _)| *p).collect();
        basis.sort_unstable();
        let mut columns: Vec<Vec<(usize, i128)>> = vec![Vec::new(); h];
        for (i, &d) in basis.iter().enumerate() {
            for (a, &x) in self.q[d].iter().enumerate() {
                if x != 0 {
                    columns[a].push((i, x as i128));
                }
            }
        }
        let mut given: HashMap<[usize; 3], i128> = HashMap::new();
        for ([a, b, c], v) in seen {
            let mut orders = vec![
                [a, b, c],
                [a, c, b],
                [b, a, c],
                [b, c, a],
                [c, a, b],
                [c, b, a],
            ];
            orders.sort_unstable();
            orders.dedup();
            for [x, y, z] in orders {
                for &(i, qx) in &columns[x] {
                    for &(j, qy) in columns[y].iter().filter(|(j, _)| *j >= i) {
                        for &(k, qz) in columns[z].iter().filter(|(k, _)| *k >= j) {
                            let e = given.entry([i, j, k]).or_default();
                            *e = (v as i128)
                                .checked_mul(qx * qy)
                                .and_then(|p| p.checked_mul(qz))
                                .and_then(|p| e.checked_add(p))
                                .ok_or(FanError::TooLarge)?;
                        }
                    }
                }
            }
        }
        // The triples on the faces of the fan, with each face's divisors at least once.
        let position: HashMap<usize, usize> =
            basis.iter().enumerate().map(|(i, &d)| (d, i)).collect();
        let mut triples: Vec<[usize; 3]> = Vec::new();
        for f in &self.faces {
            let Some(p) = f
                .iter()
                .map(|d| position.get(d).copied())
                .collect::<Option<Vec<usize>>>()
            else {
                continue;
            };
            match p[..] {
                [i] => triples.push([i, i, i]),
                [i, j] => triples.extend([[i, i, j], [i, j, j]]),
                [i, j, k] => triples.push([i, j, k]),
                _ => {}
            }
        }
        triples.extend(given.keys().copied());
        triples.sort_unstable();
        triples.dedup();
        let mut ring = Ring {
            fan: self,
            memo: HashMap::new(),
        };
        for [i, j, k] in triples {
            let (a, b, c) = (basis[i], basis[j], basis[k]);
            let mut distinct = vec![a, b, c];
            distinct.dedup();
            let mut value = Rational::from(0);
            if self.faces.contains(&distinct) {
                // D_l only meets the face if it is a or shares a 2-cone with it.
                for &l in std::iter::once(&a).chain(&self.neighbors[a]) {
                    value += ring.number([a, b, c, l]);
                }
            }
            if value != given.get(&[i, j, k]).copied().unwrap_or(0) {
                return Err(FanError::IntnumsMismatch([a, b, c]));
            }
        }
        Ok(())
    }
}

/// The intersection numbers $D_{i_1} D_{i_2} D_{i_3} D_{i_4}$ of the toric variety, memoized.
struct Ring<'a> {
    fan: &'a Fan,
    memo: HashMap<[usize; 4], Rational>,
}

impl Ring<'_> {
    fn number(&mut self, mut m: [usize; 4]) -> Rational {
        m.sort_unstable();
        if let Some(v) = self.memo.get(&m) {
            return v.clone();
        }
        let fan = self.fan;
        let mut distinct = m.to_vec();
        distinct.dedup();
        if !fan.faces.contains(&distinct) {
            return Rational::from(0);
        }
        let value = if distinct.len() == 4 {
            let rows: Vec<Vec<i64>> = m.iter().map(|&i| fan.rays[i].to_vec()).collect();
            let det = det(&rows).expect("the rays are bounded");
            Rational::from_integers(Integer::from(1), Integer::from(det).abs())
        } else {
            // A repeated divisor D_i is replaced by a linearly equivalent combination of the
            // divisors outside the face: with <u, v_i> = 1 and <u, v_j> = 0 on the face's other
            // rays, 0 = sum_j <u, v_j> D_j gives D_i = -sum_{j outside} <u, v_j> D_j.
            let k = (0..3).find(|&k| m[k] == m[k + 1]).unwrap();
            let i = m[k];
            let face: Vec<usize> = std::iter::once(i)
                .chain(distinct.iter().copied().filter(|&j| j != i))
                .collect();
            let a: Vec<Vec<Rational>> = face
                .iter()
                .map(|&j| fan.rays[j].iter().map(|&x| Rational::from(x)).collect())
                .collect();
            let b: Vec<Rational> = (0..face.len())
                .map(|r| Rational::from((r == 0) as i64))
                .collect();
            let u = solve(&a, &b, 4).expect("the rays of a cone are independent");
            let mut rest = m;
            let mut value = Rational::from(0);
            // Only the divisors that share a 2-cone with D_i meet the rest.
            for &j in fan.neighbors[i].iter().filter(|j| !distinct.contains(j)) {
                let uj: Rational = u
                    .iter()
                    .zip(fan.rays[j])
                    .map(|(x, y)| x * Rational::from(y))
                    .sum();
                if uj != 0u32 {
                    rest[k] = j;
                    value -= uj * self.number(rest);
                }
            }
            value
        };
        self.memo.insert(m, value.clone());
        value
    }
}

/// A solution of $a x = b$ over the rationals, with $n$ unknowns, the free ones zero, or `None`
/// if there is none.
fn solve(a: &[Vec<Rational>], b: &[Rational], n: usize) -> Option<Vec<Rational>> {
    let mut m: Vec<Vec<Rational>> = a
        .iter()
        .zip(b)
        .map(|(r, x)| {
            r.iter()
                .cloned()
                .chain(std::iter::once(x.clone()))
                .collect()
        })
        .collect();
    let mut pivots = Vec::new();
    for c in 0..n {
        let r = pivots.len();
        let Some(p) = (r..m.len()).find(|&i| m[i][c] != 0u32) else {
            continue;
        };
        m.swap(r, p);
        let inv = Rational::from(1) / &m[r][c];
        m[r].iter_mut().for_each(|x| *x *= &inv);
        let pivot = m[r].clone();
        for (i, row) in m.iter_mut().enumerate() {
            if i != r && row[c] != 0u32 {
                let f = row[c].clone();
                row.iter_mut().zip(&pivot).for_each(|(x, y)| *x -= &f * y);
            }
        }
        pivots.push(c);
    }
    if m[pivots.len()..].iter().any(|r| r[n] != 0u32) {
        return None;
    }
    let mut x = vec![Rational::from(0); n];
    for (r, &c) in pivots.iter().enumerate() {
        x[c] = m[r][n].clone();
    }
    Some(x)
}

/// Pivot divisors $p_k$ and vectors $t_k$ with $n = \sum_k (Q n)_{p_k} t_k$ for every $n$: the
/// inverse of the rows $p_k$ of the charges `q` (one row per divisor), which have rank `h`.
/// The divisors with the fewest nonzero charges are tried first, so that an identity block,
/// which CYTools' charges have, takes no arithmetic.
fn left_inverse(q: &[Vec<i64>], h: usize) -> Vec<(usize, Vec<Rational>)> {
    // Row a is column a of the charges, then row a of the identity; reducing it to T Q^T | T
    // leaves t_k = (Q_P^{-1})^T e_k in row r_k, the row of the k-th pivot.
    let n = q.len();
    let mut m: Vec<Vec<Rational>> = (0..h)
        .map(|a| {
            (0..n + h)
                .map(|c| {
                    if c < n {
                        Rational::from(q[c][a])
                    } else {
                        Rational::from(i64::from(c - n == a))
                    }
                })
                .collect()
        })
        .collect();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&i| q[i].iter().filter(|&&x| x != 0).count());
    let mut pivots = Vec::with_capacity(h);
    let mut used = vec![false; h];
    for c in order {
        if pivots.len() == h {
            break;
        }
        let Some(r) = (0..h).find(|&r| !used[r] && m[r][c] != 0u32) else {
            continue;
        };
        used[r] = true;
        let inv = Rational::from(1) / &m[r][c];
        if inv != 1u32 {
            m[r].iter_mut().for_each(|x| *x *= &inv);
        }
        let pivot = m[r].clone();
        for (i, row) in m.iter_mut().enumerate() {
            if i != r && row[c] != 0u32 {
                let f = row[c].clone();
                row.iter_mut()
                    .zip(&pivot)
                    .filter(|(_, y)| **y != 0u32)
                    .for_each(|(x, y)| *x -= &f * y);
            }
        }
        pivots.push((c, r));
    }
    assert_eq!(pivots.len(), h, "the charges have rank h");
    pivots
        .into_iter()
        .map(|(c, r)| (c, m[r][n..].to_vec()))
        .collect()
}

/// Whether the grading vector is positive on every column of `classes`.
fn positive(grading_vector: &RowDVector<i32>, classes: &DMatrix<i32>) -> bool {
    columns(classes).all(|c| degree(grading_vector, c) > 0)
}

/// The degree of a class under the grading vector.
fn degree(grading_vector: &RowDVector<i32>, c: &[i32]) -> i64 {
    grading_vector
        .iter()
        .zip(c)
        .map(|(&w, &x)| w as i64 * x as i64)
        .sum()
}

fn gcd(a: i128, b: i128) -> i128 {
    if b == 0 {
        a.abs()
    } else {
        gcd(b, a % b)
    }
}

/// The primitive integer vector on the ray of a rational vector, or zero.
fn primitive(v: &[Rational]) -> Result<Vec<i32>, FanError> {
    let fractions: Vec<(i64, i64)> = v.iter().map(fraction).collect::<Result<_, _>>()?;
    let mut den = 1i128;
    for &(_, d) in &fractions {
        den = (den / gcd(den, d as i128))
            .checked_mul(d as i128)
            .ok_or(FanError::TooLarge)?;
    }
    let ints: Vec<i128> = fractions
        .iter()
        .map(|&(n, d)| {
            (n as i128)
                .checked_mul(den / d as i128)
                .ok_or(FanError::TooLarge)
        })
        .collect::<Result<_, _>>()?;
    let g = ints.iter().fold(0, |g, &x| gcd(g, x)).max(1);
    ints.iter()
        .map(|&x| i32::try_from(x / g).map_err(|_| FanError::TooLarge))
        .collect()
}

/// A rational as a fraction (numerator, positive denominator) of 64-bit integers.
fn fraction(x: &Rational) -> Result<(i64, i64), FanError> {
    let num = i64::try_from(&Integer::from_sign_and_abs_ref(
        *x >= 0u32,
        x.numerator_ref(),
    ));
    let den = i64::try_from(x.denominator_ref());
    num.ok().zip(den.ok()).ok_or(FanError::TooLarge)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The fan with maximal cones `cones` and the charges with the given rows, one per Kähler
    /// parameter.
    fn fan(cones: &[&[usize]], q: &[&[i32]]) -> Result<Fan, FanError> {
        let cones: Vec<Vec<usize>> = cones.iter().map(|c| c.to_vec()).collect();
        let q = DMatrix::from_fn(q[0].len(), q.len(), |i, a| q[a][i]);
        Fan::new(&cones, &q)
    }

    /// The cones of the fan of $\mathbb{P}^4$, or of a weighted one: every four of five rays.
    const SIMPLEX: &[&[usize]] = &[
        &[0, 1, 2, 3],
        &[0, 1, 2, 4],
        &[0, 1, 3, 4],
        &[0, 2, 3, 4],
        &[1, 2, 3, 4],
    ];

    #[test]
    fn test_p4() {
        let fan = fan(SIMPLEX, &[&[1, 1, 1, 1, 1]]).unwrap();
        // The rays sum to zero and span the lattice.
        let rays = fan.rays();
        assert_eq!(rays.column_sum(), nalgebra::DVector::zeros(4));
        assert_eq!(
            det(&(0..4)
                .map(|i| rays.column(i).iter().map(|&x| x as i64).collect())
                .collect::<Vec<_>>())
            .unwrap()
            .abs(),
            1
        );
        assert_eq!(fan.walls().len(), 10);
        assert_eq!(fan.wall_class([0, 1, 2]), vec![Rational::from(1)]);
        assert!(fan.vex_cones().is_empty());
        let mori = fan.mori_cone().unwrap();
        assert_eq!(
            mori.inequalities.hyperplanes,
            DMatrix::from_element(1, 1, 1)
        );
        assert_eq!(mori.grading_vector, RowDVector::from_element(1, 1));
        assert_eq!(fan.strata(), Ok(Vec::new()));
        assert_eq!(fan.check_grading(&RowDVector::from_element(1, 2)), Ok(()));
        assert!(matches!(
            fan.check_grading(&RowDVector::from_element(1, -1)),
            Err(FanError::NonpositiveGrading(Some(_)))
        ));
        assert_eq!(
            fan.check_grading(&RowDVector::zeros(2)),
            Err(FanError::NonpositiveGrading(None))
        );
        // The quintic.
        let mut intnums = HashMap::from([((0, 0, 0), 5)]);
        assert_eq!(fan.check_intnums(&intnums), Ok(()));
        intnums.insert((0, 0, 0), 4);
        assert_eq!(
            fan.check_intnums(&intnums),
            Err(FanError::IntnumsMismatch([0, 0, 0]))
        );
        intnums.insert((0, 0, 1), 4);
        assert_eq!(fan.check_intnums(&intnums), Err(FanError::InvalidIntnums));
    }

    #[test]
    fn test_invalid_fans() {
        let p4: &[&[i32]] = &[&[1, 1, 1, 1, 1]];
        assert_eq!(
            fan(SIMPLEX, &[&[1, 1, 1, 1]]).unwrap_err(),
            FanError::NotFourfold
        );
        assert_eq!(
            fan(SIMPLEX, &[&[1, 1, 1, 1, 1], &[1, 1, 1, 1, 1]]).unwrap_err(),
            FanError::NotFourfold
        );
        assert_eq!(fan(&SIMPLEX[1..], p4).unwrap_err(), FanError::NotComplete);
        for cone in [&[0, 1, 2][..], &[0, 1, 2, 2], &[0, 1, 2, 5]] {
            let cones = [SIMPLEX, &[cone]].concat();
            assert_eq!(
                fan(&cones, p4).unwrap_err(),
                FanError::InvalidCone(cone.to_vec())
            );
        }
        let cones = [SIMPLEX, &[&[3, 2, 1, 0][..]]].concat();
        assert_eq!(
            fan(&cones, p4).unwrap_err(),
            FanError::InvalidCone(vec![3, 2, 1, 0])
        );
        // P(1,1,1,1,3), which is not Gorenstein.
        assert_eq!(
            fan(SIMPLEX, &[&[1, 1, 1, 1, 3]]).unwrap_err(),
            FanError::NotReflexive
        );
        // A fan of P(1,1,1,1,4) refined by the point interior to the facet opposite the ray of
        // weight 4.
        let cones: &[&[usize]] = &[
            &[0, 1, 2, 3],
            &[0, 1, 2, 4],
            &[0, 1, 3, 4],
            &[0, 2, 3, 5],
            &[0, 2, 4, 5],
            &[0, 3, 4, 5],
            &[1, 2, 3, 4],
            &[2, 3, 4, 5],
        ];
        assert_eq!(
            fan(cones, &[&[1, 0, 1, 1, 1, -4], &[0, 1, 0, 0, 0, 1]]).unwrap_err(),
            FanError::InteriorToFacet(5)
        );
    }
}
