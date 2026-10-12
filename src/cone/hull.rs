//! Conversion between the two descriptions of a polyhedral cone, by spanning vectors and by
//! inequalities, with the double description method in exact integer arithmetic.
//!
//! Both directions compute the dual cone $\{a : a \cdot v \ge 0 \text{ for all } v\}$ of a
//! finite set of vectors. Its extreme rays and lineality space are the facets and equations of
//! the cone the vectors span ([`facets`]), and the dual of the inequalities and of both signs
//! of the equations of a cone is the cone itself ([`rays`]). The dual is built by Motzkin's
//! iteration as Fukuda and Prodon formulate it (and cdd implements it): start from the whole
//! space, all of it lineality, and add the inequalities $a \cdot v \ge 0$ one at a time. While
//! the lineality space has a vector $l$ with $l \cdot v \neq 0$, $l$ becomes a ray and the rest
//! is made orthogonal to $v$ with it. Otherwise every pair of adjacent rays on opposite sides
//! of the hyperplane gives a new ray on it, and the rays on the negative side go. Adjacency is
//! tested combinatorially: two rays are adjacent if no other ray is tight on every inequality
//! both are tight on.
//!
//! The order of the inequalities only changes the run time, but by a lot: the iteration takes
//! the short vectors first, which kept the intermediate cones small on every Mori cone tried
//! (`PLAN.md`, phase 7a).

use super::columns;
use nalgebra::DMatrix;
use std::fmt;

/// A cone by inequalities: $\{C : h \cdot C \ge 0,\ e \cdot C = 0\}$ for the columns $h$ of
/// `hyperplanes` and $e$ of `equations`.
///
/// As [`facets`] returns it, the description is irredundant and canonical: the hyperplanes are
/// the facets, and the equations a basis of the orthogonal complement of the cone's span. The
/// equations are in reduced row echelon form scaled to primitive integer vectors, and each
/// hyperplane is the one primitive integer vector that defines its facet and vanishes at
/// the pivot coordinates of the equations; both are sorted lexicographically.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Inequalities {
    /// The inequalities $h \cdot C \ge 0$, one per column.
    pub hyperplanes: DMatrix<i32>,
    /// The equations $e \cdot C = 0$, one per column.
    pub equations: DMatrix<i32>,
}

/// A cone by generators: the non-negative combinations of the columns of `rays` plus the
/// linear combinations of those of `lineality`.
///
/// As [`rays`] returns it, the description is irredundant and canonical in the same way as
/// [`Inequalities`]: `lineality` is a basis of the largest linear subspace in the cone, and
/// `rays` holds one primitive integer vector for each of its extreme rays modulo that subspace.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Generators {
    /// The extreme rays (modulo the lineality space), one per column.
    pub rays: DMatrix<i32>,
    /// A basis of the lineality space, one vector per column.
    pub lineality: DMatrix<i32>,
}

/// The result, or an intermediate value of the computation, overflowed: an entry of the result
/// does not fit 32 bits, or one of an intermediate ray 128.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TooLarge;

impl fmt::Display for TooLarge {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "the entries of the cone's description are too large")
    }
}

impl std::error::Error for TooLarge {}

/// The facets and equations of the cone spanned by the columns of `generators`, in any
/// dimension. Redundant, repeated and zero columns are allowed; no columns give the cone
/// $\{0\}$, which has no facets.
pub fn facets(generators: &DMatrix<i32>) -> Result<Inequalities, TooLarge> {
    let (rays, lineality) = dual(generators.nrows(), columns(generators).map(widen).collect())?;
    Ok(Inequalities {
        hyperplanes: to_matrix(generators.nrows(), rays)?,
        equations: to_matrix(generators.nrows(), lineality)?,
    })
}

/// The extreme rays and the lineality space of the cone cut out by the columns of
/// `hyperplanes` and `equations`, which must have the same number of rows (or no columns).
pub fn rays(hyperplanes: &DMatrix<i32>, equations: &DMatrix<i32>) -> Result<Generators, TooLarge> {
    let n = if hyperplanes.ncols() > 0 {
        hyperplanes.nrows()
    } else {
        equations.nrows()
    };
    assert!(equations.ncols() == 0 || equations.nrows() == n);
    let mut vectors: Vectors = columns(hyperplanes).map(widen).collect();
    for e in columns(equations) {
        vectors.push(widen(e));
        vectors.push(e.iter().map(|&x| -(x as i128)).collect());
    }
    let (rays, lineality) = dual(n, vectors)?;
    Ok(Generators {
        rays: to_matrix(n, rays)?,
        lineality: to_matrix(n, lineality)?,
    })
}

/// Integer vectors, all of the same length.
type Vectors = Vec<Vec<i128>>;

fn widen(v: &[i32]) -> Vec<i128> {
    v.iter().map(|&x| x as i128).collect()
}

/// The vectors as the columns of an `n`-row matrix.
fn to_matrix(n: usize, vectors: Vectors) -> Result<DMatrix<i32>, TooLarge> {
    let ncols = vectors.len();
    let entries = vectors
        .into_iter()
        .flatten()
        .map(|x| i32::try_from(x).map_err(|_| TooLarge))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(DMatrix::from_vec(n, ncols, entries))
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

/// Divides `v` by the gcd of its entries.
fn make_primitive(v: &mut [i128]) {
    let g = v.iter().fold(0, |g, x| gcd(g, x.unsigned_abs()));
    if g > 1 {
        // g is at most the largest absolute value, which fits since no value is i128::MIN.
        let g = g as i128;
        v.iter_mut().for_each(|x| *x /= g);
    }
}

/// Treats `i128::MIN` as an overflow, so that every value can be negated and its absolute value
/// taken. The inputs, from `i32`, are never it, and the other values come from here.
fn checked(x: Option<i128>) -> Result<i128, TooLarge> {
    x.filter(|&x| x != i128::MIN).ok_or(TooLarge)
}

fn dot(a: &[i128], b: &[i128]) -> Result<i128, TooLarge> {
    a.iter().zip(b).try_fold(0i128, |s, (&x, &y)| {
        checked(x.checked_mul(y).and_then(|p| s.checked_add(p)))
    })
}

/// The primitive vector in the direction of $s a - t b$.
fn combine(s: i128, a: &[i128], t: i128, b: &[i128]) -> Result<Vec<i128>, TooLarge> {
    let g = gcd(s.unsigned_abs(), t.unsigned_abs());
    let g = i128::try_from(g).map_err(|_| TooLarge)?;
    let (s, t) = (s / g, t / g);
    let mut v = a
        .iter()
        .zip(b)
        .map(|(&x, &y)| {
            let (sx, ty) = (s.checked_mul(x), t.checked_mul(y));
            checked(sx.zip(ty).and_then(|(sx, ty)| sx.checked_sub(ty)))
        })
        .collect::<Result<Vec<_>, _>>()?;
    make_primitive(&mut v);
    Ok(v)
}

/// A ray of the dual cone under construction, with the set of inequalities it is tight on, as a
/// bit set over their indices.
struct Ray {
    v: Vec<i128>,
    tight: Vec<u64>,
}

fn insert(set: &mut [u64], i: usize) {
    set[i / 64] |= 1 << (i % 64);
}

/// The extreme rays and a basis of the lineality space of the dual of the cone spanned by the
/// `n`-dimensional `vectors`, both canonical as [`Inequalities`] describes.
fn dual(n: usize, mut vectors: Vectors) -> Result<(Vectors, Vectors), TooLarge> {
    vectors.retain(|v| v.iter().any(|&x| x != 0));
    vectors.sort_by_cached_key(|v| (v.iter().map(|x| x.unsigned_abs()).sum::<u128>(), v.clone()));
    vectors.dedup();
    let words = vectors.len().div_ceil(64);

    let mut lineality: Vectors = (0..n)
        .map(|i| (0..n).map(|j| (i == j) as i128).collect())
        .collect();
    let mut rays: Vec<Ray> = Vec::new();
    for (k, v) in vectors.iter().enumerate() {
        let mut first = None;
        for (i, l) in lineality.iter().enumerate() {
            let s = dot(l, v)?;
            if s != 0 {
                first = Some((i, s));
                break;
            }
        }
        if let Some((p, mut s)) = first {
            // The lineality vector becomes a ray, on the positive side, and the rest of the
            // lineality space and the rays are moved onto the hyperplane along it. That keeps
            // every earlier inequality tight where it was, since the vector is orthogonal to
            // all of them.
            let mut l = lineality.swap_remove(p);
            if s < 0 {
                l.iter_mut().for_each(|x| *x = -*x);
                s = -s;
            }
            for m in lineality.iter_mut() {
                let t = dot(m, v)?;
                if t != 0 {
                    *m = combine(s, m, t, &l)?;
                }
            }
            for r in rays.iter_mut() {
                let t = dot(&r.v, v)?;
                if t != 0 {
                    r.v = combine(s, &r.v, t, &l)?;
                }
                insert(&mut r.tight, k);
            }
            let mut tight = vec![0; words];
            (0..k).for_each(|j| insert(&mut tight, j));
            rays.push(Ray { v: l, tight });
            continue;
        }

        let sides = rays
            .iter()
            .map(|r| dot(&r.v, v))
            .collect::<Result<Vec<_>, _>>()?;
        let positive: Vec<usize> = (0..rays.len()).filter(|&i| sides[i] > 0).collect();
        let negative: Vec<usize> = (0..rays.len()).filter(|&i| sides[i] < 0).collect();
        let mut new = Vec::new();
        for &i in &positive {
            for &j in &negative {
                let common: Vec<u64> = (rays[i].tight.iter().zip(&rays[j].tight))
                    .map(|(a, b)| a & b)
                    .collect();
                let adjacent = (0..rays.len()).all(|m| {
                    m == i || m == j || common.iter().zip(&rays[m].tight).any(|(c, t)| c & !t != 0)
                });
                if adjacent {
                    // sides[i] > 0 > sides[j], so this is a positive combination.
                    let v = combine(sides[i], &rays[j].v, sides[j], &rays[i].v)?;
                    let mut tight = common;
                    insert(&mut tight, k);
                    new.push(Ray { v, tight });
                }
            }
        }
        rays = rays
            .into_iter()
            .zip(&sides)
            .filter(|(_, &s)| s >= 0)
            .map(|(mut r, &s)| {
                if s == 0 {
                    insert(&mut r.tight, k);
                }
                r
            })
            .chain(new)
            .collect();
    }

    let lineality = echelon(lineality)?;
    let mut rays = rays
        .into_iter()
        .map(|r| reduce(r.v, &lineality))
        .collect::<Result<Vec<_>, _>>()?;
    rays.sort();
    Ok((rays, lineality))
}

/// The pivot column of a row of an echelon form, its first nonzero entry.
fn pivot(row: &[i128]) -> usize {
    row.iter().position(|&x| x != 0).unwrap()
}

/// A basis of the span of the linearly independent `rows` in reduced row echelon form, each
/// row scaled to a primitive integer vector with a positive pivot, sorted lexicographically.
/// The reduced row echelon form is unique, so the result only depends on the span.
fn echelon(mut rows: Vectors) -> Result<Vectors, TooLarge> {
    let n = rows.first().map_or(0, Vec::len);
    let mut r = 0;
    for c in 0..n {
        let Some(p) = (r..rows.len()).find(|&i| rows[i][c] != 0) else {
            continue;
        };
        rows.swap(r, p);
        if rows[r][c] < 0 {
            rows[r].iter_mut().for_each(|x| *x = -*x);
        }
        make_primitive(&mut rows[r]);
        for i in 0..rows.len() {
            if i != r && rows[i][c] != 0 {
                let (s, t) = (rows[r][c], rows[i][c]);
                rows[i] = combine(s, &rows[i], t, &rows[r])?;
            }
        }
        r += 1;
    }
    rows.truncate(r);
    rows.sort();
    Ok(rows)
}

/// `v` moved along the rows of `echelon` to vanish at their pivot columns, which leaves its
/// inner products with the vectors orthogonal to them unchanged and makes it the unique such
/// representative, then made primitive.
fn reduce(mut v: Vec<i128>, echelon: &[Vec<i128>]) -> Result<Vec<i128>, TooLarge> {
    for e in echelon {
        let p = pivot(e);
        if v[p] != 0 {
            // e[p] > 0 keeps v on its side of every hyperplane.
            v = combine(e[p], &v, v[p], e)?;
        }
    }
    make_primitive(&mut v);
    Ok(v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::dmatrix;

    /// The columns of a matrix as vectors.
    fn cols(m: &DMatrix<i32>) -> Vec<Vec<i32>> {
        columns(m).map(<[i32]>::to_vec).collect()
    }

    fn sorted(mut v: Vec<Vec<i32>>) -> Vec<Vec<i32>> {
        v.sort();
        v
    }

    fn det(m: &[Vec<i64>]) -> i64 {
        match m.len() {
            0 => 1,
            n => (0..n)
                .map(|j| {
                    let minor: Vec<Vec<i64>> = m[1..]
                        .iter()
                        .map(|r| [&r[..j], &r[j + 1..]].concat())
                        .collect();
                    let sign = if j % 2 == 0 { 1 } else { -1 };
                    sign * m[0][j] * det(&minor)
                })
                .sum(),
        }
    }

    /// The facets of a full-dimensional cone by brute force: the normal of every $n - 1$ of
    /// its generators (the vector of signed maximal minors) that the generators all lie on
    /// one side of, and are not all tight on.
    fn brute_force_facets(gens: &[Vec<i32>]) -> Vec<Vec<i32>> {
        let n = gens[0].len();
        let mut facets = Vec::new();
        let mut subset = vec![0usize; n - 1];
        fn rec(
            gens: &[Vec<i32>],
            start: usize,
            depth: usize,
            subset: &mut Vec<usize>,
            facets: &mut Vec<Vec<i32>>,
        ) {
            let n = gens[0].len();
            if depth == n - 1 {
                let normal: Vec<i64> = (0..n)
                    .map(|j| {
                        let m: Vec<Vec<i64>> = subset
                            .iter()
                            .map(|&s| {
                                (0..n)
                                    .filter(|&c| c != j)
                                    .map(|c| gens[s][c] as i64)
                                    .collect()
                            })
                            .collect();
                        let sign = if (n - 1 + j).is_multiple_of(2) { 1 } else { -1 };
                        sign * det(&m)
                    })
                    .collect();
                for sign in [1, -1] {
                    let mut h: Vec<i128> = normal.iter().map(|&x| (sign * x) as i128).collect();
                    let d: Vec<i128> = gens.iter().map(|g| dot(&h, &widen(g)).unwrap()).collect();
                    if h.iter().any(|&x| x != 0)
                        && d.iter().all(|&x| x >= 0)
                        && d.iter().any(|&x| x > 0)
                    {
                        make_primitive(&mut h);
                        facets.push(h.into_iter().map(|x| x as i32).collect());
                    }
                }
                return;
            }
            for s in start..gens.len() {
                subset[depth] = s;
                rec(gens, s + 1, depth + 1, subset, facets);
            }
        }
        rec(gens, 0, 0, &mut subset, &mut facets);
        facets.sort();
        facets.dedup();
        facets
    }

    #[test]
    fn test_simple_cones() {
        // The positive orthant, with a redundant, a repeated and a zero generator.
        let f = facets(&dmatrix![1, 0, 0, 1, 1, 0; 0, 1, 0, 1, 1, 0; 0, 0, 1, 1, 1, 0]).unwrap();
        assert_eq!(f.hyperplanes, dmatrix![0, 0, 1; 0, 1, 0; 1, 0, 0]);
        assert_eq!(f.equations.shape(), (3, 0));

        // The cone over a square has four facets.
        let f = facets(&dmatrix![1, 1, 1, 1; 1, -1, 1, -1; 1, 1, -1, -1]).unwrap();
        assert_eq!(
            cols(&f.hyperplanes),
            [[1, -1, 0], [1, 0, -1], [1, 0, 1], [1, 1, 0]]
        );

        // A two-dimensional cone in three dimensions: its span is x = y, and each facet is
        // represented by the vector vanishing at the pivot of the equation, x.
        let f = facets(&dmatrix![1, 1; 1, 1; 0, 2]).unwrap();
        assert_eq!(cols(&f.equations), [[1, -1, 0]]);
        assert_eq!(cols(&f.hyperplanes), [[0, 0, 1], [0, 2, -1]]);

        // A ray, a half-plane (with a line), a whole plane, and the cone {0}.
        let f = facets(&dmatrix![2; 0; 4]).unwrap();
        assert_eq!(cols(&f.hyperplanes), [[0, 0, 1]]);
        assert_eq!(cols(&f.equations), [[0, 1, 0], [2, 0, -1]]);
        let f = facets(&dmatrix![1, -1, 0; 0, 0, 1]).unwrap();
        assert_eq!(cols(&f.hyperplanes), [[0, 1]]);
        assert_eq!(f.equations.shape(), (2, 0));
        let f = facets(&dmatrix![1, -1, 0, 0; 0, 0, 1, -1]).unwrap();
        assert_eq!(f.hyperplanes.shape(), (2, 0));
        assert_eq!(f.equations.shape(), (2, 0));
        let f = facets(&DMatrix::zeros(2, 0)).unwrap();
        assert_eq!(f.hyperplanes.shape(), (2, 0));
        assert_eq!(cols(&f.equations), [[0, 1], [1, 0]]);
    }

    #[test]
    fn test_rays() {
        // The positive orthant, with a redundant inequality.
        let g = rays(&dmatrix![1, 0, 1; 0, 1, 1], &DMatrix::zeros(0, 0)).unwrap();
        assert_eq!(g.rays, dmatrix![0, 1; 1, 0]);
        assert_eq!(g.lineality.shape(), (2, 0));

        // The cone over a square, back from its facets, and its faces on the planes z = 0,
        // y = 0 and x = 0, the last of which is {0}.
        let square = dmatrix![1, 1, 1, 1; 1, -1, 1, -1; 1, 1, -1, -1];
        let f = facets(&square).unwrap();
        let g = rays(&f.hyperplanes, &f.equations).unwrap();
        assert_eq!(sorted(cols(&g.rays)), sorted(cols(&square)));
        let g = rays(&f.hyperplanes, &dmatrix![0; 0; 1]).unwrap();
        assert_eq!(cols(&g.rays), [[1, -1, 0], [1, 1, 0]]);
        assert_eq!(g.lineality.shape(), (3, 0));
        let g = rays(&f.hyperplanes, &dmatrix![0; 1; 0]).unwrap();
        assert_eq!(cols(&g.rays), [[1, 0, -1], [1, 0, 1]]);
        let g = rays(&f.hyperplanes, &dmatrix![1; 0; 0]).unwrap();
        assert_eq!(g.rays.shape(), (3, 0));
        assert_eq!(g.lineality.shape(), (3, 0));

        // A half-plane, with a line, and the whole space.
        let g = rays(&dmatrix![0; 1], &DMatrix::zeros(2, 0)).unwrap();
        assert_eq!(cols(&g.rays), [[0, 1]]);
        assert_eq!(cols(&g.lineality), [[1, 0]]);
        let g = rays(&DMatrix::zeros(3, 0), &DMatrix::zeros(3, 0)).unwrap();
        assert_eq!(g.rays.shape(), (3, 0));
        assert_eq!(cols(&g.lineality), [[0, 0, 1], [0, 1, 0], [1, 0, 0]]);
    }

    /// Random cones against the brute force, and back. A third of them lie on the hyperplane
    /// $x_n = x_1$; their facets are those of their projection that drops $x_n$, which is
    /// full-dimensional, lifted to vanish at $x_1$, the pivot of the equation $x_1 - x_n = 0$.
    #[test]
    fn test_random_cones() {
        let mut state = 0x2545f4914f6cdd1du64;
        let mut next = |m: u64| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % m
        };
        for n in 2..=5 {
            for _ in 0..40 {
                let count = n + next(6) as usize;
                let on_plane = n > 2 && next(3) == 0;
                let gens: Vec<Vec<i32>> = (0..count)
                    .map(|_| {
                        let mut g: Vec<i32> = (0..n).map(|_| next(7) as i32 - 3).collect();
                        if on_plane {
                            g[n - 1] = g[0];
                        }
                        g
                    })
                    .collect();
                let m = DMatrix::from_vec(n, count, gens.concat());
                let f = facets(&m).unwrap();
                if on_plane {
                    let projected: Vec<Vec<i32>> =
                        gens.iter().map(|g| g[..n - 1].to_vec()).collect();
                    if facets(&DMatrix::from_vec(n - 1, count, projected.concat()))
                        .unwrap()
                        .equations
                        .ncols()
                        > 0
                    {
                        continue;
                    }
                    let mut e = vec![0; n];
                    (e[0], e[n - 1]) = (1, -1);
                    assert_eq!(cols(&f.equations), [e], "{m}");
                    let lifted: Vec<Vec<i32>> = brute_force_facets(&projected)
                        .into_iter()
                        .map(|h| [&[0], &h[1..], &h[..1]].concat())
                        .collect();
                    assert_eq!(cols(&f.hyperplanes), sorted(lifted), "{m}");
                } else if f.equations.ncols() == 0 {
                    assert_eq!(cols(&f.hyperplanes), brute_force_facets(&gens), "{m}");
                }
                // The rays they give back span the same cone.
                let g = rays(&f.hyperplanes, &f.equations).unwrap();
                let mut back = cols(&g.rays);
                for l in cols(&g.lineality) {
                    back.push(l.iter().map(|x| -x).collect());
                    back.push(l);
                }
                let back = DMatrix::from_vec(n, back.len(), back.concat());
                assert_eq!(facets(&back).unwrap(), f, "{m}");
            }
        }
    }

    #[test]
    fn test_too_large() {
        let big = i32::MAX;
        assert!(facets(&dmatrix![big, 1; 1, big]).is_ok());
        // The facet through (big, 1, 0) and (0, big, 1) has the normal (1, -big, big^2).
        assert_eq!(
            facets(&dmatrix![big, 0, 1; 1, big, 0; 0, 1, big]),
            Err(TooLarge)
        );
    }
}
