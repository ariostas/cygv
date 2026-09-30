//! A module for properties of polynomials.

use crate::polynomial::coefficient::PolynomialCoeff;
use crate::semigroup::Semigroup;
use rustc_hash::FxHashMap as HashMap;

/// The fingerprint of a monomial. See [`PolynomialProperties::fingerprint`].
///
/// Sixty-four bits are enough because no lookup relies on distinct monomials
/// having distinct fingerprints by chance: the elements are checked to have
/// distinct ones when they are built, and a product either cannot be a
/// non-element or is checked against the monomial it stands for. The width only
/// sets how often that first check has to draw a new set of coefficients.
pub type Fingerprint = u64;

/// A structure containing properties of a polynomial.
///
/// This structure contains the data common to all polynomials and is needed
/// for most operations.
#[derive(Clone, Debug)]
pub struct PolynomialProperties<'a, T> {
    pub semigroup: &'a Semigroup,
    /// The fingerprint of each element of the semigroup, indexed the same way
    /// its columns are.
    ///
    /// Since the fingerprint is linear, the one of a product of two monomials
    /// is the sum of theirs, so multiplication never has to look at the
    /// monomials themselves.
    pub fingerprints: Vec<Fingerprint>,
    /// The coefficients of the fingerprint, one per coordinate.
    coefficients: Vec<Fingerprint>,
    /// Map from the fingerprint of a monomial to its index.
    monomial_map: HashMap<Fingerprint, u32>,
    /// Whether a product found by its fingerprint has to be checked against the
    /// monomial it is meant to be.
    ///
    /// False when the sums of elements were shown to stay in the semigroup, in
    /// which case a product is always an element and its fingerprint always
    /// finds it. True otherwise, which costs a pass over the coordinates for
    /// each product, and is what a pruned set of elements such as
    /// [`Semigroup::with_target_points`] builds pays.
    pub verify_products: bool,
    pub zero_cutoff: T,
    /// A zero coefficient, kept around as the template new coefficients are
    /// cloned from. For [`rug::Float`] coefficients it carries the precision
    /// they all have to be created with.
    pub zero: T,
}

impl<'a, T: PolynomialCoeff<T>> PolynomialProperties<'a, T> {
    /// Create a new PolynomialProperties structure.
    pub fn new(semigroup: &'a Semigroup, zero_cutoff: &T) -> Self {
        let mut zero = zero_cutoff.clone();
        zero.assign(0i32);

        let n_coords = semigroup.elements.nrows();
        assert!(
            semigroup.elements.ncols() <= u32::MAX as usize,
            "the semigroup has more elements than an index can hold"
        );
        // The elements must have distinct fingerprints, or the map would lose
        // one of them and a lookup could answer with the wrong element; that is
        // checked here rather than assumed, and the coefficients are redrawn
        // until it holds.
        let (coefficients, fingerprints, monomial_map) = (0..)
            .find_map(|attempt| {
                let coefficients = fingerprint_coefficients(n_coords, attempt);
                let fingerprints: Vec<Fingerprint> = semigroup
                    .elements
                    .column_iter()
                    .map(|c| fingerprint_with(&coefficients, c.as_slice()))
                    .collect();
                let monomial_map: HashMap<Fingerprint, u32> = fingerprints
                    .iter()
                    .enumerate()
                    .map(|(i, &f)| (f, i as u32))
                    .collect();
                (monomial_map.len() == fingerprints.len()).then_some((
                    coefficients,
                    fingerprints,
                    monomial_map,
                ))
            })
            .expect("the search for a collision-free fingerprint never gives up");

        // Whether a product can be trusted to the fingerprint alone is settled
        // once, here, rather than reasoned about at each lookup.
        let verify_products =
            !sums_stay_in_the_semigroup(semigroup, &coefficients, &fingerprints, &monomial_map);

        Self {
            semigroup,
            fingerprints,
            coefficients,
            monomial_map,
            verify_products,
            zero_cutoff: zero_cutoff.clone(),
            zero,
        }
    }

    /// The fingerprint of a monomial.
    ///
    /// The fingerprint is the dot product of the monomial with a fixed vector
    /// of random coefficients, computed modulo $2^{64}$. It is therefore
    /// linear: the fingerprint of a product of monomials is the (wrapping) sum
    /// of their fingerprints, and the one of a ratio is their difference. That
    /// is what lets [`Polynomial::mul`](crate::polynomial::Polynomial::mul)
    /// look a product up in constant time, instead of assembling the product
    /// monomial and hashing it, which costs a pass over all $h^{1,1}$
    /// coordinates.
    #[inline]
    pub fn fingerprint(&self, monomial: &[i32]) -> Fingerprint {
        fingerprint_with(&self.coefficients, monomial)
    }

    /// The index of the semigroup element with a given fingerprint, if there is
    /// one.
    ///
    /// The elements are checked to have distinct fingerprints when they are
    /// built, so an element always finds its own index. A monomial that is not
    /// an element can share a fingerprint with one, and is then answered with
    /// that element, so every caller either cannot ask about a non-element (the
    /// products, when [`PolynomialProperties::verify_products`] is false) or
    /// checks the answer against the monomial itself.
    #[inline]
    pub fn index_of_fingerprint(&self, fingerprint: Fingerprint) -> Option<usize> {
        self.monomial_map.get(&fingerprint).map(|&i| i as usize)
    }

    /// The index of a monomial in the semigroup, if it is one of its elements.
    #[inline]
    pub fn monomial_index(&self, monomial: &[i32]) -> Option<usize> {
        self.index_of_fingerprint(self.fingerprint(monomial))
    }
}

/// The fingerprint of a monomial under the given coefficients.
#[inline]
fn fingerprint_with(coefficients: &[Fingerprint], monomial: &[i32]) -> Fingerprint {
    monomial.iter().zip(coefficients).fold(0, |acc, (&x, &c)| {
        acc.wrapping_add(c.wrapping_mul(x as i64 as Fingerprint))
    })
}

/// Draw the coefficients of the fingerprint.
///
/// They are odd so that a coordinate cannot be lost to the low bits of the
/// product, and they come from a fixed sequence rather than from the system
/// entropy so that a computation gives the same answer every time it is run.
fn fingerprint_coefficients(n_coords: usize, attempt: u64) -> Vec<Fingerprint> {
    // SplitMix64, which is enough to spread a counter into unrelated words.
    let mut state =
        0x243f_6a88_85a3_08d3_u64.wrapping_add(attempt.wrapping_mul(0x9e37_79b9_7f4a_7c15));
    let mut next = || {
        state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    };
    (0..n_coords).map(|_| next() | 1).collect()
}

/// Check that the sum of any two elements is an element as well, whenever its
/// degree is within the semigroup's maximum.
///
/// That is the property that makes a product looked up by fingerprint exact:
/// [`Polynomial::mul`](crate::polynomial::Polynomial::mul) skips the products
/// whose degree is too large, so if this holds it only ever looks up elements,
/// and the elements have distinct fingerprints.
///
/// Checking it over every pair of elements would cost a pass over
/// `n_elements^2` pairs. It is enough to check it against the generators, which
/// is a pass over `n_elements * n_generators`:
///
/// - every nonzero element `e` has a generator `g` with `e - g` an element, so
///   by induction on the degree every element is a sum of generators;
/// - for every element `e` and generator `g` with `deg(e) + deg(g)` within the
///   maximum, `e + g` is an element.
///
/// Given both, write `b` as `g_1 + ... + g_m`. Adding those to an `a` one at a
/// time stays within the maximum degree as long as `deg(a) + deg(b)` does, so
/// each partial sum is an element by the second point, and `a + b` is the last
/// of them.
///
/// The cost is a pass over `n_elements * n_generators` pairs, which is why it
/// is done for every semigroup rather than only where it is expected to hold.
/// It returns false if the semigroup records no generators, or if either point
/// fails, as a pruned set of elements makes it; the products are then checked
/// one by one as they are looked up, which is bounded by the work the
/// multiplication does anyway, unlike a pass over every pair of elements.
fn sums_stay_in_the_semigroup(
    semigroup: &Semigroup,
    coefficients: &[Fingerprint],
    fingerprints: &[Fingerprint],
    monomial_map: &HashMap<Fingerprint, u32>,
) -> bool {
    let Some(generators) = &semigroup.generators else {
        return false;
    };
    let elements = &semigroup.elements;
    let grading_vector = &semigroup.grading_vector;

    // The generators are not elements of the semigroup by index, so their
    // fingerprints and degrees are computed here.
    let generator_data: Vec<(Fingerprint, u32, usize)> = generators
        .column_iter()
        .enumerate()
        .map(|(g, c)| {
            let degree = (grading_vector * c)[0] as u32;
            (fingerprint_with(coefficients, c.as_slice()), degree, g)
        })
        .collect();

    for (i, col_i) in elements.column_iter().enumerate() {
        let degree = semigroup.degrees[i];
        let mut has_predecessor = degree == 0;
        for &(generator_fingerprint, generator_degree, g) in generator_data.iter() {
            let col_g = generators.column(g);
            if !has_predecessor && generator_degree <= degree {
                let difference = fingerprints[i].wrapping_sub(generator_fingerprint);
                if let Some(&candidate) = monomial_map.get(&difference) {
                    has_predecessor = elements
                        .column(candidate as usize)
                        .iter()
                        .zip(col_i.iter().zip(col_g.iter()))
                        .all(|(&c, (&a, &b))| c == a - b);
                }
            }
            if degree + generator_degree <= semigroup.max_degree {
                let sum = fingerprints[i].wrapping_add(generator_fingerprint);
                let Some(&candidate) = monomial_map.get(&sum) else {
                    return false;
                };
                if elements
                    .column(candidate as usize)
                    .iter()
                    .zip(col_i.iter().zip(col_g.iter()))
                    .any(|(&c, (&a, &b))| c != a + b)
                {
                    return false;
                }
            }
        }
        if !has_predecessor {
            return false;
        }
    }
    true
}
