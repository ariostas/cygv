#![doc = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/README.md"))]

pub mod cone;
pub mod dispatch;
#[cfg(feature = "cli")]
pub mod io;
pub mod modular;
pub mod phase;
#[cfg(feature = "python")]
pub mod python;
pub mod reference;

/// Which invariants to compute.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum InvariantKind {
    /// Gopakumar-Vafa invariants, which are integers.
    GV,
    /// Gromov-Witten invariants, which are rational.
    GW,
}

impl InvariantKind {
    /// The name of these invariants in the input and output documents.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::GV => "gv",
            Self::GW => "gw",
        }
    }
}

/// The kind of Calabi-Yau that the invariants are computed for.
///
/// Threefolds are treated separately throughout. Their intersection numbers are
/// canonicalized as unordered triplets, and their invariants are labeled by a curve
/// class alone. In higher dimensions the first index of the intersection numbers is
/// a reference surface, which labels the resulting invariants alongside the curve
/// class.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CYKind {
    /// A CY threefold.
    Threefold,
    /// A CY of dimension greater than three.
    Nfold,
}

impl CYKind {
    /// Selects the variant from the `is_threefold` flag of the input formats.
    pub fn from_is_threefold(is_threefold: bool) -> Self {
        if is_threefold {
            Self::Threefold
        } else {
            Self::Nfold
        }
    }

    /// Whether the CY is a threefold.
    pub fn is_threefold(self) -> bool {
        matches!(self, Self::Threefold)
    }
}

/// A stratum of a vex phase of the ambient toric variety.
///
/// In a vex phase (MacFadden–Sheridan, arXiv:2512.14817) some curves have negative
/// anticanonical degree, where the HKTY coefficient has a pole. It cancels against the class of
/// the stratum $V_S$ that such a curve forces into the CY, where $S$ is the set of toric
/// divisors the curve meets negatively (McAllister's prescription, via Wang,
/// arXiv:1910.14440). For fine fans on the boundary points of reflexive 4-polytopes $S$ is a
/// 3-cone of the fan and $V_S$ a toric curve, so each stratum is given by the three divisors
/// and the class $d_S$ of $V_S$, with $d_{S,t} = J_t \cdot V_S$ in the basis of the charges.
/// Only threefold hypersurfaces support strata.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Stratum {
    /// The indices of the three toric divisors, i.e. of rows of the charge matrix.
    pub divisors: [usize; 3],
    /// The class $d_S$, as fractions (numerator, positive denominator), one per Kähler
    /// parameter.
    pub class: Vec<(i64, i64)>,
}

/// Checks that the strata fit a threefold hypersurface with `n_divisors` toric divisors and
/// `h11` Kähler parameters, and that no two of them have the same divisors, returning what is
/// wrong otherwise.
pub(crate) fn check_strata(
    strata: &[Stratum],
    n_divisors: usize,
    h11: usize,
    cy_kind: CYKind,
    is_hypersurface: bool,
) -> Result<(), &'static str> {
    if strata.is_empty() {
        return Ok(());
    }
    if !cy_kind.is_threefold() || !is_hypersurface {
        return Err("strata are only supported for threefold hypersurfaces");
    }
    let mut seen = std::collections::HashSet::new();
    for s in strata {
        let mut d = s.divisors;
        d.sort_unstable();
        if d[2] >= n_divisors || d[0] == d[1] || d[1] == d[2] {
            return Err("the divisors of a stratum must be three distinct rows of the charges");
        }
        if !seen.insert(d) {
            return Err("two strata have the same divisors");
        }
        if s.class.len() != h11 {
            return Err("the class of a stratum must have one entry per Kähler parameter");
        }
        if s.class.iter().any(|&(_, den)| den <= 0) {
            return Err("the denominators of the class of a stratum must be positive");
        }
    }
    Ok(())
}

/// Checks that the nef partition, unless it is empty (a hypersurface), splits the
/// `n_divisors` toric divisors, i.e. the rows of the charges, into nonempty parts, which
/// makes its anticanonical rows sum to the anticanonical class, returning what is wrong
/// otherwise.
pub(crate) fn check_nefpart(
    nefpart: &[nalgebra::DVector<i32>],
    n_divisors: usize,
) -> Result<(), &'static str> {
    if nefpart.is_empty() {
        return Ok(());
    }
    let mut seen = vec![false; n_divisors];
    for part in nefpart {
        if part.is_empty() {
            return Err("the parts of the nef partition must not be empty");
        }
        for &r in part.iter() {
            let r = usize::try_from(r)
                .ok()
                .filter(|&r| r < n_divisors)
                .ok_or("the nef partition must consist of indices of rows of the charges")?;
            if core::mem::replace(&mut seen[r], true) {
                return Err("the parts of the nef partition must be disjoint");
            }
        }
    }
    if seen.contains(&false) {
        return Err("the nef partition must include every row of the charges");
    }
    Ok(())
}

/// Checks that the intersection numbers are indexed by the `h11` divisors of the basis (but
/// the first index of an n-fold, a reference surface), that each one is given in one
/// ordering of its indices only, and that some of them do not vanish, returning what is
/// wrong otherwise.
pub(crate) fn check_intnums(
    intnums: &std::collections::HashMap<(usize, usize, usize), i32>,
    h11: usize,
    cy_kind: CYKind,
) -> Result<(), &'static str> {
    let mut seen = std::collections::HashSet::new();
    for (&(i, j, k), &v) in intnums {
        let idx = match cy_kind {
            CYKind::Threefold => {
                let mut idx = [i, j, k];
                idx.sort_unstable();
                idx
            }
            CYKind::Nfold => [i, j.min(k), j.max(k)],
        };
        if idx[1..].iter().any(|&x| x >= h11) || (cy_kind.is_threefold() && idx[0] >= h11) {
            return Err(match cy_kind {
                CYKind::Threefold => "the indices of the intersection numbers must be below h11",
                CYKind::Nfold => {
                    "the last two indices of the intersection numbers must be below h11"
                }
            });
        }
        if v != 0 && !seen.insert(idx) {
            return Err("each intersection number must be given in one ordering of its indices");
        }
    }
    if seen.is_empty() {
        return Err("there must be at least one nonzero intersection number");
    }
    Ok(())
}

/// The likely causes of GV invariants that are not integers, for the error messages of both
/// pipelines.
pub(crate) const NON_INTEGER_CAUSES: &str = "this usually means that the input is \
    inconsistent (e.g. charges and intersection numbers of different geometries, or a nef \
    partition that is missing or wrong), or that the generators do not span the Mori cone or \
    one of its faces";

/// Which pipeline computes the invariants.
///
/// Both compute the same invariants exactly. The modular pipeline is much faster, but only
/// supports truncation by a maximum degree or by target points, not by a minimum number of
/// points.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Method {
    /// The modular pipeline when it supports the request, and the reference one otherwise.
    #[default]
    Auto,
    /// The modular pipeline ([`mod@modular`]). A request it does not support is an error.
    Modular,
    /// The reference pipeline ([`mod@reference`]), a faithful implementation of the HKTY
    /// procedure.
    Reference,
}

impl Method {
    /// The name of the method in the input formats.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Modular => "modular",
            Self::Reference => "reference",
        }
    }

    /// The method with the given name, ignoring case.
    pub fn from_name(name: &str) -> Option<Self> {
        [Self::Auto, Self::Modular, Self::Reference]
            .into_iter()
            .find(|m| m.as_str().eq_ignore_ascii_case(name))
    }
}

// Re-export the main structs
#[doc(inline)]
pub use reference::polynomial::{
    prettyprint::PrettyPrintPolynomial, properties::PolynomialProperties, Polynomial,
};
#[doc(inline)]
pub use reference::semigroup::Semigroup;

// Re-export main trait
#[doc(inline)]
pub use reference::polynomial::coefficient::PolynomialCoeff;

// Re-export python module
#[cfg(feature = "python")]
#[doc(inline)]
pub use python::cygv;

// Re-export the cone data
#[doc(inline)]
pub use cone::{CandidateCone, ConeData};

// Re-export the big-number crate of the results, so that callers can name its types and
// traits without depending on a matching version themselves.
pub use malachite;

// Re-export main functions
#[doc(inline)]
pub use dispatch::*;
#[doc(inline)]
pub use reference::run_hkty;

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::dvector;
    use std::collections::HashMap;

    #[test]
    fn test_check_nefpart() {
        assert_eq!(check_nefpart(&[], 6), Ok(()));
        assert_eq!(
            check_nefpart(&[dvector![3, 0, 1], dvector![2, 5, 4]], 6),
            Ok(())
        );
        for nefpart in [
            vec![dvector![0, 1, 2], dvector![3, 4]],
            vec![dvector![0, 1, 2], dvector![2, 3, 4, 5]],
            vec![dvector![0, 1, 2], dvector![3, 4, 6]],
            vec![dvector![0, 1, 2], dvector![3, 4, -5]],
            vec![dvector![0, 1, 2, 3, 4, 5], dvector![]],
        ] {
            assert!(check_nefpart(&nefpart, 6).is_err(), "{nefpart:?}");
        }
    }

    #[test]
    fn test_check_intnums() {
        let ok = HashMap::from([((0, 0, 1), 2), ((1, 1, 1), 0), ((0, 1, 1), 3)]);
        assert_eq!(check_intnums(&ok, 2, CYKind::Threefold), Ok(()));
        // The same triple in two orderings, unless one of them vanishes.
        let twice = HashMap::from([((0, 0, 1), 2), ((1, 0, 0), 2)]);
        assert!(check_intnums(&twice, 2, CYKind::Threefold).is_err());
        let zero = HashMap::from([((0, 0, 1), 2), ((1, 0, 0), 0)]);
        assert_eq!(check_intnums(&zero, 2, CYKind::Threefold), Ok(()));
        // For an n-fold the first index is a reference surface, and only the other two are
        // unordered.
        let nfold = HashMap::from([((0, 0, 1), 2), ((1, 0, 1), 3), ((5, 1, 1), 1)]);
        assert_eq!(check_intnums(&nfold, 2, CYKind::Nfold), Ok(()));
        assert!(check_intnums(&nfold, 2, CYKind::Threefold).is_err());
        let twice = HashMap::from([((0, 0, 1), 2), ((0, 1, 0), 2)]);
        assert!(check_intnums(&twice, 2, CYKind::Nfold).is_err());
        assert!(check_intnums(&HashMap::from([((0, 0, 2), 1)]), 2, CYKind::Nfold).is_err());
        assert!(check_intnums(&HashMap::from([((0, 0, 0), 0)]), 2, CYKind::Nfold).is_err());
        assert!(check_intnums(&HashMap::new(), 2, CYKind::Threefold).is_err());
    }
}
