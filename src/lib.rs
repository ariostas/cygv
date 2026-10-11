#![doc = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/README.md"))]

pub mod cone;
pub mod dispatch;
#[cfg(feature = "cli")]
pub mod io;
pub mod modular;
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

/// The likely causes of GV invariants that are not integers, for the error messages of both
/// pipelines.
pub(crate) const NON_INTEGER_CAUSES: &str = "this usually means that the input is \
    inconsistent (e.g. charges and intersection numbers of different geometries, or a nef \
    partition that is missing or wrong), or that the generators do not span the Mori cone or \
    one of its faces";

/// Which pipeline computes the invariants.
///
/// Both compute the same invariants exactly. The modular pipeline is much faster, but only
/// supports truncation by a maximum degree.
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
