#![doc = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/README.md"))]

pub mod cone;
mod dispatch;
mod interrupt;
#[cfg(feature = "cli")]
pub mod io;
pub mod modular;
pub mod phase;
mod problem;
#[cfg(feature = "python")]
pub mod python;
pub mod reference;

use malachite::Rational;

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
    /// Selects the variant by whether the CY is a threefold.
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
    /// The class $d_S$, one entry per Kähler parameter.
    pub class: Vec<Rational>,
}

impl Stratum {
    /// The class as fractions of 64-bit integers (numerator, positive denominator), which
    /// [`check_strata`] makes sure they are.
    pub(crate) fn class_i64(&self) -> Vec<(i64, i64)> {
        self.class
            .iter()
            .map(|x| fraction_i64(x).expect("the strata were checked"))
            .collect()
    }
}

/// The class of a stratum from fractions (numerator, denominator), for the tests.
#[cfg(test)]
pub(crate) fn class(fractions: &[(i64, i64)]) -> Vec<Rational> {
    fractions
        .iter()
        .map(|&(n, d)| Rational::from(n) / Rational::from(d))
        .collect()
}

/// A rational as a numerator and a positive denominator, if both fit in 64 bits.
pub(crate) fn fraction_i64(x: &Rational) -> Option<(i64, i64)> {
    let num = i64::try_from(x.numerator_ref()).ok()?;
    let den = i64::try_from(x.denominator_ref()).ok()?;
    Some((if *x < 0 { -num } else { num }, den))
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
        if s.class.iter().any(|x| fraction_i64(x).is_none()) {
            return Err(
                "the entries of the class of a stratum must be fractions of 64-bit integers",
            );
        }
    }
    Ok(())
}

/// Checks that the nef partition, unless it is empty (a hypersurface), splits the
/// `n_divisors` toric divisors, i.e. the rows of the charges, into nonempty parts, which
/// makes its anticanonical rows sum to the anticanonical class, returning what is wrong
/// otherwise.
pub(crate) fn check_nefpart(nefpart: &[Vec<usize>], n_divisors: usize) -> Result<(), &'static str> {
    if nefpart.is_empty() {
        return Ok(());
    }
    let mut seen = vec![false; n_divisors];
    for part in nefpart {
        if part.is_empty() {
            return Err("the parts of the nef partition must not be empty");
        }
        for &r in part {
            if r >= n_divisors {
                return Err("the nef partition must consist of indices of rows of the charges");
            }
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
/// the first index of an n-fold, a reference surface, numbered from 0 without gaps), that
/// each one is given in one
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
    // An n-fold's first index numbers the reference surfaces, and every table holds one
    // component per number up to the largest.
    if !cy_kind.is_threefold() {
        let surfaces: std::collections::BTreeSet<usize> = seen.iter().map(|idx| idx[0]).collect();
        if surfaces.last() != Some(&(surfaces.len() - 1)) {
            return Err(
                "the reference surfaces, the first indices of the intersection numbers, \
                 must be numbered from 0 without gaps",
            );
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
/// Both compute the same invariants exactly. The modular pipeline is much faster, but does
/// not support [`Selection::MinPoints`], and only takes curve classes whose coordinates fit
/// its keys; the reference pipeline cannot finish at high
/// $h^{1,1}$ anyway.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Method {
    /// The modular pipeline when the selection allows it, and the reference one otherwise.
    /// Curve classes too large for the modular keys are an error, which names
    /// [`Method::Reference`], rather than a slow run of it.
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

#[doc(inline)]
pub use dispatch::{
    compute, compute_gv, compute_gw, Computed, Error, Invariant, Invariants, Report,
};
#[doc(inline)]
pub use interrupt::Interrupt;
pub use problem::{Ambient, Geometry, MoriCone, Options, Problem, Selection};

#[doc(inline)]
pub use cone::CandidateCone;

// Re-export python module
#[cfg(feature = "python")]
#[doc(inline)]
pub use python::cygv;

/// Makes the allocator return large blocks to the system as soon as they are freed, which
/// keeps the resident memory of a run close to what it holds: about a third less in the
/// low-memory mode ([`Options::low_memory`]), for a few percent of time. It changes a setting
/// of the whole process, glibc's `M_MMAP_THRESHOLD`, which otherwise rises as large blocks
/// are freed, so that later ones come from a heap that keeps their memory. That also makes
/// every later allocation of 128 KiB or more a system call, which slows code that makes many
/// (numpy's, five times on arrays of 1.6 MB). The command line and the Python module call it
/// for a run in the low-memory mode, and the Python module undoes it afterwards as far as
/// glibc allows; it does nothing on other platforms than Linux with glibc.
pub fn return_freed_memory() {
    #[cfg(all(target_os = "linux", target_env = "gnu"))]
    // SAFETY: mallopt only changes a setting of glibc's allocator, which is safe to do at any
    // time.
    unsafe {
        glibc::mallopt(glibc::M_MMAP_THRESHOLD, 128 * 1024);
    }
}

/// Undoes [`return_freed_memory`] as far as glibc allows: its threshold no longer moves once
/// set, so this sets it, and the threshold above which the heap is trimmed, to the most they
/// rise to by themselves (32 MiB, and twice that).
#[cfg(feature = "python")]
pub(crate) fn keep_freed_memory() {
    #[cfg(all(target_os = "linux", target_env = "gnu"))]
    // SAFETY: as in return_freed_memory.
    unsafe {
        glibc::mallopt(glibc::M_MMAP_THRESHOLD, 32 << 20);
        glibc::mallopt(glibc::M_TRIM_THRESHOLD, 64 << 20);
    }
}

/// Gives the memory that the allocator holds free back to the system. A computation leaves
/// much of it in the arenas of its threads, which have exited: on a run that peaks at 1.26 GB,
/// 578 MB stay resident once its results are dropped, and 274 MB after this, which takes 25
/// ms. The command line calls it after each document and the Python module after each
/// computation; it does nothing on other platforms than Linux with glibc.
pub fn trim_freed_memory() {
    #[cfg(all(target_os = "linux", target_env = "gnu"))]
    // SAFETY: malloc_trim only releases free memory, and may be called at any time.
    unsafe {
        glibc::malloc_trim(0);
    }
}

/// The settings of glibc's allocator.
#[cfg(all(target_os = "linux", target_env = "gnu"))]
mod glibc {
    use core::ffi::c_int;

    extern "C" {
        pub fn mallopt(param: c_int, value: c_int) -> c_int;
        pub fn malloc_trim(pad: usize) -> c_int;
    }

    #[cfg(feature = "python")]
    pub const M_TRIM_THRESHOLD: c_int = -1;
    pub const M_MMAP_THRESHOLD: c_int = -3;
}

/// Writes the progress messages of the computations to the standard error, each line
/// stamped with the seconds since this was first called: with a `verbosity` of 1 the stages
/// and passes, with 2 also the degrees of each pass, and with 0 nothing.
///
/// The messages go through the [`log`](https://docs.rs/log) facade, at the levels `info`
/// and `debug`, so library users can route them through any logger instead. This one is
/// what the command line's `--verbose` installs (the Python keyword `verbose` has a logger of
/// its own); once a logger is set, later calls only change the verbosity, which is that of
/// the whole process (`log::set_max_level`).
pub fn log_to_stderr(verbosity: u8) {
    use std::io::Write;
    use std::time::Instant;

    struct Stderr(Instant);

    impl log::Log for Stderr {
        fn enabled(&self, metadata: &log::Metadata) -> bool {
            metadata.target().starts_with("cygv")
        }

        fn log(&self, record: &log::Record) {
            if self.enabled(record.metadata()) {
                let secs = self.0.elapsed().as_secs_f64();
                // A failed write to the standard error has nowhere to be reported.
                let _ = writeln!(std::io::stderr(), "[{secs:9.2}s] {}", record.args());
            }
        }

        fn flush(&self) {}
    }

    static LOGGER: std::sync::OnceLock<Stderr> = std::sync::OnceLock::new();
    let logger = LOGGER.get_or_init(|| Stderr(Instant::now()));
    // Fails only when a logger is already set, which keeps it.
    let _ = log::set_logger(logger);
    log::set_max_level(match verbosity {
        0 => log::LevelFilter::Off,
        1 => log::LevelFilter::Info,
        _ => log::LevelFilter::Debug,
    });
}

// Re-export the big-number crate of the results and the matrix crate of the input, so that
// callers can name their types and traits without depending on a matching version themselves.
pub use malachite;
pub use nalgebra;

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_check_nefpart() {
        assert_eq!(check_nefpart(&[], 6), Ok(()));
        assert_eq!(check_nefpart(&[vec![3, 0, 1], vec![2, 5, 4]], 6), Ok(()));
        for nefpart in [
            vec![vec![0, 1, 2], vec![3, 4]],
            vec![vec![0, 1, 2], vec![2, 3, 4, 5]],
            vec![vec![0, 1, 2], vec![3, 4, 6]],
            vec![vec![0, 1, 2, 3, 4, 5], vec![]],
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
        let nfold = HashMap::from([((0, 0, 1), 2), ((1, 0, 1), 3), ((2, 1, 1), 1)]);
        assert_eq!(check_intnums(&nfold, 2, CYKind::Nfold), Ok(()));
        // The surfaces are numbered without gaps, since each number takes a component.
        let gap = HashMap::from([((0, 0, 1), 2), ((5, 1, 1), 1), ((3, 1, 1), 0)]);
        assert!(check_intnums(&gap, 2, CYKind::Nfold).is_err());
        assert!(check_intnums(&nfold, 2, CYKind::Threefold).is_err());
        let twice = HashMap::from([((0, 0, 1), 2), ((0, 1, 0), 2)]);
        assert!(check_intnums(&twice, 2, CYKind::Nfold).is_err());
        assert!(check_intnums(&HashMap::from([((0, 0, 2), 1)]), 2, CYKind::Nfold).is_err());
        assert!(check_intnums(&HashMap::from([((0, 0, 0), 0)]), 2, CYKind::Nfold).is_err());
        assert!(check_intnums(&HashMap::new(), 2, CYKind::Threefold).is_err());
    }
}
