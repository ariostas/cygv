//! What to compute: a [`Problem`], which is a Calabi-Yau ([`Geometry`]), the cone its curve
//! classes live in ([`MoriCone`]) and which of them to compute ([`Selection`]), and the
//! [`Options`] of the computation.
//!
//! Every field has the name of the Python keyword and of the YAML key that set it, and every
//! choice between alternatives is an enum, so that only one of them can be given.

use crate::cone::{check_cone_data, columns, CandidateCone, ConeData};
use crate::dispatch::Error;
use crate::phase::{Fan, FanError};
use crate::{check_intnums, check_nefpart, check_strata, modular, CYKind, Method, Stratum};
use nalgebra::DMatrix;
use std::collections::HashMap;

/// A Calabi-Yau hypersurface or complete intersection in a toric variety.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Geometry {
    /// The GLSM charges: one row per toric divisor, one column per Kähler parameter. (In
    /// Python and YAML, `q` is given as one charge vector per Kähler parameter, i.e. by its
    /// columns.) A divisor index, as in `nefpart`, the fan or the strata, is a row of `q`.
    pub q: DMatrix<i32>,
    /// The intersection numbers $\kappa_{ijk}$ in the basis of the Kähler parameters, each
    /// given in one ordering of its indices (the others are implied, and those not given are
    /// zero). For an n-fold, the first index is a reference surface instead, and only the last
    /// two are unordered.
    pub intnums: HashMap<(usize, usize, usize), i32>,
    /// The nef partition of a complete intersection: disjoint parts that cover every toric
    /// divisor. Empty for a hypersurface.
    pub nefpart: Vec<Vec<usize>>,
    /// What else is known about the ambient toric variety.
    pub ambient: Ambient,
}

/// What is known about the ambient toric variety beyond its charges.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum Ambient {
    /// Nothing: a phase without strata is assumed.
    #[default]
    Plain,
    /// The maximal cones of its fan, as divisor indices (threefold hypersurfaces only; see
    /// [`Fan`]). The strata of a vex phase are derived from it, the Mori cone can be
    /// ([`MoriCone::FromFan`]), and so is the default grading vector; the intersection
    /// numbers and a given grading vector are checked against it.
    Fan(Vec<Vec<usize>>),
    /// The strata of a vex phase, for callers that have them but not the fan (threefold
    /// hypersurfaces only).
    Strata(Vec<Stratum>),
}

/// The Mori cone, or a face of it: the cone the curve classes live in. Every lattice point of
/// it is used, even when the generators do not generate them all.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MoriCone {
    /// Vectors that span it, one per column: its rays, or rays and lattice points in it.
    Generators(DMatrix<i32>),
    /// The inequalities $h \cdot C \ge 0$ that cut it out, one per column (redundant ones are
    /// allowed), and, for a face, the equations $e \cdot C = 0$ of its span, one per column.
    Inequalities {
        hyperplanes: DMatrix<i32>,
        equations: DMatrix<i32>,
    },
    /// The cone spanned by the classes of the walls of [`Ambient::Fan`].
    FromFan,
}

/// Which curve classes to compute the invariants of.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Selection {
    /// Those of degree at most this.
    MaxDeg(u32),
    /// Those in the backward lightcones of the target points (one per column), i.e. the $C$
    /// with $P - C$ in the Mori cone for some target $P$, up to the largest degree of the
    /// targets, or `max_deg` if it is lower.
    TargetPoints {
        points: DMatrix<i32>,
        max_deg: Option<u32>,
    },
    /// The lowest degrees that hold at least this many classes (reference pipeline only).
    MinPoints(u32),
    /// Exactly these classes, one per column, which must lie in the Mori cone and have
    /// positive degree; zero is added if it is missing (reference pipeline only). The results
    /// are only right if every class below one of them is listed too.
    CurveClasses(DMatrix<i32>),
}

impl Selection {
    /// Whether the modular pipeline supports this selection.
    pub(crate) fn modular_supported(&self) -> bool {
        matches!(self, Self::MaxDeg(_) | Self::TargetPoints { .. })
    }

    /// The name of this selection in the input formats.
    pub(crate) fn name(&self) -> &'static str {
        match self {
            Self::MaxDeg(_) => "max_deg",
            Self::TargetPoints { .. } => "target_points",
            Self::MinPoints(_) => "min_points",
            Self::CurveClasses(_) => "curve_classes",
        }
    }
}

/// A computation of invariants: the Calabi-Yau, its Mori cone, which curve classes, and the
/// grading vector that defines their degrees.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Problem {
    pub geometry: Geometry,
    pub cone: MoriCone,
    pub selection: Selection,
    /// The degree of a curve class $C$ is $g \cdot C$, which must be positive on the nonzero
    /// points of the Mori cone. By default, the primitive sum of the inequalities of the cone,
    /// which is positive on it because it is pointed.
    pub grading_vector: Option<Vec<i32>>,
}

/// How to compute the invariants. None of these change the results.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Options {
    /// Which pipeline to use.
    pub method: Method,
    /// The number of worker threads, at least one, or `None` for one per core.
    pub n_threads: Option<usize>,
    /// Whether to prove the results exact (see [`modular::Options::certify`]); the reference
    /// pipeline is exact anyway. Not supported in vex phases yet.
    pub certify: bool,
    /// Whether to prefer memory to speed (see [`modular::Options::low_memory`]), for the runs
    /// that do not fit otherwise: about a third of the memory for about two and a half times
    /// the time. Only the modular pipeline has such a mode; the reference ignores it.
    pub low_memory: bool,
    /// The Hilbert bases of the cones $K_T$ (see [`crate::cone::ConeData`]), computed by an
    /// external program such as Normaliz, which the modular pipeline then uses instead of
    /// computing them (the reference pipeline does not need them). They are trusted: an
    /// element missing from them drops curve classes.
    pub hilbert_bases: Option<Vec<CandidateCone>>,
    /// The tuning of the modular pipeline.
    pub tuning: modular::Tuning,
}

/// What the checks of a problem find out on the way.
pub(crate) struct Checked {
    pub cy_kind: CYKind,
    /// The strata: those given, those of the fan, or none.
    pub strata: Vec<Stratum>,
    pub fan: Option<Fan>,
}

/// An invalid input, for the field named.
pub(crate) fn input_err(field: &'static str, reason: impl Into<String>) -> Error {
    Error::Input {
        field,
        reason: reason.into(),
    }
}

/// Checks the Hilbert bases of [`Options::hilbert_bases`] against the Mori cone, the charges,
/// the generators (if the cone is given by them) and the strata.
pub(crate) fn check_hilbert_bases(
    bases: &[CandidateCone],
    hyperplanes: &DMatrix<i32>,
    equations: &DMatrix<i32>,
    q: &DMatrix<i32>,
    generators: Option<&DMatrix<i32>>,
    strata: &[Stratum],
) -> Result<ConeData, Error> {
    let cone_data = ConeData {
        hyperplanes: hyperplanes.clone(),
        equations: if equations.ncols() == 0 {
            DMatrix::zeros(q.ncols(), 0)
        } else {
            equations.clone()
        },
        cones: bases.to_vec(),
    };
    let none = DMatrix::zeros(q.ncols(), 0);
    check_cone_data(&cone_data, q, generators.unwrap_or(&none), strata).map_err(
        |(field, reason)| {
            let field = match field {
                "cones" => "hilbert_bases",
                f => f,
            };
            input_err(field, reason)
        },
    )?;
    Ok(cone_data)
}

/// Checks that the target points have positive degree under the grading vector.
pub(crate) fn check_target_degrees(
    points: &DMatrix<i32>,
    grading_vector: &[i32],
) -> Result<(), Error> {
    let positive = columns(points).all(|p| {
        p.iter()
            .zip(grading_vector)
            .map(|(&x, &w)| x as i64 * w as i64)
            .sum::<i64>()
            > 0
    });
    if !positive {
        return Err(input_err(
            "target_points",
            "the target points must have positive degree under the grading vector",
        ));
    }
    Ok(())
}

/// The input field of the Mori cone, for the error messages.
pub(crate) fn cone_field(cone: &MoriCone) -> &'static str {
    match cone {
        MoriCone::Inequalities { .. } => "mori_hyperplanes",
        _ => "generators",
    }
}

impl Geometry {
    /// The number of Kähler parameters.
    pub fn h11(&self) -> usize {
        self.q.ncols()
    }

    /// The dimension of the Calabi-Yau that the charges and the nef partition imply.
    pub fn dimension(&self) -> i64 {
        self.q.nrows() as i64 - self.q.ncols() as i64 - self.nefpart.len().max(1) as i64
    }

    /// The kind of Calabi-Yau, once [`Problem::validate`] has accepted it.
    pub(crate) fn cy_kind(&self) -> CYKind {
        CYKind::from_is_threefold(self.dimension() == 3)
    }
}

impl Problem {
    /// Checks everything about the problem and the options that can be checked before the
    /// computation, returning the first error, which names the field at fault. The
    /// computation checks them anyway; this is for callers that want to know before they
    /// start any of several computations, like the command line interface.
    pub fn validate(&self, options: &Options) -> Result<(), Error> {
        self.check(options).map(drop)
    }

    pub(crate) fn check(&self, options: &Options) -> Result<Checked, Error> {
        let Geometry {
            q,
            intnums,
            nefpart,
            ambient,
        } = &self.geometry;
        let (n_divisors, h11) = q.shape();
        if n_divisors == 0 || h11 == 0 {
            return Err(input_err("q", "the charges must not be empty"));
        }
        check_nefpart(nefpart, n_divisors).map_err(|r| input_err("nefpart", r))?;
        let dim = self.geometry.dimension();
        if dim < 3 {
            // One charge vector per divisor instead of one per Kähler parameter.
            let transposed = n_divisors as i64 - h11 as i64;
            let hint = if n_divisors < h11 && -transposed - nefpart.len().max(1) as i64 >= 3 {
                format!(
                    "; the charges have {h11} vectors of {n_divisors} entries, which looks \
                     transposed: they need one vector per Kähler parameter, with one entry per \
                     toric divisor"
                )
            } else {
                String::new()
            };
            return Err(input_err(
                "q",
                format!("the dimension of the CY must be at least three, but it is {dim}{hint}"),
            ));
        }
        let cy_kind = self.geometry.cy_kind();
        check_intnums(intnums, h11, cy_kind).map_err(|r| input_err("intnums", r))?;
        let (strata, fan) = match ambient {
            Ambient::Plain => (Vec::new(), None),
            Ambient::Strata(strata) => {
                check_strata(strata, n_divisors, h11, cy_kind, nefpart.is_empty())
                    .map_err(|r| input_err("strata", r))?;
                (strata.clone(), None)
            }
            Ambient::Fan(cones) => {
                let fan_err = |e: FanError| match e {
                    FanError::InvalidIntnums | FanError::IntnumsMismatch(_) => {
                        input_err("intnums", e.to_string())
                    }
                    e => input_err("fan", e.to_string()),
                };
                if !cy_kind.is_threefold() || !nefpart.is_empty() {
                    return Err(fan_err(FanError::NotHypersurface));
                }
                let fan = Fan::new(cones, q).map_err(fan_err)?;
                fan.check_intnums(intnums).map_err(fan_err)?;
                (fan.strata().map_err(fan_err)?, Some(fan))
            }
        };
        match &self.cone {
            MoriCone::Generators(g) => {
                if g.ncols() == 0 {
                    return Err(input_err("generators", "the generators must not be empty"));
                }
                if g.nrows() != h11 {
                    return Err(input_err(
                        "generators",
                        "the generators must have one entry per Kähler parameter",
                    ));
                }
            }
            MoriCone::Inequalities {
                hyperplanes,
                equations,
            } => {
                if hyperplanes.ncols() == 0 {
                    return Err(input_err(
                        "mori_hyperplanes",
                        "the Mori cone needs at least one hyperplane",
                    ));
                }
                if hyperplanes.nrows() != h11 {
                    return Err(input_err(
                        "mori_hyperplanes",
                        "the hyperplanes must have one entry per Kähler parameter",
                    ));
                }
                if equations.ncols() > 0 && equations.nrows() != h11 {
                    return Err(input_err(
                        "mori_equations",
                        "the equations must have one entry per Kähler parameter",
                    ));
                }
            }
            MoriCone::FromFan if fan.is_none() => {
                return Err(input_err(
                    "generators",
                    "the Mori cone is needed: give its generators or its inequalities \
                     (mori_hyperplanes), or the fan to derive it from",
                ))
            }
            MoriCone::FromFan => {}
        }
        if let Some(g) = &self.grading_vector {
            if g.len() != h11 {
                return Err(input_err(
                    "grading_vector",
                    "the grading vector must have one entry per Kähler parameter",
                ));
            }
            if let Some(fan) = &fan {
                fan.check_grading(&nalgebra::RowDVector::from_row_slice(g))
                    .map_err(|e| input_err("grading_vector", e.to_string()))?;
            }
            if let MoriCone::Generators(gens) = &self.cone {
                let degree = |c: &[i32]| -> i64 {
                    c.iter().zip(g).map(|(&x, &w)| x as i64 * w as i64).sum()
                };
                if columns(gens).any(|c| c.iter().any(|&x| x != 0) && degree(c) <= 0) {
                    return Err(input_err(
                        "grading_vector",
                        "the grading vector must be positive on every nonzero generator",
                    ));
                }
            }
        }
        match &self.selection {
            Selection::MaxDeg(_) | Selection::MinPoints(_) => {}
            Selection::TargetPoints { points, .. } => {
                if points.ncols() == 0 {
                    return Err(input_err(
                        "target_points",
                        "there must be at least one target point",
                    ));
                }
                if points.nrows() != h11 {
                    return Err(input_err(
                        "target_points",
                        "the target points must have one entry per Kähler parameter",
                    ));
                }
            }
            Selection::CurveClasses(classes) => {
                if classes.ncols() == 0 || classes.nrows() != h11 {
                    return Err(input_err(
                        "curve_classes",
                        "the curve classes must be nonempty, with one entry per Kähler \
                         parameter",
                    ));
                }
            }
        }
        if options.method == Method::Modular && !self.selection.modular_supported() {
            return Err(Error::Unsupported(format!(
                "the modular method does not support {}; use the auto or reference method",
                self.selection.name()
            )));
        }
        if options.n_threads == Some(0) {
            return Err(input_err(
                "n_threads",
                "the number of threads must be positive",
            ));
        }
        let modular = options.method != Method::Reference && self.selection.modular_supported();
        if options.certify && !strata.is_empty() && modular {
            return Err(Error::Unsupported(
                "the invariants of a vex phase cannot be certified yet; use the reference \
                 method, which is exact, or leave out certify"
                    .to_owned(),
            ));
        }
        if let Some(bases) = &options.hilbert_bases {
            if bases
                .iter()
                .any(|c| c.hilbert_basis.ncols() > 0 && c.hilbert_basis.nrows() != h11)
            {
                return Err(input_err(
                    "hilbert_bases",
                    "the Hilbert basis elements must have one entry per Kähler parameter",
                ));
            }
            // With the inequalities at hand they are checked here, and otherwise once the
            // computation has them.
            if let MoriCone::Inequalities {
                hyperplanes,
                equations,
            } = &self.cone
            {
                check_hilbert_bases(bases, hyperplanes, equations, q, None, &strata)?;
            }
        }
        if let (Selection::TargetPoints { points, .. }, Some(g)) =
            (&self.selection, &self.grading_vector)
        {
            check_target_degrees(points, g)?;
        }
        Ok(Checked {
            cy_kind,
            strata,
            fan,
        })
    }
}
