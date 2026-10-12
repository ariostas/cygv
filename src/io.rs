//! Input and output functions
//!
//! This module implements the YAML-based data format used by the `cygv` command line
//! interface. An input stream contains one or more YAML documents, each of which
//! specifies a single computation.
//!
//! ```yaml
//! ---
//! # Optional label that is echoed back in the output.
//! name: quintic-like threefold
//! # The generators of the Mori cone, given as a list of vectors.
//! generators: [[0, -1], [1, 2]]
//! # A vector that has a positive inner product with every generator.
//! grading_vector: [3, -1]
//! # The GLSM charge matrix, given as a list of charge vectors.
//! q: [[1, 1, 1, 0, 1, 2], [0, 0, -1, 1, 1, -1]]
//! # The intersection numbers, given as a list of [i, j, k, value] entries.
//! intnums:
//!   - [0, 0, 0, 2]
//!   - [0, 0, 1, 1]
//!   - [0, 1, 1, -1]
//!   - [1, 1, 1, 5]
//! # Whether to compute GV or GW invariants. Defaults to gv.
//! invariants: gv
//! # How many curve classes to compute: those up to a maximum degree (max_deg), at
//! # least a number of them (min_points), or those below target points
//! # (target_points), i.e. that leave a target in the Mori cone when subtracted.
//! # max_deg may be given with target_points, and lowers the largest degree of the
//! # targets. When none of them is given, the generators are taken to be the
//! # complete list of curve classes to use.
//! min_points: 100
//! # Which pipeline to use: auto, modular or reference. Defaults to auto, which
//! # uses the modular one, which is much faster, whenever max_deg or target_points
//! # is given.
//! method: auto
//! # Whether to prove the results of the modular pipeline exact, taking more primes if
//! # needed (not in vex phases). Defaults to false.
//! certify: false
//! ```
//!
//! A vex phase of the ambient toric variety also needs its strata (see [`Stratum`]), each
//! given by three divisors and a class whose entries may be fractions.
//!
//! ```yaml
//! strata:
//!   - divisors: [1, 2, 4]
//!     class: [1, "-1/2"]
//! ```
//!
//! They can instead be derived from the fan of the ambient toric variety (see [`Fan`]), given by
//! its maximal cones as indices of the divisors. The fan also gives the Mori cone and a grading
//! vector, so `generators` and `grading_vector` can be left out (the generators only when the
//! computation is truncated), and the intersection numbers and a grading vector that is given
//! are checked against it.
//!
//! ```yaml
//! fan: [[0, 1, 2, 3], [0, 1, 2, 4], [0, 1, 3, 4], [0, 2, 3, 4], [1, 2, 3, 5], [1, 2, 4, 5],
//!   [1, 3, 4, 5], [2, 3, 4, 5]]
//! ```
//!
//! The Mori cone can also be given by cone data (see [`ConeData`]): its inequalities, the
//! equations of its span if it is a face, and the Hilbert basis of the cone $K_T$ of every set
//! $T$ of one or two divisors, which an external program such as normaliz computes. Then the
//! generators can be left out. Without it, the cone data is computed from the generators.
//!
//! ```yaml
//! mori_hyperplanes: [[1, 0], [0, 1]]
//! cones:
//!   - divisors: [0]
//!     hilbert_basis: [[1, 0], [2, 1]]
//!   - divisors: [0, 5]
//!     hilbert_basis: [[1, 0], [0, 1]]
//!   # ... one entry per set of one or two divisors
//! ```
//!
//! The results are written back as one YAML document per input document.
//!
//! ```yaml
//! ---
//! name: quintic-like threefold
//! invariants: gv
//! is_threefold: true
//! grading_vector: [3, -1]
//! results:
//!   - {curve_class: [0, -1], degree: 1, gv: -4}
//!   - {curve_class: [1, 0], degree: 3, gv: 7524}
//! ```
//!
//! See [`Input`] for the full list of fields.

pub mod error;

use crate::cone::{check_cone_data, CandidateCone, ConeData};
use crate::dispatch::{
    check_target_points, compute_gvgw_strings, modular_unsupported, ComputeError,
};
use crate::phase::{Fan, FanError};
use crate::reference::misc::process_int_nums;
use crate::{check_strata, CYKind, InvariantKind, Method, Stratum};
use error::IoError;
use nalgebra::{DMatrix, DVector, RowDVector};
use std::collections::HashMap;
use std::io::Write;
use yaml_rust2::yaml::{Hash, Yaml};
use yaml_rust2::YamlLoader;

/// The fields that are accepted in an input document.
const INPUT_FIELDS: [&str; 19] = [
    "name",
    "generators",
    "grading_vector",
    "q",
    "intnums",
    "nefpart",
    "invariants",
    "max_deg",
    "min_points",
    "target_points",
    "is_threefold",
    "n_threads",
    "method",
    "certify",
    "strata",
    "fan",
    "mori_hyperplanes",
    "mori_equations",
    "cones",
];

/// The specification of a single computation of GV or GW invariants.
///
/// The matrices are stored using the same convention as the rest of the crate, i.e.
/// the vectors of the input data are the *columns* of the matrices.
#[derive(Debug, Clone, PartialEq)]
pub struct Input {
    /// An optional label that is echoed back in the output.
    pub name: Option<String>,
    /// The generators of the Mori cone. It has no columns when they are left out, which the
    /// cone data allows.
    pub generators: DMatrix<i32>,
    /// A vector that has a positive inner product with every generator.
    pub grading_vector: RowDVector<i32>,
    /// The GLSM charge matrix.
    pub q: DMatrix<i32>,
    /// The nef partition. It is empty for hypersurfaces.
    pub nefpart: Vec<DVector<i32>>,
    /// The triple intersection numbers.
    pub intnums: HashMap<(usize, usize, usize), i32>,
    /// The strata of a vex phase. It is empty otherwise.
    pub strata: Vec<Stratum>,
    /// The maximal cones of the fan of the ambient toric variety, as indices of rows of `q`,
    /// from which the strata are derived (see [`Fan`]).
    pub fan: Option<Vec<Vec<usize>>>,
    /// The cone data for the modular method, if any.
    pub cone_data: Option<ConeData>,
    /// Which invariants to compute.
    pub invariant_kind: InvariantKind,
    /// The kind of CY the invariants are computed for.
    pub cy_kind: CYKind,
    /// The maximum degree of the curve classes that are computed.
    pub max_deg: Option<u32>,
    /// The minimum number of curve classes that are computed.
    pub min_points: Option<u32>,
    /// The curve classes that must be included in the computation.
    pub target_points: Option<DMatrix<i32>>,
    /// The number of threads to use. It is deduced from the machine when [`None`].
    pub n_threads: Option<u32>,
    /// Which pipeline to use.
    pub method: Method,
    /// Whether to prove the results of the modular pipeline exact (see
    /// [`crate::modular::Options::certify`]).
    pub certify: bool,
}

/// A computed invariant, consisting of a curve class, a reference surface index that
/// is only meaningful when the CY is not a threefold, and the invariant itself.
pub type InvariantResult = ((DVector<i32>, usize), String);

impl Input {
    /// Parse every YAML document of the given string into an [`Input`].
    pub fn load_all(data: &str) -> Result<Vec<Self>, IoError> {
        let docs =
            YamlLoader::load_from_str(data).map_err(|e| IoError::InvalidYaml(e.to_string()))?;
        if docs.is_empty() {
            return Err(IoError::EmptyInput);
        }
        docs.iter().map(Input::from_yaml).collect()
    }

    /// Parse a single YAML document into an [`Input`].
    pub fn from_yaml(doc: &Yaml) -> Result<Self, IoError> {
        let hash = doc.as_hash().ok_or(IoError::NotAMapping)?;
        for key in hash.keys() {
            match key.as_str() {
                Some(k) if INPUT_FIELDS.contains(&k) => {}
                Some(k) => return Err(IoError::UnknownField(k.to_owned())),
                None => return Err(IoError::UnknownField(format!("{key:?}"))),
            }
        }

        let name = match get_field(hash, "name") {
            None => None,
            Some(y) => Some(
                y.as_str()
                    .ok_or_else(|| IoError::invalid_field("name", "expected a string"))?
                    .to_owned(),
            ),
        };
        let generators = get_field(hash, "generators")
            .map(|y| as_matrix(y, "generators"))
            .transpose()?;
        let grading_vector = get_field(hash, "grading_vector")
            .map(|y| as_int_vec(y, "grading_vector"))
            .transpose()?;
        let fan = get_field(hash, "fan").map(as_fan).transpose()?;
        if grading_vector.is_none() && fan.is_none() {
            return Err(IoError::MissingField("grading_vector"));
        }
        let q = as_matrix(require_field(hash, "q")?, "q")?;
        // Without a grading vector, which the fan gives a default for, the charges give h11.
        let h11 = grading_vector.as_ref().map_or(q.ncols(), Vec::len);
        if let Some(generators) = &generators {
            if generators.nrows() != h11 {
                let field = if grading_vector.is_some() {
                    "grading_vector"
                } else {
                    "q"
                };
                return Err(IoError::invalid_field(
                    field,
                    format!(
                        "expected {} entries to match the generators",
                        generators.nrows()
                    ),
                ));
            }
        }
        if q.ncols() != h11 {
            return Err(IoError::invalid_field(
                "q",
                format!("expected {h11} charge vectors to match the grading vector"),
            ));
        }
        let intnums = as_intnums(require_field(hash, "intnums")?)?;
        let nefpart = match get_field(hash, "nefpart") {
            None => Vec::new(),
            Some(y) => as_vectors(y, "nefpart")?,
        };
        let strata = match get_field(hash, "strata") {
            None => Vec::new(),
            Some(y) => as_strata(y)?,
        };
        let cone_data = match (
            get_field(hash, "mori_hyperplanes"),
            get_field(hash, "cones"),
        ) {
            (Some(hyperplanes), Some(cones)) => Some(ConeData {
                hyperplanes: as_matrix(hyperplanes, "mori_hyperplanes")?,
                equations: match get_field(hash, "mori_equations") {
                    None => DMatrix::zeros(h11, 0),
                    Some(y) => as_matrix(y, "mori_equations")?,
                },
                cones: as_cones(cones, h11)?,
            }),
            (None, None) if get_field(hash, "mori_equations").is_none() => None,
            (None, _) => return Err(IoError::MissingField("mori_hyperplanes")),
            (_, None) => return Err(IoError::MissingField("cones")),
        };

        let invariant_kind = match get_field(hash, "invariants") {
            None => InvariantKind::GV,
            Some(y) => match y.as_str().map(str::to_lowercase).as_deref() {
                Some("gv") => InvariantKind::GV,
                Some("gw") => InvariantKind::GW,
                _ => {
                    return Err(IoError::invalid_field(
                        "invariants",
                        "expected either \"gv\" or \"gw\"",
                    ))
                }
            },
        };
        let max_deg = get_field(hash, "max_deg")
            .map(|y| as_u32(y, "max_deg"))
            .transpose()?;
        let min_points = get_field(hash, "min_points")
            .map(|y| as_u32(y, "min_points"))
            .transpose()?;
        let target_points = get_field(hash, "target_points")
            .map(|y| as_target_points(y, h11))
            .transpose()?;
        let n_threads = get_field(hash, "n_threads")
            .map(|y| as_u32(y, "n_threads"))
            .transpose()?;
        let certify = match get_field(hash, "certify") {
            None => false,
            Some(y) => y
                .as_bool()
                .ok_or_else(|| IoError::invalid_field("certify", "expected a boolean"))?,
        };
        let method = match get_field(hash, "method") {
            None => Method::default(),
            Some(y) => y.as_str().and_then(Method::from_name).ok_or_else(|| {
                IoError::invalid_field(
                    "method",
                    "expected one of \"auto\", \"modular\" or \"reference\"",
                )
            })?,
        };

        let cy_kind = match get_field(hash, "is_threefold") {
            None => infer_cy_kind(&q, &nefpart),
            Some(y) => CYKind::from_is_threefold(
                y.as_bool()
                    .ok_or_else(|| IoError::invalid_field("is_threefold", "expected a boolean"))?,
            ),
        };
        if min_points.is_some() && (max_deg.is_some() || target_points.is_some()) {
            return Err(IoError::ConflictingTruncations);
        }
        // The fan is checked again by the computation, but the default grading vector is
        // needed here.
        let grading_vector = match &fan {
            None => RowDVector::from_vec(grading_vector.unwrap_or_default()),
            Some(cones) => {
                let fan_err = |e: FanError| IoError::invalid_field("fan", e.to_string());
                if !cy_kind.is_threefold() || !nefpart.is_empty() {
                    return Err(fan_err(FanError::NotHypersurface));
                }
                if !strata.is_empty() {
                    return Err(IoError::invalid_field(
                        "strata",
                        "the strata are derived from the fan, so they cannot be given with it",
                    ));
                }
                let fan = Fan::new(cones, &q).map_err(fan_err)?;
                fan.check_intnums(&intnums).map_err(|e| match e {
                    FanError::InvalidIntnums => IoError::invalid_field("intnums", e.to_string()),
                    e => fan_err(e),
                })?;
                fan.strata().map_err(fan_err)?;
                match grading_vector {
                    Some(g) => {
                        let g = RowDVector::from_vec(g);
                        fan.check_grading(&g).map_err(|e| match e {
                            FanError::NonpositiveGrading => {
                                IoError::invalid_field("grading_vector", e.to_string())
                            }
                            e => fan_err(e),
                        })?;
                        g
                    }
                    None => fan.mori_cone().map_err(fan_err)?.grading_vector,
                }
            }
        };
        if let Some(target_points) = &target_points {
            check_target_points(target_points, &grading_vector)
                .map_err(|reason| IoError::invalid_field("target_points", reason))?;
        }
        let unsupported = modular_unsupported(max_deg, min_points, target_points.is_some());
        if let (Method::Modular, Some(what)) = (method, unsupported) {
            return Err(IoError::invalid_field(
                "method",
                ComputeError::Unsupported(what).to_string(),
            ));
        }
        // The cone data describes the Mori cone without them, but without any truncation the
        // generators are the curve classes to use.
        let truncated = max_deg.is_some() || min_points.is_some() || target_points.is_some();
        let generators = match generators {
            Some(generators) => generators,
            None if (cone_data.is_some() || fan.is_some()) && truncated => DMatrix::zeros(h11, 0),
            None => return Err(IoError::MissingField("generators")),
        };
        // These are checked again, deeper in the computation, but the errors there are
        // not recoverable, so they are caught here to report them properly.
        let n_rays = q.nrows() as i64;
        if nefpart
            .iter()
            .flat_map(|p| p.iter())
            .any(|c| *c < 0 || (*c as i64) >= n_rays)
        {
            return Err(IoError::invalid_field(
                "nefpart",
                format!("the indices must be indices of the {n_rays} rays of \"q\""),
            ));
        }
        let cy_dim = n_rays - q.ncols() as i64 - nefpart.len().max(1) as i64;
        if cy_dim < 3 {
            return Err(IoError::invalid_field(
                "q",
                format!("the dimension of the CY must be at least three, but it is {cy_dim}"),
            ));
        }
        check_strata(&strata, q.nrows(), h11, cy_kind, nefpart.is_empty())
            .map_err(|reason| IoError::invalid_field("strata", reason))?;
        if let Some(cone_data) = &cone_data {
            check_cone_data(cone_data, &q, &generators, &strata)
                .map_err(|(field, reason)| IoError::invalid_field(field, reason))?;
        }
        // The intersection numbers are canonicalized differently depending on the
        // dimension of the CY, so they can only be checked once it is known.
        process_int_nums(intnums.clone(), cy_kind)
            .map_err(|e| IoError::invalid_field("intnums", e.to_string()))?;

        Ok(Input {
            name,
            generators,
            grading_vector,
            q,
            nefpart,
            intnums,
            strata,
            fan,
            cone_data,
            invariant_kind,
            cy_kind,
            max_deg,
            min_points,
            target_points,
            n_threads,
            method,
            certify,
        })
    }

    /// Run the computation specified by this input.
    ///
    /// The results are sorted by degree, and then by curve class and reference surface
    /// index. The last two are only used to break ties, since the order in which
    /// curve classes of equal degree are computed is not deterministic.
    pub fn compute(&self) -> Result<Vec<InvariantResult>, ComputeError> {
        let mut results = self.compute_unsorted()?;
        results.sort_unstable_by(|a, b| {
            let (curve_a, curve_b) = (&a.0 .0, &b.0 .0);
            self.degree(curve_a)
                .cmp(&self.degree(curve_b))
                .then_with(|| curve_a.as_slice().cmp(curve_b.as_slice()))
                .then_with(|| a.0 .1.cmp(&b.0 .1))
        });
        Ok(results)
    }

    /// The degree of a curve class under the grading vector.
    fn degree(&self, curve: &DVector<i32>) -> i64 {
        self.grading_vector
            .iter()
            .zip(curve.iter())
            .map(|(g, c)| (*g as i64) * (*c as i64))
            .sum()
    }

    /// Run the computation specified by this input, leaving the results in the order
    /// in which they were computed.
    fn compute_unsorted(&self) -> Result<Vec<InvariantResult>, ComputeError> {
        compute_gvgw_strings(
            self.generators.clone(),
            self.grading_vector.clone(),
            self.q.clone(),
            self.nefpart.clone(),
            self.intnums.clone(),
            self.strata.clone(),
            self.fan.clone(),
            self.cone_data.clone(),
            self.invariant_kind,
            self.cy_kind,
            self.max_deg,
            self.min_points,
            self.target_points.clone(),
            self.n_threads,
            self.method,
            self.certify,
        )
    }

    /// Write the results of a computation as a YAML document.
    ///
    /// The degree of each curve class under the grading vector is written along with
    /// it, since it is what the results are sorted by. GV invariants are written as
    /// integers, while GW invariants are written as quoted strings, since YAML has no
    /// way to represent fractions.
    pub fn write_results(
        &self,
        writer: &mut impl Write,
        results: &[InvariantResult],
    ) -> std::io::Result<()> {
        let key = self.invariant_kind.as_str();
        writeln!(writer, "---")?;
        if let Some(name) = &self.name {
            writeln!(writer, "name: {}", quote_string(name))?;
        }
        writeln!(writer, "invariants: {key}")?;
        writeln!(writer, "is_threefold: {}", self.cy_kind.is_threefold())?;
        write!(writer, "grading_vector: ")?;
        write_int_seq(writer, self.grading_vector.iter())?;
        writeln!(writer)?;
        if results.is_empty() {
            return writeln!(writer, "results: []");
        }
        writeln!(writer, "results:")?;
        for ((curve, idx), value) in results {
            write!(writer, "  - {{curve_class: ")?;
            write_int_seq(writer, curve.iter())?;
            write!(writer, ", degree: {}", self.degree(curve))?;
            if !self.cy_kind.is_threefold() {
                write!(writer, ", surface_index: {idx}")?;
            }
            match self.invariant_kind {
                InvariantKind::GV => writeln!(writer, ", {key}: {value}}}")?,
                InvariantKind::GW => writeln!(writer, ", {key}: '{value}'}}")?,
            }
        }
        Ok(())
    }
}

/// Deduce the kind of CY from the shape of the GLSM charge matrix and the nef
/// partition.
fn infer_cy_kind(q: &DMatrix<i32>, nefpart: &[DVector<i32>]) -> CYKind {
    let ambient_dim = q.nrows() as i64 - q.ncols() as i64;
    let cy_codim = if nefpart.is_empty() {
        1
    } else {
        nefpart.len() as i64
    };
    CYKind::from_is_threefold(ambient_dim - cy_codim == 3)
}

/// Look up a field, treating an explicitly null value as a missing one.
fn get_field<'a>(hash: &'a Hash, field: &str) -> Option<&'a Yaml> {
    match hash.get(&Yaml::String(field.to_owned())) {
        Some(y) if !y.is_null() && !y.is_badvalue() => Some(y),
        _ => None,
    }
}

/// Look up a field that must be present.
fn require_field<'a>(hash: &'a Hash, field: &'static str) -> Result<&'a Yaml, IoError> {
    get_field(hash, field).ok_or(IoError::MissingField(field))
}

/// Read a value that must be an integer.
fn as_i64(data: &Yaml, field: &str) -> Result<i64, IoError> {
    data.as_i64()
        .ok_or_else(|| IoError::invalid_field(field, "expected an integer"))
}

/// Read a value that must be a 32-bit integer.
fn as_i32(data: &Yaml, field: &str) -> Result<i32, IoError> {
    i32::try_from(as_i64(data, field)?)
        .map_err(|_| IoError::invalid_field(field, "the integers must fit in 32 bits"))
}

/// Read a value that must be a non-negative 32-bit integer.
fn as_u32(data: &Yaml, field: &str) -> Result<u32, IoError> {
    u32::try_from(as_i64(data, field)?)
        .map_err(|_| IoError::invalid_field(field, "expected a non-negative integer"))
}

/// Read a value that must be a list of integers.
fn as_int_vec(data: &Yaml, field: &str) -> Result<Vec<i32>, IoError> {
    let entries = data
        .as_vec()
        .ok_or_else(|| IoError::invalid_field(field, "expected a list of integers"))?;
    entries.iter().map(|c| as_i32(c, field)).collect()
}

/// Read a value that must be a list of vectors, which are allowed to have different
/// lengths.
fn as_vectors(data: &Yaml, field: &str) -> Result<Vec<DVector<i32>>, IoError> {
    let vectors = data
        .as_vec()
        .ok_or_else(|| IoError::invalid_field(field, "expected a list of vectors"))?;
    vectors
        .iter()
        .map(|v| Ok(DVector::from_vec(as_int_vec(v, field)?)))
        .collect()
}

/// Read a value that must be a list of vectors of equal length, and turn it into a
/// matrix whose columns are the given vectors.
fn as_matrix(data: &Yaml, field: &str) -> Result<DMatrix<i32>, IoError> {
    let columns = data
        .as_vec()
        .ok_or_else(|| IoError::invalid_field(field, "expected a list of vectors"))?;
    let columns: Vec<_> = columns
        .iter()
        .map(|c| as_int_vec(c, field))
        .collect::<Result<_, _>>()?;
    let n_cols = columns.len();
    let n_rows = columns.first().map_or(0, Vec::len);
    if n_cols == 0 || n_rows == 0 {
        return Err(IoError::invalid_field(
            field,
            "expected a non-empty list of non-empty vectors",
        ));
    }
    if columns.iter().any(|c| c.len() != n_rows) {
        return Err(IoError::invalid_field(
            field,
            "all of the vectors must have the same length",
        ));
    }
    Ok(DMatrix::from_iterator(
        n_rows,
        n_cols,
        columns.into_iter().flatten(),
    ))
}

/// Read the maximal cones of a fan, a list of lists of indices.
fn as_fan(data: &Yaml) -> Result<Vec<Vec<usize>>, IoError> {
    let err = || IoError::invalid_field("fan", "expected a list of lists of indices of rows of q");
    data.as_vec()
        .ok_or_else(err)?
        .iter()
        .map(|c| {
            c.as_vec()
                .ok_or_else(err)?
                .iter()
                .map(|i| {
                    i.as_i64()
                        .and_then(|i| usize::try_from(i).ok())
                        .ok_or_else(err)
                })
                .collect()
        })
        .collect()
}

/// Read the strata of a vex phase, a list of mappings with the three `divisors` and the
/// `class`, whose entries are integers or fractions written as strings.
fn as_strata(data: &Yaml) -> Result<Vec<Stratum>, IoError> {
    let field = "strata";
    let entries = data.as_vec().ok_or_else(|| {
        IoError::invalid_field(field, "expected a list of mappings with divisors and class")
    })?;
    entries
        .iter()
        .map(|e| {
            let hash = e.as_hash().ok_or_else(|| {
                IoError::invalid_field(
                    field,
                    "each stratum must be a mapping with divisors and class",
                )
            })?;
            if hash
                .keys()
                .any(|k| !matches!(k.as_str(), Some("divisors" | "class")))
            {
                return Err(IoError::invalid_field(
                    field,
                    "the fields of a stratum are divisors and class",
                ));
            }
            let divisors = as_int_vec(require_field(hash, "divisors")?, field)?;
            let divisors: [usize; 3] = divisors
                .iter()
                .map(|&d| usize::try_from(d).ok())
                .collect::<Option<Vec<_>>>()
                .and_then(|d| d.try_into().ok())
                .ok_or_else(|| {
                    IoError::invalid_field(field, "the divisors must be three non-negative indices")
                })?;
            let class = require_field(hash, "class")?
                .as_vec()
                .ok_or_else(|| IoError::invalid_field(field, "the class must be a list"))?
                .iter()
                .map(as_fraction)
                .collect::<Result<_, _>>()?;
            Ok(Stratum { divisors, class })
        })
        .collect()
}

/// Read the cones of the cone data, a list of mappings with the `divisors` and the
/// `hilbert_basis`, which may be empty.
fn as_cones(data: &Yaml, h11: usize) -> Result<Vec<CandidateCone>, IoError> {
    let field = "cones";
    let entries = data.as_vec().ok_or_else(|| {
        IoError::invalid_field(
            field,
            "expected a list of mappings with divisors and hilbert_basis",
        )
    })?;
    entries
        .iter()
        .map(|e| {
            let hash = e.as_hash().ok_or_else(|| {
                IoError::invalid_field(
                    field,
                    "each cone must be a mapping with divisors and hilbert_basis",
                )
            })?;
            if hash
                .keys()
                .any(|k| !matches!(k.as_str(), Some("divisors" | "hilbert_basis")))
            {
                return Err(IoError::invalid_field(
                    field,
                    "the fields of a cone are divisors and hilbert_basis",
                ));
            }
            let divisors = as_int_vec(require_field(hash, "divisors")?, field)?
                .into_iter()
                .map(|d| usize::try_from(d).ok())
                .collect::<Option<_>>()
                .ok_or_else(|| {
                    IoError::invalid_field(field, "the divisors must be non-negative indices")
                })?;
            let basis = require_field(hash, "hilbert_basis")?;
            let hilbert_basis = match basis.as_vec() {
                Some(v) if v.is_empty() => DMatrix::zeros(h11, 0),
                _ => as_matrix(basis, field)?,
            };
            Ok(CandidateCone {
                divisors,
                hilbert_basis,
            })
        })
        .collect()
}

/// Read a fraction, given as an integer or as a string of the form `"p/q"` or `"p"`.
fn as_fraction(data: &Yaml) -> Result<(i64, i64), IoError> {
    let invalid = || {
        IoError::invalid_field(
            "strata",
            "the entries of a class must be integers or fractions such as \"-1/2\"",
        )
    };
    match data {
        Yaml::Integer(n) => Ok((*n, 1)),
        Yaml::String(s) => {
            let (num, den) = s.split_once('/').unwrap_or((s, "1"));
            let num = num.trim().parse().map_err(|_| invalid())?;
            let den = den.trim().parse().map_err(|_| invalid())?;
            Ok((num, den))
        }
        _ => Err(invalid()),
    }
}

/// Read the target points, which can be given either as a single vector or as a list
/// of vectors.
fn as_target_points(data: &Yaml, h11: usize) -> Result<DMatrix<i32>, IoError> {
    let field = "target_points";
    let is_single_point = matches!(
        data.as_vec().and_then(|v| v.first()),
        Some(Yaml::Integer(_))
    );
    let points = if is_single_point {
        let point = as_int_vec(data, field)?;
        DMatrix::from_column_slice(point.len(), 1, &point)
    } else {
        as_matrix(data, field)?
    };
    if points.nrows() != h11 {
        return Err(IoError::invalid_field(
            field,
            format!("expected vectors with {h11} entries to match the grading vector"),
        ));
    }
    Ok(points)
}

/// Read the intersection numbers, which can be given either as a list of
/// `[i, j, k, value]` entries, or as a mapping from `"i,j,k"` keys to values.
fn as_intnums(data: &Yaml) -> Result<HashMap<(usize, usize, usize), i32>, IoError> {
    let field = "intnums";
    let mut intnums = HashMap::new();
    let entries: Vec<((usize, usize, usize), i32)> = match data {
        Yaml::Array(entries) => entries
            .iter()
            .map(|e| {
                let entry = as_int_vec(e, field)?;
                if entry.len() != 4 {
                    return Err(IoError::invalid_field(
                        field,
                        "each entry must be of the form [i, j, k, value]",
                    ));
                }
                Ok((as_indices(&entry[..3])?, entry[3]))
            })
            .collect::<Result<_, _>>()?,
        Yaml::Hash(entries) => entries
            .iter()
            .map(|(k, v)| {
                let indices = match k {
                    Yaml::String(s) => s
                        .split([',', ' '])
                        .filter(|c| !c.is_empty())
                        .map(|c| c.parse::<i32>().ok())
                        .collect::<Option<Vec<_>>>()
                        .ok_or_else(|| {
                            IoError::invalid_field(
                                field,
                                "the keys must be indices of the form \"i,j,k\"",
                            )
                        })?,
                    Yaml::Array(_) => as_int_vec(k, field)?,
                    _ => {
                        return Err(IoError::invalid_field(
                            field,
                            "the keys must be indices of the form \"i,j,k\"",
                        ))
                    }
                };
                if indices.len() != 3 {
                    return Err(IoError::invalid_field(
                        field,
                        "the keys must be indices of the form \"i,j,k\"",
                    ));
                }
                Ok((as_indices(&indices)?, as_i32(v, field)?))
            })
            .collect::<Result<_, _>>()?,
        _ => {
            return Err(IoError::invalid_field(
                field,
                "expected a list of [i, j, k, value] entries",
            ))
        }
    };
    for (indices, value) in entries {
        if intnums.insert(indices, value).is_some() {
            return Err(IoError::invalid_field(
                field,
                format!(
                    "the indices {},{},{} appear more than once",
                    indices.0, indices.1, indices.2
                ),
            ));
        }
    }
    if intnums.is_empty() {
        return Err(IoError::invalid_field(
            field,
            "there must be at least one intersection number",
        ));
    }
    Ok(intnums)
}

/// Turn a triplet of intersection number indices into a tuple.
fn as_indices(indices: &[i32]) -> Result<(usize, usize, usize), IoError> {
    if indices.iter().any(|c| *c < 0) {
        return Err(IoError::invalid_field(
            "intnums",
            "the indices must be non-negative",
        ));
    }
    Ok((
        indices[0] as usize,
        indices[1] as usize,
        indices[2] as usize,
    ))
}

/// Write a sequence of integers as a YAML flow sequence.
fn write_int_seq<'a>(
    writer: &mut impl Write,
    values: impl Iterator<Item = &'a i32>,
) -> std::io::Result<()> {
    write!(writer, "[")?;
    for (i, c) in values.enumerate() {
        if i > 0 {
            write!(writer, ", ")?;
        }
        write!(writer, "{c}")?;
    }
    write!(writer, "]")
}

/// Turn a string into a YAML scalar, quoting it only when it is needed.
fn quote_string(data: &str) -> String {
    let is_plain = data.starts_with(|c: char| c.is_ascii_alphabetic())
        && !data.ends_with(' ')
        && data
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || " _-.".contains(c))
        && !matches!(
            data.to_lowercase().as_str(),
            "true" | "false" | "null" | "yes" | "no" | "on" | "off"
        );
    if is_plain {
        return data.to_owned();
    }
    let mut res = String::with_capacity(data.len() + 2);
    res.push('"');
    for c in data.chars() {
        match c {
            '"' => res.push_str("\\\""),
            '\\' => res.push_str("\\\\"),
            '\n' => res.push_str("\\n"),
            '\r' => res.push_str("\\r"),
            '\t' => res.push_str("\\t"),
            _ => res.push(c),
        }
    }
    res.push('"');
    res
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{dmatrix, dvector};

    const THREEFOLD: &str = "\
name: test
generators: [[0, -1], [1, 2]]
grading_vector: [3, -1]
q: [[1, 1, 1, 0, 1, 2], [0, 0, -1, 1, 1, -1]]
intnums: [[0, 0, 0, 2], [0, 0, 1, 1], [0, 1, 1, -1], [1, 1, 1, 5]]
min_points: 20
";

    fn load_one(data: &str) -> Result<Input, IoError> {
        Input::load_all(data).map(|mut v| v.remove(0))
    }

    #[test]
    fn test_parse_threefold() {
        let input = load_one(THREEFOLD).unwrap();
        assert_eq!(input.name.as_deref(), Some("test"));
        assert_eq!(input.generators, dmatrix![0, 1; -1, 2]);
        assert_eq!(input.grading_vector, RowDVector::from_row_slice(&[3, -1]));
        assert_eq!(input.q.shape(), (6, 2));
        assert_eq!(input.q[(2, 1)], -1);
        assert_eq!(input.intnums[&(0, 1, 1)], -1);
        assert!(input.nefpart.is_empty());
        assert_eq!(input.invariant_kind, InvariantKind::GV);
        assert_eq!(input.cy_kind, CYKind::Threefold);
        assert_eq!(input.min_points, Some(20));
        assert_eq!(input.max_deg, None);
        assert_eq!(input.target_points, None);
        assert_eq!(input.method, Method::default());
    }

    #[test]
    fn test_parse_strata() {
        let base = THREEFOLD.replace("min_points: 20\n", "max_deg: 5\n");
        let input = load_one(&format!(
            "{base}strata:\n  - {{divisors: [3, 1, 2], class: [1, \"-1/2\"]}}\n  \
             - {{divisors: [0, 1, 5], class: [\"3\", 0]}}\n"
        ))
        .unwrap();
        assert_eq!(
            input.strata,
            [
                Stratum {
                    divisors: [3, 1, 2],
                    class: vec![(1, 1), (-1, 2)]
                },
                Stratum {
                    divisors: [0, 1, 5],
                    class: vec![(3, 1), (0, 1)]
                }
            ]
        );
        assert!(load_one(&base).unwrap().strata.is_empty());

        for bad in [
            "[{divisors: [1, 2], class: [1, 1]}]",
            "[{divisors: [1, 2, 6], class: [1, 1]}]",
            "[{divisors: [1, 2, 2], class: [1, 1]}]",
            "[{divisors: [1, 2, 3], class: [1]}]",
            "[{divisors: [1, 2, 3], class: [1, \"1/0\"]}]",
            "[{divisors: [1, 2, 3], class: [1, \"a\"]}]",
            "[{divisors: [1, 2, 3], class: [1, 1], extra: 1}]",
            "[{divisors: [1, 2, 3], class: [1, 1]}, {divisors: [3, 2, 1], class: [0, 1]}]",
            "[[1, 2, 3], [1, 1]]",
        ] {
            assert!(
                matches!(
                    load_one(&format!("{base}strata: {bad}\n")),
                    Err(IoError::InvalidField { field, .. }) if field == "strata"
                ),
                "{bad}"
            );
        }
        assert_eq!(
            load_one(&format!("{base}strata: [{{divisors: [1, 2, 3]}}]\n")),
            Err(IoError::MissingField("class"))
        );
        // Strata are only supported for threefold hypersurfaces, by both methods.
        let strata = "strata: [{divisors: [1, 2, 3], class: [1, 1]}]\n";
        assert!(matches!(
            load_one(&format!("{base}{strata}is_threefold: false\n")),
            Err(IoError::InvalidField { field, .. }) if field == "strata"
        ));
        assert!(load_one(&format!("{base}{strata}method: modular\n")).is_ok());
    }

    /// A vex phase of P(1,1,2,2,2)[8] given by its fan.
    const FAN: &str = "\
q: [[1, 0, -1, 1, 0, 1], [0, 1, 1, 1, 1, 0]]
intnums: [[0, 1, 1, 4], [1, 1, 1, 2]]
max_deg: 6
fan: [[0, 1, 2, 3], [0, 1, 2, 4], [0, 1, 3, 4], [0, 2, 3, 4], [1, 2, 3, 5], [1, 2, 4, 5],
  [1, 3, 4, 5], [2, 3, 4, 5]]
";

    #[test]
    fn test_parse_fan() {
        // The grading vector is derived from the fan; the strata and the generators are
        // left to the computation.
        let input = load_one(FAN).unwrap();
        assert_eq!(input.grading_vector, RowDVector::from_row_slice(&[2, 1]));
        assert_eq!(input.generators.shape(), (2, 0));
        assert!(input.strata.is_empty());
        assert_eq!(input.fan.as_ref().map(Vec::len), Some(8));
        let given = load_one(&format!("{FAN}grading_vector: [3, 1]\n")).unwrap();
        assert_eq!(given.grading_vector, RowDVector::from_row_slice(&[3, 1]));

        let field_of = |data: &str| match load_one(data) {
            Err(IoError::InvalidField { field, .. }) => field,
            res => panic!("{res:?}"),
        };
        let strata = "strata: [{divisors: [1, 2, 4], class: [1, -1]}]\n";
        assert_eq!(field_of(&format!("{FAN}{strata}")), "strata");
        assert_eq!(field_of(&format!("{FAN}is_threefold: false\n")), "fan");
        assert_eq!(
            field_of(&FAN.replace("[1, 1, 1, 2]", "[1, 1, 1, 3]")),
            "fan"
        );
        assert_eq!(
            field_of(&FAN.replace("[1, 1, 1, 2]", "[1, 1, 2, 2]")),
            "intnums"
        );
        assert_eq!(field_of(&FAN.replace("[0, 1, 2, 3], ", "")), "fan");
        assert_eq!(
            field_of(&FAN.replace("[0, 1, 2, 3]", "[0, 1, 2, -3]")),
            "fan"
        );
        assert_eq!(field_of(&format!("{FAN}generators: [[1, 0, 0]]\n")), "q");
        // Without any truncation the generators are the curve classes to use.
        let untruncated = FAN.replace("max_deg: 6\n", "method: reference\n");
        assert_eq!(
            load_one(&untruncated),
            Err(IoError::MissingField("generators"))
        );
        assert!(load_one(&format!("{untruncated}generators: [[0, 1], [1, -1]]\n")).is_ok());
        // Without a fan, the grading vector is needed.
        let without = FAN.replace("fan: ", "strata: []\nignored: ");
        let without = &without[..without.find("ignored").unwrap()];
        assert_eq!(
            load_one(without),
            Err(IoError::MissingField("grading_vector"))
        );
    }

    /// P(1,1,2,2,2)[8] with its cone data and no generators.
    fn p11222_cones() -> String {
        let mut data = "\
grading_vector: [1, 1]
q: [[0, 0, 1, 1, 1, 1], [1, 1, 0, 0, 0, -2]]
intnums: [[0, 0, 0, 8], [0, 0, 1, 4]]
max_deg: 6
mori_hyperplanes: [[1, 0], [0, 1]]
cones:
"
        .to_owned();
        for s in 0..6 {
            for t in s..6 {
                let divisors = if s == t { vec![s] } else { vec![s, t] };
                let basis = if divisors.contains(&5) {
                    "[[1, 0], [0, 1]]"
                } else {
                    "[[1, 0], [2, 1]]"
                };
                data += &format!("  - {{divisors: {divisors:?}, hilbert_basis: {basis}}}\n");
            }
        }
        data
    }

    #[test]
    fn test_parse_cone_data() {
        let data = p11222_cones();
        let input = load_one(&data).unwrap();
        assert_eq!(input.generators.shape(), (2, 0));
        let cone_data = input.cone_data.unwrap();
        assert_eq!(cone_data.hyperplanes, dmatrix![1, 0; 0, 1]);
        assert_eq!(cone_data.equations.shape(), (2, 0));
        assert_eq!(cone_data.cones.len(), 21);
        assert_eq!(cone_data.cones[1].divisors, [0, 1]);
        assert_eq!(cone_data.cones[1].hilbert_basis, dmatrix![1, 2; 0, 1]);
        assert!(load_one(&format!("{data}generators: [[1, 0], [0, 1]]\n"))
            .unwrap()
            .cone_data
            .is_some());

        // The face spanned by (1, 0), which meets no divisor negatively, so it is every K_T.
        let face = data
            .replace("max_deg: 6\n", "max_deg: 6\nmori_equations: [[0, 1]]\n")
            .replace("[[1, 0], [0, 1]]}", "[[1, 0]]}")
            .replace("[[1, 0], [2, 1]]}", "[[1, 0]]}");
        let cone_data = load_one(&face).unwrap().cone_data.unwrap();
        assert_eq!(cone_data.equations, dmatrix![0; 1]);
        assert_eq!(cone_data.cones[0].hilbert_basis, dmatrix![1; 0]);
        // An empty Hilbert basis is read too (a K_T can be just the origin).
        let empty = data.replacen("[[1, 0], [2, 1]]}", "[]}", 1);
        let cone_data = load_one(&empty).unwrap().cone_data.unwrap();
        assert_eq!(cone_data.cones[0].hilbert_basis.shape(), (2, 0));

        // The cone data needs both the hyperplanes and the cones, and replaces the generators
        // for every method, unless nothing truncates the computation.
        let no_cones = &data[..data.find("cones:").unwrap()];
        assert_eq!(load_one(no_cones), Err(IoError::MissingField("cones")));
        let no_hyperplanes = data.replace("mori_hyperplanes: [[1, 0], [0, 1]]\n", "");
        assert_eq!(
            load_one(&no_hyperplanes),
            Err(IoError::MissingField("mori_hyperplanes"))
        );
        let head = &data[..data.find("mori_hyperplanes").unwrap()];
        assert_eq!(
            load_one(&format!("{head}mori_equations: [[0, 1]]\n")),
            Err(IoError::MissingField("mori_hyperplanes"))
        );
        assert_eq!(load_one(head), Err(IoError::MissingField("generators")));
        for extra in ["method: reference\nmax_deg: 6\n", "min_points: 5\n"] {
            let data = data.replace("max_deg: 6\n", extra);
            assert_eq!(load_one(&data).unwrap().generators.shape(), (2, 0));
        }
        let untruncated = data.replace("max_deg: 6\n", "method: reference\n");
        assert_eq!(
            load_one(&untruncated),
            Err(IoError::MissingField("generators"))
        );

        // The cone data is checked.
        let field_of = |data: &str| match load_one(data) {
            Err(IoError::InvalidField { field, .. }) => field,
            res => panic!("{res:?}"),
        };
        let last = data.rfind("  - ").unwrap();
        assert_eq!(field_of(&data[..last]), "cones");
        for (from, to) in [
            (
                "{divisors: [0], hilbert_basis: [[1, 0], [2, 1]]}",
                "{divisors: [0], hilbert_basis: [[0, 1]]}",
            ),
            (
                "{divisors: [0], hilbert_basis: [[1, 0], [2, 1]]}",
                "{divisors: [-1], hilbert_basis: []}",
            ),
            (
                "{divisors: [0], hilbert_basis: [[1, 0], [2, 1]]}",
                "{divisors: [0], basis: []}",
            ),
            ("{divisors: [0], hilbert_basis: [[1, 0], [2, 1]]}", "[0]"),
            // Three divisors only for a stratum.
            (
                "{divisors: [0], hilbert_basis: [[1, 0], [2, 1]]}",
                "{divisors: [0, 1, 2], hilbert_basis: []}",
            ),
        ] {
            assert_eq!(field_of(&data.replacen(from, to, 1)), "cones", "{to}");
        }
        assert_eq!(
            field_of(&data.replace(
                "mori_hyperplanes: [[1, 0], [0, 1]]",
                "mori_hyperplanes: [[1, 0, 0]]"
            )),
            "mori_hyperplanes"
        );
        assert_eq!(
            field_of(&format!("{data}generators: [[1, 0], [-1, 1]]\n")),
            "generators"
        );
        assert_eq!(
            field_of(&format!("{data}mori_equations: [[1]]\n")),
            "mori_equations"
        );
    }

    #[test]
    fn test_parse_method() {
        let base = THREEFOLD.replace("min_points: 20\n", "max_deg: 5\n");
        for (name, method) in [
            ("auto", Method::Auto),
            ("Modular", Method::Modular),
            ("reference", Method::Reference),
        ] {
            let input = load_one(&format!("{base}method: {name}\n")).unwrap();
            assert_eq!(input.method, method);
        }
        assert!(matches!(
            load_one(&format!("{base}method: fast\n")),
            Err(IoError::InvalidField { field, .. }) if field == "method"
        ));
        // The modular method does not support min_points, which is checked eagerly.
        assert!(load_one(&format!("{THREEFOLD}method: auto\n")).is_ok());
        assert!(matches!(
            load_one(&format!("{THREEFOLD}method: modular\n")),
            Err(IoError::InvalidField { field, reason })
                if field == "method" && reason.contains("min_points")
        ));
    }

    #[test]
    fn test_parse_certify() {
        let base = THREEFOLD.replace("min_points: 20\n", "max_deg: 5\n");
        assert!(!load_one(&base).unwrap().certify);
        assert!(load_one(&format!("{base}certify: true\n")).unwrap().certify);
        assert!(matches!(
            load_one(&format!("{base}certify: 1\n")),
            Err(IoError::InvalidField { field, .. }) if field == "certify"
        ));
    }

    #[test]
    fn test_parse_multiple_documents() {
        let data = format!("---\n{THREEFOLD}---\n{THREEFOLD}invariants: gw\n");
        let inputs = Input::load_all(&data).unwrap();
        assert_eq!(inputs.len(), 2);
        assert_eq!(inputs[0].invariant_kind, InvariantKind::GV);
        assert_eq!(inputs[1].invariant_kind, InvariantKind::GW);
    }

    #[test]
    fn test_parse_intnums_mapping() {
        let data = THREEFOLD.replace(
            "intnums: [[0, 0, 0, 2], [0, 0, 1, 1], [0, 1, 1, -1], [1, 1, 1, 5]]",
            "intnums: {\"0,0,0\": 2, \"0 0 1\": 1, \"0,1,1\": -1, \"1,1,1\": 5}",
        );
        let input = load_one(&data).unwrap();
        assert_eq!(input.intnums, load_one(THREEFOLD).unwrap().intnums);
    }

    #[test]
    fn test_parse_target_points() {
        let base = THREEFOLD.replace("min_points: 20\n", "");
        let input = load_one(&format!("{base}target_points: [1, 2]\n")).unwrap();
        assert_eq!(input.target_points, Some(dmatrix![1; 2]));

        let input = load_one(&format!("{base}target_points: [[1, 2], [3, 4]]\n")).unwrap();
        assert_eq!(input.target_points, Some(dmatrix![1, 3; 2, 4]));

        // With a maximum degree, which the modular method supports too.
        let input = load_one(&format!(
            "{base}target_points: [1, 2]\nmax_deg: 3\nmethod: modular\n"
        ))
        .unwrap();
        assert_eq!(input.max_deg, Some(3));

        assert_eq!(
            load_one(&format!("{THREEFOLD}target_points: [1, 2]\n")),
            Err(IoError::ConflictingTruncations)
        );
        assert_eq!(
            load_one(&format!(
                "{}target_points: [1, 2]\n",
                THREEFOLD.replace("min_points: 20\n", "max_deg: 3\n")
            ))
            .unwrap()
            .max_deg,
            Some(3)
        );
        for bad in ["[[1, 2, 3]]", "[0, 0]", "[[1, 2], [-1, -1]]"] {
            assert!(matches!(
                load_one(&format!("{base}target_points: {bad}\n")),
                Err(IoError::InvalidField { field, .. }) if field == "target_points"
            ));
        }
    }

    #[test]
    fn test_infer_cy_kind() {
        assert_eq!(load_one(THREEFOLD).unwrap().cy_kind, CYKind::Threefold);

        // Eight rays, two Kähler parameters, and a codimension-two CY, i.e. a fourfold.
        let data = THREEFOLD.replace(
            "q: [[1, 1, 1, 0, 1, 2], [0, 0, -1, 1, 1, -1]]",
            "q: [[1, 1, 1, 0, 1, 2, 1, 0], [0, 0, -1, 1, 1, -1, 0, 1]]",
        );
        let data = format!("{data}nefpart: [[0, 1], [2, 3, 4, 5, 6, 7]]\n");
        let input = load_one(&data).unwrap();
        assert_eq!(input.nefpart.len(), 2);
        assert_eq!(input.cy_kind, CYKind::Nfold);

        assert_eq!(
            load_one(&format!("{data}is_threefold: true\n"))
                .unwrap()
                .cy_kind,
            CYKind::Threefold
        );
        // The nef partition must consist of indices of the rays.
        assert!(matches!(
            load_one(&data.replace("[2, 3, 4, 5, 6, 7]", "[2, 3, 4, 5, 6, 8]")),
            Err(IoError::InvalidField { field, .. }) if field == "nefpart"
        ));
        // The resulting CY must be at least three-dimensional.
        assert!(matches!(
            load_one(&format!("{THREEFOLD}nefpart: [[0, 1], [2, 3, 4, 5]]\n")),
            Err(IoError::InvalidField { field, .. }) if field == "q"
        ));
    }

    #[test]
    fn test_parse_errors() {
        assert_eq!(Input::load_all(""), Err(IoError::EmptyInput));
        assert_eq!(Input::load_all("- 1"), Err(IoError::NotAMapping));
        assert_eq!(
            load_one("generators: [[1]]"),
            Err(IoError::MissingField("grading_vector"))
        );
        assert_eq!(
            load_one(&format!("{THREEFOLD}bad_field: 1\n")),
            Err(IoError::UnknownField("bad_field".to_owned()))
        );
        assert!(matches!(
            load_one("generators: [[1, 2], [3]]"),
            Err(IoError::InvalidField { field, .. }) if field == "generators"
        ));
        assert!(matches!(
            load_one(&THREEFOLD.replace("[3, -1]", "[3, -1, 0]")),
            Err(IoError::InvalidField { field, .. }) if field == "grading_vector"
        ));
        assert!(matches!(
            load_one(&THREEFOLD.replace("[0, 1, 1, -1]", "[0, 1, 1, -1], [1, 0, 1, 2]")),
            Err(IoError::InvalidField { field, .. }) if field == "intnums"
        ));
        // Intersection numbers that are all zero must be rejected, not panic.
        assert!(matches!(
            load_one(&THREEFOLD.replace(
                "intnums: [[0, 0, 0, 2], [0, 0, 1, 1], [0, 1, 1, -1], [1, 1, 1, 5]]",
                "intnums: [[0, 0, 0, 0], [0, 0, 1, 0]]",
            )),
            Err(IoError::InvalidField { field, .. }) if field == "intnums"
        ));
        assert!(Input::load_all("generators: [[1]], bad").is_err());
    }

    #[test]
    fn test_write_results() {
        let mut input = load_one(THREEFOLD).unwrap();
        let results = vec![
            ((dvector![0, -1], 0), "-4".to_owned()),
            ((dvector![1, 0], 2), "7524".to_owned()),
        ];

        let mut buf = Vec::new();
        input.write_results(&mut buf, &results).unwrap();
        assert_eq!(
            String::from_utf8(buf).unwrap(),
            "---\nname: test\ninvariants: gv\nis_threefold: true\n\
             grading_vector: [3, -1]\nresults:\n  \
             - {curve_class: [0, -1], degree: 1, gv: -4}\n  \
             - {curve_class: [1, 0], degree: 3, gv: 7524}\n"
        );

        input.name = None;
        input.invariant_kind = InvariantKind::GW;
        input.cy_kind = CYKind::Nfold;
        let mut buf = Vec::new();
        input.write_results(&mut buf, &results).unwrap();
        assert_eq!(
            String::from_utf8(buf).unwrap(),
            "---\ninvariants: gw\nis_threefold: false\n\
             grading_vector: [3, -1]\nresults:\n  \
             - {curve_class: [0, -1], degree: 1, surface_index: 0, gw: '-4'}\n  \
             - {curve_class: [1, 0], degree: 3, surface_index: 2, gw: '7524'}\n"
        );

        let mut buf = Vec::new();
        input.write_results(&mut buf, &[]).unwrap();
        assert_eq!(
            String::from_utf8(buf).unwrap(),
            "---\ninvariants: gw\nis_threefold: false\ngrading_vector: [3, -1]\nresults: []\n"
        );
    }

    #[test]
    fn test_quote_string() {
        assert_eq!(quote_string("a name"), "a name");
        assert_eq!(quote_string("CY-3.1_a"), "CY-3.1_a");
        assert_eq!(quote_string("yes"), "\"yes\"");
        assert_eq!(quote_string("2 fold"), "\"2 fold\"");
        assert_eq!(quote_string("a: b"), "\"a: b\"");
        assert_eq!(quote_string("a\"b\\c\nd"), "\"a\\\"b\\\\c\\nd\"");
        // Every name must round-trip through the emitted document.
        for name in ["a name", "yes", "2 fold", "a: b", "a\"b\\c\nd", "", "#x"] {
            let data = format!("name: {}\n", quote_string(name));
            let doc = &YamlLoader::load_from_str(&data).unwrap()[0];
            assert_eq!(doc["name"].as_str(), Some(name), "{data}");
        }
    }

    #[test]
    fn test_results_are_sorted() {
        let input = load_one(THREEFOLD).unwrap();
        let results = input.compute().unwrap();
        assert!(results.len() > 1);
        let degrees: Vec<_> = results.iter().map(|((c, _), _)| input.degree(c)).collect();
        assert!(degrees.windows(2).all(|w| w[0] <= w[1]));
        // The curve classes of equal degree are sorted lexicographically.
        assert!(results
            .windows(2)
            .all(|w| input.degree(&w[0].0 .0) < input.degree(&w[1].0 .0)
                || w[0].0 .0.as_slice() <= w[1].0 .0.as_slice()));
    }

    #[test]
    fn test_output_is_valid_yaml() {
        let input = load_one(THREEFOLD).unwrap();
        let results = input.compute().unwrap();
        assert!(!results.is_empty());

        let mut buf = Vec::new();
        input.write_results(&mut buf, &results).unwrap();
        let output = YamlLoader::load_from_str(&String::from_utf8(buf).unwrap()).unwrap();
        assert_eq!(output.len(), 1);
        assert_eq!(output[0]["name"].as_str(), Some("test"));
        assert_eq!(output[0]["invariants"].as_str(), Some("gv"));
        assert_eq!(output[0]["is_threefold"].as_bool(), Some(true));
        assert_eq!(
            as_int_vec(&output[0]["grading_vector"], "").unwrap(),
            input.grading_vector.as_slice()
        );
        let parsed = output[0]["results"].as_vec().unwrap();
        assert_eq!(parsed.len(), results.len());
        for (entry, ((curve, _), value)) in parsed.iter().zip(results.iter()) {
            assert_eq!(
                as_int_vec(&entry["curve_class"], "").unwrap(),
                curve.as_slice()
            );
            assert_eq!(entry["degree"].as_i64(), Some(input.degree(curve)));
            assert_eq!(entry["gv"].as_i64().unwrap().to_string(), *value);
        }
    }
}
