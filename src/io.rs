//! Input and output functions
//!
//! This module implements the YAML-based data format used by the `cygv` command line
//! interface. An input stream contains one or more YAML documents, each of which
//! specifies a single computation. The fields have the names of the fields of [`Problem`] and
//! [`Options`], and of the keywords of the Python package.
//!
//! ```yaml
//! ---
//! # Optional label that is echoed back in the output.
//! name: quintic-like threefold
//! # The GLSM charges, one charge vector per Kähler parameter.
//! q: [[1, 1, 1, 0, 1, 2], [0, 0, -1, 1, 1, -1]]
//! # The intersection numbers, given as a list of [i, j, k, value] entries.
//! intnums:
//!   - [0, 0, 0, 2]
//!   - [0, 0, 1, 1]
//!   - [0, 1, 1, -1]
//!   - [1, 1, 1, 5]
//! # The Mori cone: vectors that span it (or mori_hyperplanes, or nothing with a fan).
//! generators: [[0, -1], [1, 2]]
//! # The degree is the inner product with this vector, positive on the Mori cone. It
//! # defaults to the primitive sum of the inequalities of the Mori cone.
//! grading_vector: [3, -1]
//! # Which curve classes: those up to a maximum degree (max_deg), those below target
//! # points (target_points, which max_deg may accompany), at least a number of them
//! # (min_points), or exactly the ones listed (curve_classes).
//! max_deg: 10
//! # Whether to compute GV or GW invariants. Defaults to gv.
//! invariants: gv
//! # Which pipeline to use: auto, modular or reference. Defaults to auto, which uses the
//! # modular one, which is much faster, for max_deg and target_points.
//! method: auto
//! # Whether to prove the results exact, taking more primes if needed (not in vex phases).
//! # Defaults to false.
//! certify: false
//! ```
//!
//! A vex phase of the ambient toric variety needs its strata (see [`Stratum`]), each given by
//! three divisors and a class whose entries may be fractions.
//!
//! ```yaml
//! strata:
//!   - divisors: [1, 2, 4]
//!     class: [1, "-1/2"]
//! ```
//!
//! They can instead be derived from the fan of the ambient toric variety (see
//! [`crate::phase::Fan`]), given by its maximal cones as indices of the divisors. The fan also
//! gives the Mori cone, so `generators` can be left out, and the default grading vector; the
//! intersection numbers and a grading vector that is given are checked against it.
//!
//! ```yaml
//! fan: [[0, 1, 2, 3], [0, 1, 2, 4], [0, 1, 3, 4], [0, 2, 3, 4], [1, 2, 3, 5], [1, 2, 4, 5],
//!   [1, 3, 4, 5], [2, 3, 4, 5]]
//! ```
//!
//! The Mori cone can also be given by its inequalities, and the equations of its span if it is
//! a face, and the modular pipeline can take the Hilbert basis of the cone $K_T$ of every set
//! $T$ of one or two divisors from an external program such as Normaliz (see
//! [`crate::cone::ConeData`]) instead of computing them.
//!
//! ```yaml
//! mori_hyperplanes: [[1, 0], [0, 1]]
//! hilbert_bases:
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
//! dimension: 3
//! grading_vector: [3, -1]
//! report: {method: modular, certified: false, grading_vector_derived: false}
//! results:
//!   - {curve_class: [0, -1], degree: 1, gv: -4}
//!   - {curve_class: [1, 0], degree: 3, gv: 7524}
//! ```
//!
//! See [`Input`] for the full list of fields.

pub mod error;

use crate::cone::CandidateCone;
use crate::{
    compute, Ambient, Computed, Error, Geometry, InvariantKind, Invariants, Method, MoriCone,
    Options, Problem, Selection, Stratum,
};
use error::IoError;
use malachite::Rational;
use nalgebra::DMatrix;
use std::collections::HashMap;
use std::io::Write;
use std::str::FromStr;
use yaml_rust2::yaml::{Hash, Yaml};
use yaml_rust2::YamlLoader;

/// The fields that are accepted in an input document.
const INPUT_FIELDS: [&str; 20] = [
    "name",
    "q",
    "intnums",
    "nefpart",
    "fan",
    "strata",
    "generators",
    "mori_hyperplanes",
    "mori_equations",
    "grading_vector",
    "max_deg",
    "target_points",
    "min_points",
    "curve_classes",
    "invariants",
    "method",
    "n_threads",
    "certify",
    "hilbert_bases",
    // Rejected below with what to write instead.
    "cones",
];

/// The specification of a single computation of GV or GW invariants.
#[derive(Debug, Clone, PartialEq)]
pub struct Input {
    /// An optional label that is echoed back in the output.
    pub name: Option<String>,
    /// What to compute.
    pub problem: Problem,
    /// Which invariants to compute.
    pub kind: InvariantKind,
    /// How to compute them.
    pub options: Options,
}

impl From<Error> for IoError {
    fn from(e: Error) -> Self {
        match e {
            Error::Input { field, reason } => IoError::invalid_field(field, reason),
            e => IoError::Unsupported(e.to_string()),
        }
    }
}

impl Input {
    /// Parse every YAML document of the given string into an [`Input`], and check each.
    pub fn load_all(data: &str) -> Result<Vec<Self>, IoError> {
        let docs =
            YamlLoader::load_from_str(data).map_err(|e| IoError::InvalidYaml(e.to_string()))?;
        if docs.is_empty() {
            return Err(IoError::EmptyInput);
        }
        docs.iter().map(Input::from_yaml).collect()
    }

    /// Parse a single YAML document into an [`Input`], and check it ([`Input::validate`]).
    pub fn from_yaml(doc: &Yaml) -> Result<Self, IoError> {
        let input = Self::parse(doc)?;
        input.validate()?;
        Ok(input)
    }

    /// Checks the problem and the options, as far as that can be done before computing.
    pub fn validate(&self) -> Result<(), IoError> {
        Ok(self.problem.validate(&self.options)?)
    }

    /// Parse a single YAML document into an [`Input`], checking the shapes of the fields but
    /// not how they fit together.
    fn parse(doc: &Yaml) -> Result<Self, IoError> {
        let hash = doc.as_hash().ok_or(IoError::NotAMapping)?;
        for key in hash.keys() {
            match key.as_str() {
                Some("cones") => {
                    return Err(IoError::invalid_field(
                        "cones",
                        "it is called hilbert_bases now, and the Mori cone it came with is \
                         given by mori_hyperplanes alone",
                    ))
                }
                Some("is_threefold") => {
                    return Err(IoError::invalid_field(
                        "is_threefold",
                        "the dimension is derived from q and nefpart, so leave it out",
                    ))
                }
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
        let q = as_matrix(require_field(hash, "q")?, "q")?;
        let h11 = q.ncols();
        let intnums = as_intnums(require_field(hash, "intnums")?)?;
        let nefpart = match get_field(hash, "nefpart") {
            None => Vec::new(),
            Some(y) => as_index_lists(y, "nefpart")?,
        };
        let ambient = match (get_field(hash, "fan"), get_field(hash, "strata")) {
            (Some(_), Some(_)) => return Err(IoError::Conflicting("fan", "strata")),
            (Some(fan), None) => Ambient::Fan(as_index_lists(fan, "fan")?),
            (None, Some(strata)) => Ambient::Strata(as_strata(strata)?),
            (None, None) => Ambient::Plain,
        };
        let cone = match (
            get_field(hash, "generators"),
            get_field(hash, "mori_hyperplanes"),
            get_field(hash, "mori_equations"),
        ) {
            (Some(_), Some(_), _) => {
                return Err(IoError::Conflicting("generators", "mori_hyperplanes"))
            }
            (Some(_), None, Some(_)) => {
                return Err(IoError::Conflicting("generators", "mori_equations"))
            }
            (Some(g), None, None) => MoriCone::Generators(as_matrix(g, "generators")?),
            (None, Some(h), e) => MoriCone::Inequalities {
                hyperplanes: as_matrix(h, "mori_hyperplanes")?,
                equations: match e {
                    None => DMatrix::zeros(h11, 0),
                    Some(e) => as_matrix_or_empty(e, "mori_equations", h11)?,
                },
            },
            (None, None, Some(_)) => return Err(IoError::MissingField("mori_hyperplanes")),
            (None, None, None) if matches!(ambient, Ambient::Fan(_)) => MoriCone::FromFan,
            (None, None, None) => {
                return Err(IoError::MissingOneOf(&[
                    "generators",
                    "mori_hyperplanes",
                    "fan",
                ]))
            }
        };
        let grading_vector = get_field(hash, "grading_vector")
            .map(|y| as_int_vec(y, "grading_vector"))
            .transpose()?;
        let max_deg = get_field(hash, "max_deg")
            .map(|y| as_u32(y, "max_deg"))
            .transpose()?;
        let target_points = get_field(hash, "target_points")
            .map(|y| as_target_points(y, h11))
            .transpose()?;
        let min_points = get_field(hash, "min_points")
            .map(|y| as_u32(y, "min_points"))
            .transpose()?;
        let curve_classes = get_field(hash, "curve_classes")
            .map(|y| as_matrix(y, "curve_classes"))
            .transpose()?;
        let selection = match (max_deg, target_points, min_points, curve_classes) {
            (max_deg, Some(points), None, None) => Selection::TargetPoints { points, max_deg },
            (Some(d), None, None, None) => Selection::MaxDeg(d),
            (None, None, Some(n), None) => Selection::MinPoints(n),
            (None, None, None, Some(classes)) => Selection::CurveClasses(classes),
            (None, None, None, None) => {
                return Err(IoError::MissingOneOf(&[
                    "max_deg",
                    "target_points",
                    "min_points",
                    "curve_classes",
                ]))
            }
            (max_deg, target_points, min_points, curve_classes) => {
                let given: Vec<&'static str> = [
                    ("max_deg", max_deg.is_some()),
                    ("target_points", target_points.is_some()),
                    ("min_points", min_points.is_some()),
                    ("curve_classes", curve_classes.is_some()),
                ]
                .into_iter()
                .filter_map(|(f, given)| given.then_some(f))
                .collect();
                return Err(IoError::Conflicting(given[0], given[1]));
            }
        };

        let kind = match get_field(hash, "invariants") {
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
        let n_threads = get_field(hash, "n_threads")
            .map(|y| as_u32(y, "n_threads").map(|n| n as usize))
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
        let hilbert_bases = get_field(hash, "hilbert_bases")
            .map(|y| as_hilbert_bases(y, h11))
            .transpose()?;

        Ok(Input {
            name,
            problem: Problem {
                geometry: Geometry {
                    q,
                    intnums,
                    nefpart,
                    ambient,
                },
                cone,
                selection,
                grading_vector,
            },
            kind,
            options: Options {
                method,
                n_threads,
                certify,
                hilbert_bases,
                ..Options::default()
            },
        })
    }

    /// Sets the maximum degree, as the command line option does: of a [`Selection::MaxDeg`],
    /// or alongside [`Selection::TargetPoints`]. The other selections do not take one.
    pub fn set_max_deg(&mut self, max_deg: u32) -> Result<(), IoError> {
        match &mut self.problem.selection {
            Selection::MaxDeg(d) => *d = max_deg,
            Selection::TargetPoints { max_deg: d, .. } => *d = Some(max_deg),
            s => return Err(IoError::Conflicting("max_deg", s.name())),
        }
        Ok(())
    }

    /// Run the computation specified by this input. The results are sorted by degree, and
    /// then by curve class and reference surface.
    pub fn compute(&self) -> Result<Computed, Error> {
        compute(&self.problem, self.kind, &self.options)
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
        results: &Computed,
    ) -> std::io::Result<()> {
        match results {
            Computed::Gv(gv) => self.write(writer, gv, InvariantKind::GV),
            Computed::Gw(gw) => self.write(writer, gw, InvariantKind::GW),
        }
    }

    fn write<T: std::fmt::Display>(
        &self,
        writer: &mut impl Write,
        invariants: &Invariants<T>,
        kind: InvariantKind,
    ) -> std::io::Result<()> {
        let key = kind.as_str();
        writeln!(writer, "---")?;
        if let Some(name) = &self.name {
            writeln!(writer, "name: {}", quote_string(name))?;
        }
        writeln!(writer, "invariants: {key}")?;
        writeln!(writer, "dimension: {}", invariants.dimension)?;
        write!(writer, "grading_vector: ")?;
        write_int_seq(writer, invariants.grading_vector.iter())?;
        writeln!(writer)?;
        let report = &invariants.report;
        writeln!(
            writer,
            "report: {{method: {}, certified: {}, grading_vector_derived: {}}}",
            report.method.as_str(),
            report.certified,
            report.grading_vector_derived
        )?;
        if invariants.entries.is_empty() {
            return writeln!(writer, "results: []");
        }
        writeln!(writer, "results:")?;
        for e in &invariants.entries {
            write!(writer, "  - {{curve_class: ")?;
            write_int_seq(writer, e.curve.iter())?;
            write!(writer, ", degree: {}", e.degree)?;
            if let Some(surface) = e.surface {
                write!(writer, ", surface: {surface}")?;
            }
            match kind {
                InvariantKind::GV => writeln!(writer, ", {key}: {}}}", e.value)?,
                InvariantKind::GW => writeln!(writer, ", {key}: '{}'}}", e.value)?,
            }
        }
        Ok(())
    }
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

/// Read a value that must be a list of lists of non-negative indices, which are allowed to
/// have different lengths.
fn as_index_lists(data: &Yaml, field: &'static str) -> Result<Vec<Vec<usize>>, IoError> {
    let err = || IoError::invalid_field(field, "expected a list of lists of indices of divisors");
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

/// Read a list of vectors like [`as_matrix`], but allowing an empty list, which gives a
/// matrix with `n_rows` rows and no columns.
fn as_matrix_or_empty(data: &Yaml, field: &str, n_rows: usize) -> Result<DMatrix<i32>, IoError> {
    match data.as_vec() {
        Some(v) if v.is_empty() => Ok(DMatrix::zeros(n_rows, 0)),
        _ => as_matrix(data, field),
    }
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

/// Read the Hilbert bases of the cones $K_T$, a list of mappings with the `divisors` and the
/// `hilbert_basis`, which may be empty.
fn as_hilbert_bases(data: &Yaml, h11: usize) -> Result<Vec<CandidateCone>, IoError> {
    let field = "hilbert_bases";
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
            let hilbert_basis = as_matrix_or_empty(basis, field, h11)?;
            Ok(CandidateCone {
                divisors,
                hilbert_basis,
            })
        })
        .collect()
}

/// Read a fraction, given as an integer or as a string of the form `"p/q"` or `"p"`.
fn as_fraction(data: &Yaml) -> Result<Rational, IoError> {
    let invalid = || {
        IoError::invalid_field(
            "strata",
            "the entries of a class must be integers or fractions such as \"-1/2\"",
        )
    };
    match data {
        Yaml::Integer(n) => Ok(Rational::from(*n)),
        Yaml::String(s) => {
            let (num, den) = s.split_once('/').unwrap_or((s, "1"));
            let num = i64::from_str(num.trim()).map_err(|_| invalid())?;
            let den = i64::from_str(den.trim()).map_err(|_| invalid())?;
            if den == 0 {
                return Err(invalid());
            }
            Ok(Rational::from(num) / Rational::from(den))
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
            format!("expected vectors with {h11} entries, one per Kähler parameter"),
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
    use crate::{Invariant, Report};
    use malachite::Integer;
    use nalgebra::dmatrix;

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

    /// The field of an invalid-field error.
    fn field_of(data: &str) -> String {
        match load_one(data) {
            Err(IoError::InvalidField { field, .. }) => field,
            res => panic!("{res:?}"),
        }
    }

    #[test]
    fn test_parse_threefold() {
        let input = load_one(THREEFOLD).unwrap();
        assert_eq!(input.name.as_deref(), Some("test"));
        let problem = &input.problem;
        assert_eq!(problem.cone, MoriCone::Generators(dmatrix![0, 1; -1, 2]));
        assert_eq!(problem.grading_vector, Some(vec![3, -1]));
        let geometry = &problem.geometry;
        assert_eq!(geometry.q.shape(), (6, 2));
        assert_eq!(geometry.q[(2, 1)], -1);
        assert_eq!(geometry.intnums[&(0, 1, 1)], -1);
        assert!(geometry.nefpart.is_empty());
        assert_eq!(geometry.ambient, Ambient::Plain);
        assert_eq!(geometry.dimension(), 3);
        assert_eq!(problem.selection, Selection::MinPoints(20));
        assert_eq!(input.kind, InvariantKind::GV);
        assert_eq!(input.options, Options::default());
    }

    #[test]
    fn test_parse_selection() {
        let base = THREEFOLD.replace("min_points: 20\n", "");
        let selection =
            |extra: &str| load_one(&format!("{base}{extra}")).map(|i| i.problem.selection);
        assert_eq!(selection("max_deg: 5\n"), Ok(Selection::MaxDeg(5)));
        assert_eq!(
            selection("curve_classes: [[0, -1], [1, 0]]\n"),
            Ok(Selection::CurveClasses(dmatrix![0, 1; -1, 0]))
        );
        assert_eq!(
            selection(""),
            Err(IoError::MissingOneOf(&[
                "max_deg",
                "target_points",
                "min_points",
                "curve_classes"
            ]))
        );
        assert_eq!(
            selection("max_deg: 5\nmin_points: 3\n"),
            Err(IoError::Conflicting("max_deg", "min_points"))
        );
        assert_eq!(
            selection("max_deg: 5\ncurve_classes: [[1, 0]]\n"),
            Err(IoError::Conflicting("max_deg", "curve_classes"))
        );
        // The modular method does not support min_points or curve_classes.
        assert!(matches!(
            selection("min_points: 3\nmethod: modular\n"),
            Err(IoError::Unsupported(reason)) if reason.contains("min_points")
        ));
        assert_eq!(
            field_of(&format!("{base}curve_classes: [[1, 0, 0]]\n")),
            "curve_classes"
        );
    }

    #[test]
    fn test_set_max_deg() {
        let mut input = load_one(&THREEFOLD.replace("min_points: 20", "max_deg: 5")).unwrap();
        input.set_max_deg(7).unwrap();
        assert_eq!(input.problem.selection, Selection::MaxDeg(7));
        let mut input =
            load_one(&THREEFOLD.replace("min_points: 20", "target_points: [1, 2]")).unwrap();
        input.set_max_deg(2).unwrap();
        assert_eq!(
            input.problem.selection,
            Selection::TargetPoints {
                points: dmatrix![1; 2],
                max_deg: Some(2)
            }
        );
        let mut input = load_one(THREEFOLD).unwrap();
        assert_eq!(
            input.set_max_deg(2),
            Err(IoError::Conflicting("max_deg", "min_points"))
        );
    }

    #[test]
    fn test_parse_strata() {
        let base = THREEFOLD.replace("min_points: 20\n", "max_deg: 5\n");
        let input = load_one(&format!(
            "{base}strata:\n  - {{divisors: [3, 1, 2], class: [1, \"-1/2\"]}}\n  \
             - {{divisors: [0, 1, 5], class: [\"3\", 0]}}\n"
        ))
        .unwrap();
        let r = |n: i64, d: i64| Rational::from(n) / Rational::from(d);
        assert_eq!(
            input.problem.geometry.ambient,
            Ambient::Strata(vec![
                Stratum {
                    divisors: [3, 1, 2],
                    class: vec![r(1, 1), r(-1, 2)]
                },
                Stratum {
                    divisors: [0, 1, 5],
                    class: vec![r(3, 1), r(0, 1)]
                }
            ])
        );
        assert_eq!(
            load_one(&base).unwrap().problem.geometry.ambient,
            Ambient::Plain
        );

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
            assert_eq!(
                field_of(&format!("{base}strata: {bad}\n")),
                "strata",
                "{bad}"
            );
        }
        assert_eq!(
            load_one(&format!("{base}strata: [{{divisors: [1, 2, 3]}}]\n")),
            Err(IoError::MissingField("class"))
        );
        // Strata are only supported for threefold hypersurfaces.
        let fourfold = base.replace(
            "q: [[1, 1, 1, 0, 1, 2], [0, 0, -1, 1, 1, -1]]",
            "q: [[1, 1, 1, 0, 1, 2, 1], [0, 0, -1, 1, 1, -1, 0]]",
        );
        let strata = "strata: [{divisors: [1, 2, 3], class: [1, 1]}]\n";
        assert_eq!(field_of(&format!("{fourfold}{strata}")), "strata");
        assert!(load_one(&format!("{base}{strata}method: modular\n")).is_ok());
        // Nor certified, by the modular method.
        assert!(matches!(
            load_one(&format!("{base}{strata}certify: true\n")),
            Err(IoError::Unsupported(_))
        ));
        assert!(load_one(&format!("{base}{strata}certify: true\nmethod: reference\n")).is_ok());
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
        // The Mori cone and the grading vector are derived from the fan.
        let input = load_one(FAN).unwrap();
        assert_eq!(input.problem.cone, MoriCone::FromFan);
        assert_eq!(input.problem.grading_vector, None);
        assert!(matches!(&input.problem.geometry.ambient, Ambient::Fan(c) if c.len() == 8));
        let given = load_one(&format!("{FAN}grading_vector: [3, 1]\n")).unwrap();
        assert_eq!(given.problem.grading_vector, Some(vec![3, 1]));

        let strata = "strata: [{divisors: [1, 2, 4], class: [1, -1]}]\n";
        assert_eq!(
            load_one(&format!("{FAN}{strata}")),
            Err(IoError::Conflicting("fan", "strata"))
        );
        assert_eq!(
            field_of(&format!("{FAN}nefpart: [[0, 1, 2], [3, 4, 5]]\n")),
            "q"
        );
        assert_eq!(
            field_of(&FAN.replace("[1, 1, 1, 2]", "[1, 1, 1, 3]")),
            "intnums"
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
        assert_eq!(
            field_of(&format!("{FAN}generators: [[1, 0, 0]]\n")),
            "generators"
        );
        assert_eq!(
            field_of(&format!("{FAN}grading_vector: [1, -5]\n")),
            "grading_vector"
        );
        // Without a fan, the Mori cone is needed.
        let without = &FAN[..FAN.find("fan:").unwrap()];
        assert_eq!(
            load_one(without),
            Err(IoError::MissingOneOf(&[
                "generators",
                "mori_hyperplanes",
                "fan"
            ]))
        );
    }

    /// P(1,1,2,2,2)[8] with the inequalities of its Mori cone and the Hilbert bases.
    fn p11222_cones() -> String {
        let mut data = "\
grading_vector: [1, 1]
q: [[0, 0, 1, 1, 1, 1], [1, 1, 0, 0, 0, -2]]
intnums: [[0, 0, 0, 8], [0, 0, 1, 4]]
max_deg: 6
mori_hyperplanes: [[1, 0], [0, 1]]
hilbert_bases:
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
    fn test_parse_hilbert_bases() {
        let data = p11222_cones();
        let input = load_one(&data).unwrap();
        assert_eq!(
            input.problem.cone,
            MoriCone::Inequalities {
                hyperplanes: dmatrix![1, 0; 0, 1],
                equations: DMatrix::zeros(2, 0)
            }
        );
        let bases = input.options.hilbert_bases.unwrap();
        assert_eq!(bases.len(), 21);
        assert_eq!(bases[1].divisors, [0, 1]);
        assert_eq!(bases[1].hilbert_basis, dmatrix![1, 2; 0, 1]);
        assert_eq!(
            load_one(&format!("{data}generators: [[1, 0], [0, 1]]\n")),
            Err(IoError::Conflicting("generators", "mori_hyperplanes"))
        );

        // The face spanned by (1, 0), which meets no divisor negatively, so it is every K_T.
        let face = data
            .replace("max_deg: 6\n", "max_deg: 6\nmori_equations: [[0, 1]]\n")
            .replace("[[1, 0], [0, 1]]}", "[[1, 0]]}")
            .replace("[[1, 0], [2, 1]]}", "[[1, 0]]}");
        let input = load_one(&face).unwrap();
        assert!(matches!(
            &input.problem.cone,
            MoriCone::Inequalities { equations, .. } if *equations == dmatrix![0; 1]
        ));
        // An empty Hilbert basis is read too (a K_T can be just the origin).
        let empty = data.replacen("[[1, 0], [2, 1]]}", "[]}", 1);
        let bases = load_one(&empty).unwrap().options.hilbert_bases.unwrap();
        assert_eq!(bases[0].hilbert_basis.shape(), (2, 0));

        // The inequalities alone describe the Mori cone; the Hilbert bases are optional.
        let head = &data[..data.find("hilbert_bases:").unwrap()];
        assert!(load_one(head).unwrap().options.hilbert_bases.is_none());
        let no_hyperplanes = head.replace("mori_hyperplanes: [[1, 0], [0, 1]]\n", "");
        assert_eq!(
            load_one(&format!("{no_hyperplanes}mori_equations: [[0, 1]]\n")),
            Err(IoError::MissingField("mori_hyperplanes"))
        );
        // The reference method, which runs for min_points, has no use for them.
        assert!(load_one(&data.replace("max_deg: 6\n", "min_points: 5\n")).is_ok());
        assert!(load_one(&format!("{data}method: reference\n")).is_ok());

        // They are checked against the Mori cone.
        let last = data.rfind("  - ").unwrap();
        assert_eq!(field_of(&data[..last]), "hilbert_bases");
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
            assert_eq!(
                field_of(&data.replacen(from, to, 1)),
                "hilbert_bases",
                "{to}"
            );
        }
        assert_eq!(
            field_of(&data.replace(
                "mori_hyperplanes: [[1, 0], [0, 1]]",
                "mori_hyperplanes: [[1, 0, 0]]"
            )),
            "mori_hyperplanes"
        );
        assert_eq!(
            field_of(&format!("{data}mori_equations: [[1]]\n")),
            "mori_equations"
        );
    }

    #[test]
    fn test_removed_fields() {
        // The fields of earlier versions say what replaced them.
        for (field, extra) in [
            ("cones", "cones: []\n"),
            ("is_threefold", "is_threefold: true\n"),
        ] {
            assert_eq!(field_of(&format!("{THREEFOLD}{extra}")), field);
        }
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
            assert_eq!(input.options.method, method);
        }
        assert_eq!(field_of(&format!("{base}method: fast\n")), "method");
        assert!(load_one(&format!("{THREEFOLD}method: auto\n")).is_ok());
    }

    #[test]
    fn test_parse_options() {
        let base = THREEFOLD.replace("min_points: 20\n", "max_deg: 5\n");
        assert!(!load_one(&base).unwrap().options.certify);
        let input = load_one(&format!("{base}certify: true\nn_threads: 3\n")).unwrap();
        assert!(input.options.certify);
        assert_eq!(input.options.n_threads, Some(3));
        assert_eq!(field_of(&format!("{base}certify: 1\n")), "certify");
        assert_eq!(field_of(&format!("{base}n_threads: 0\n")), "n_threads");
    }

    #[test]
    fn test_parse_multiple_documents() {
        let data = format!("---\n{THREEFOLD}---\n{THREEFOLD}invariants: gw\n");
        let inputs = Input::load_all(&data).unwrap();
        assert_eq!(inputs.len(), 2);
        assert_eq!(inputs[0].kind, InvariantKind::GV);
        assert_eq!(inputs[1].kind, InvariantKind::GW);
    }

    #[test]
    fn test_parse_intnums_mapping() {
        let data = THREEFOLD.replace(
            "intnums: [[0, 0, 0, 2], [0, 0, 1, 1], [0, 1, 1, -1], [1, 1, 1, 5]]",
            "intnums: {\"0,0,0\": 2, \"0 0 1\": 1, \"0,1,1\": -1, \"1,1,1\": 5}",
        );
        let input = load_one(&data).unwrap();
        assert_eq!(input.problem, load_one(THREEFOLD).unwrap().problem);
    }

    #[test]
    fn test_parse_target_points() {
        let base = THREEFOLD.replace("min_points: 20\n", "");
        let selection =
            |extra: &str| load_one(&format!("{base}{extra}")).map(|i| i.problem.selection);
        assert_eq!(
            selection("target_points: [1, 2]\n"),
            Ok(Selection::TargetPoints {
                points: dmatrix![1; 2],
                max_deg: None
            })
        );
        assert_eq!(
            selection("target_points: [[1, 2], [3, 4]]\nmax_deg: 3\nmethod: modular\n"),
            Ok(Selection::TargetPoints {
                points: dmatrix![1, 3; 2, 4],
                max_deg: Some(3)
            })
        );
        assert_eq!(
            selection("target_points: [1, 2]\nmin_points: 20\n"),
            Err(IoError::Conflicting("target_points", "min_points"))
        );
        for bad in ["[[1, 2, 3]]", "[0, 0]", "[[1, 2], [-1, -1]]"] {
            assert_eq!(
                field_of(&format!("{base}target_points: {bad}\n")),
                "target_points"
            );
        }
    }

    #[test]
    fn test_dimension() {
        // Eight rays, two Kähler parameters, and a codimension-two CY, i.e. a fourfold.
        let data = THREEFOLD.replace(
            "q: [[1, 1, 1, 0, 1, 2], [0, 0, -1, 1, 1, -1]]",
            "q: [[1, 1, 1, 0, 1, 2, 1, 0], [0, 0, -1, 1, 1, -1, 0, 1]]",
        );
        let data = format!("{data}nefpart: [[0, 1], [2, 3, 4, 5, 6, 7]]\n");
        let input = load_one(&data).unwrap();
        assert_eq!(input.problem.geometry.nefpart.len(), 2);
        assert_eq!(input.problem.geometry.dimension(), 4);
        // The nef partition must split the rays: indices of the rays, each in one part.
        for parts in [
            "[2, 3, 4, 5, 6, 8]",
            "[1, 2, 3, 4, 5, 6, 7]",
            "[2, 3, 4, 5, 6]",
            "[2, 3, 4, 5, 6, -7]",
        ] {
            assert_eq!(
                field_of(&data.replace("[2, 3, 4, 5, 6, 7]", parts)),
                "nefpart",
                "{parts}"
            );
        }
        // The resulting CY must be at least three-dimensional.
        assert_eq!(
            field_of(&format!("{THREEFOLD}nefpart: [[0, 1], [2, 3, 4, 5]]\n")),
            "q"
        );
        // A transposed q says so.
        let transposed = THREEFOLD.replace(
            "q: [[1, 1, 1, 0, 1, 2], [0, 0, -1, 1, 1, -1]]",
            "q: [[1, 0], [1, 0], [1, -1], [0, 1], [1, 1], [2, -1]]",
        );
        assert!(matches!(
            load_one(&transposed),
            Err(IoError::InvalidField { field, reason }) if field == "q" && reason.contains("transposed")
        ));
    }

    #[test]
    fn test_parse_errors() {
        assert_eq!(Input::load_all(""), Err(IoError::EmptyInput));
        assert_eq!(Input::load_all("- 1"), Err(IoError::NotAMapping));
        assert_eq!(
            load_one("generators: [[1]]"),
            Err(IoError::MissingField("q"))
        );
        assert_eq!(
            load_one(&format!("{THREEFOLD}bad_field: 1\n")),
            Err(IoError::UnknownField("bad_field".to_owned()))
        );
        assert_eq!(
            field_of(&THREEFOLD.replace("[[0, -1], [1, 2]]", "[[1, 2], [3]]")),
            "generators"
        );
        assert_eq!(
            field_of(&THREEFOLD.replace("[[0, -1], [1, 2]]", "[[1, 2, 0]]")),
            "generators"
        );
        assert_eq!(
            field_of(&THREEFOLD.replace("[3, -1]", "[3, -1, 0]")),
            "grading_vector"
        );
        assert_eq!(
            field_of(&THREEFOLD.replace("[0, 1, 1, -1]", "[0, 1, 1, -1], [1, 0, 1, 2]")),
            "intnums"
        );
        // So must indices that are not those of the divisors.
        assert_eq!(
            field_of(&THREEFOLD.replace("[1, 1, 1, 5]", "[1, 1, 2, 5]")),
            "intnums"
        );
        // Intersection numbers that are all zero must be rejected, not panic.
        assert_eq!(
            field_of(&THREEFOLD.replace(
                "intnums: [[0, 0, 0, 2], [0, 0, 1, 1], [0, 1, 1, -1], [1, 1, 1, 5]]",
                "intnums: [[0, 0, 0, 0], [0, 0, 1, 0]]",
            )),
            "intnums"
        );
        assert!(Input::load_all("generators: [[1]], bad").is_err());
    }

    /// GV invariants with the given entries.
    fn gv(entries: Vec<Invariant<Integer>>, dimension: usize) -> Computed {
        Computed::Gv(Invariants {
            dimension,
            grading_vector: vec![3, -1],
            entries,
            report: Report {
                method: Method::Modular,
                certified: false,
                grading_vector_derived: false,
                elapsed: std::time::Duration::ZERO,
            },
        })
    }

    #[test]
    fn test_write_results() {
        let mut input = load_one(THREEFOLD).unwrap();
        let entry = |curve: Vec<i32>, surface, degree, value: i32| Invariant {
            curve,
            surface,
            degree,
            value: Integer::from(value),
        };
        let results = gv(
            vec![
                entry(vec![0, -1], None, 1, -4),
                entry(vec![1, 0], None, 3, 7524),
            ],
            3,
        );
        let mut buf = Vec::new();
        input.write_results(&mut buf, &results).unwrap();
        assert_eq!(
            String::from_utf8(buf).unwrap(),
            "---\nname: test\ninvariants: gv\ndimension: 3\ngrading_vector: [3, -1]\n\
             report: {method: modular, certified: false, grading_vector_derived: false}\n\
             results:\n  \
             - {curve_class: [0, -1], degree: 1, gv: -4}\n  \
             - {curve_class: [1, 0], degree: 3, gv: 7524}\n"
        );

        input.name = None;
        let results = Computed::Gw(Invariants {
            dimension: 4,
            grading_vector: vec![3, -1],
            entries: vec![Invariant {
                curve: vec![0, -1],
                surface: Some(2),
                degree: 1,
                value: Rational::from(-4),
            }],
            report: Report {
                method: Method::Reference,
                certified: true,
                grading_vector_derived: true,
                elapsed: std::time::Duration::ZERO,
            },
        });
        let mut buf = Vec::new();
        input.write_results(&mut buf, &results).unwrap();
        assert_eq!(
            String::from_utf8(buf).unwrap(),
            "---\ninvariants: gw\ndimension: 4\ngrading_vector: [3, -1]\n\
             report: {method: reference, certified: true, grading_vector_derived: true}\n\
             results:\n  \
             - {curve_class: [0, -1], degree: 1, surface: 2, gw: '-4'}\n"
        );

        let mut buf = Vec::new();
        input.write_results(&mut buf, &gv(vec![], 3)).unwrap();
        assert!(String::from_utf8(buf).unwrap().ends_with("results: []\n"));
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
    fn test_output_is_valid_yaml() {
        let input = load_one(THREEFOLD).unwrap();
        let results = input.compute().unwrap();
        let Computed::Gv(gv) = &results else {
            panic!("GW invariants instead of GV")
        };
        assert!(gv.entries.len() > 1);
        assert_eq!(gv.report.method, Method::Reference);
        // Sorted by degree, and then by curve class.
        assert!(gv
            .entries
            .windows(2)
            .all(|w| (w[0].degree, &w[0].curve) < (w[1].degree, &w[1].curve)));

        let mut buf = Vec::new();
        input.write_results(&mut buf, &results).unwrap();
        let output = YamlLoader::load_from_str(&String::from_utf8(buf).unwrap()).unwrap();
        assert_eq!(output.len(), 1);
        assert_eq!(output[0]["name"].as_str(), Some("test"));
        assert_eq!(output[0]["invariants"].as_str(), Some("gv"));
        assert_eq!(output[0]["dimension"].as_i64(), Some(3));
        assert_eq!(output[0]["report"]["method"].as_str(), Some("reference"));
        assert_eq!(output[0]["report"]["certified"].as_bool(), Some(true));
        assert_eq!(
            as_int_vec(&output[0]["grading_vector"], "").unwrap(),
            [3, -1]
        );
        let parsed = output[0]["results"].as_vec().unwrap();
        assert_eq!(parsed.len(), gv.entries.len());
        for (entry, e) in parsed.iter().zip(&gv.entries) {
            assert_eq!(as_int_vec(&entry["curve_class"], "").unwrap(), e.curve);
            assert_eq!(entry["degree"].as_i64(), Some(e.degree as i64));
            assert_eq!(
                entry["gv"].as_i64().unwrap().to_string(),
                e.value.to_string()
            );
        }
    }
}
