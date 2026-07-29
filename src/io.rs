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
//! # How many curve classes to compute. At most one of max_deg, min_points, and
//! # target_points may be given. When none of them is given, the generators are
//! # taken to be the complete list of curve classes to use.
//! min_points: 100
//! ```
//!
//! The results are written back as one YAML document per input document.
//!
//! ```yaml
//! ---
//! name: quintic-like threefold
//! invariants: gv
//! is_threefold: true
//! results:
//!   - {curve_class: [1, 0], gv: 3}
//!   - {curve_class: [0, 1], gv: -6}
//! ```
//!
//! See [`Input`] for the full list of fields.

pub mod error;

use crate::hkty::compute_gvgw_strings;
use crate::misc::process_int_nums;
use error::IoError;
use nalgebra::{DMatrix, DVector, RowDVector};
use std::collections::HashMap;
use std::io::Write;
use yaml_rust2::yaml::{Hash, Yaml};
use yaml_rust2::YamlLoader;

/// The default number of coefficients that are kept in the number pools.
pub const DEFAULT_POOL_SIZE: usize = 1000;

/// The fields that are accepted in an input document.
const INPUT_FIELDS: [&str; 14] = [
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
    "prec",
    "n_threads",
    "pool_size",
];

/// The specification of a single computation of GV or GW invariants.
///
/// The matrices are stored using the same convention as the rest of the crate, i.e.
/// the vectors of the input data are the *columns* of the matrices.
#[derive(Debug, Clone, PartialEq)]
pub struct Input {
    /// An optional label that is echoed back in the output.
    pub name: Option<String>,
    /// The generators of the Mori cone.
    pub generators: DMatrix<i32>,
    /// A vector that has a positive inner product with every generator.
    pub grading_vector: RowDVector<i32>,
    /// The GLSM charge matrix.
    pub q: DMatrix<i32>,
    /// The nef partition. It is empty for hypersurfaces.
    pub nefpart: Vec<DVector<i32>>,
    /// The triple intersection numbers.
    pub intnums: HashMap<(usize, usize, usize), i32>,
    /// Whether to compute GV invariants instead of GW invariants.
    pub find_gv: bool,
    /// Whether the CY is a threefold.
    pub is_threefold: bool,
    /// The maximum degree of the curve classes that are computed.
    pub max_deg: Option<u32>,
    /// The minimum number of curve classes that are computed.
    pub min_points: Option<u32>,
    /// The curve classes that must be included in the computation.
    pub target_points: Option<DMatrix<i32>>,
    /// The precision, in bits, of the floating-point arithmetic. Exact rational
    /// arithmetic is used when it is [`None`].
    pub prec: Option<u32>,
    /// The number of threads to use. It is deduced from the machine when [`None`].
    pub n_threads: Option<u32>,
    /// The number of coefficients that are kept in the number pools.
    pub pool_size: usize,
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
        let generators = as_matrix(require_field(hash, "generators")?, "generators")?;
        let h11 = generators.nrows();
        let grading_vector = as_int_vec(require_field(hash, "grading_vector")?, "grading_vector")?;
        if grading_vector.len() != h11 {
            return Err(IoError::invalid_field(
                "grading_vector",
                format!("expected {h11} entries to match the generators"),
            ));
        }
        let grading_vector = RowDVector::from_row_slice(&grading_vector);
        let q = as_matrix(require_field(hash, "q")?, "q")?;
        if q.ncols() != h11 {
            return Err(IoError::invalid_field(
                "q",
                format!("expected {h11} charge vectors to match the generators"),
            ));
        }
        let intnums = as_intnums(require_field(hash, "intnums")?)?;
        let nefpart = match get_field(hash, "nefpart") {
            None => Vec::new(),
            Some(y) => as_vectors(y, "nefpart")?,
        };

        let find_gv = match get_field(hash, "invariants") {
            None => true,
            Some(y) => match y.as_str().map(str::to_lowercase).as_deref() {
                Some("gv") => true,
                Some("gw") => false,
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
        let prec = get_field(hash, "prec")
            .map(|y| as_u32(y, "prec"))
            .transpose()?;
        if prec == Some(0) {
            return Err(IoError::invalid_field(
                "prec",
                "expected a positive integer",
            ));
        }
        let n_threads = get_field(hash, "n_threads")
            .map(|y| as_u32(y, "n_threads"))
            .transpose()?;
        let pool_size = get_field(hash, "pool_size")
            .map(|y| as_usize(y, "pool_size"))
            .transpose()?
            .unwrap_or(DEFAULT_POOL_SIZE);

        let is_threefold = match get_field(hash, "is_threefold") {
            None => infer_is_threefold(&q, &nefpart),
            Some(y) => y
                .as_bool()
                .ok_or_else(|| IoError::invalid_field("is_threefold", "expected a boolean"))?,
        };
        let n_truncations = [
            max_deg.is_some(),
            min_points.is_some(),
            target_points.is_some(),
        ]
        .iter()
        .filter(|c| **c)
        .count();
        if n_truncations > 1 {
            return Err(IoError::ConflictingTruncations);
        }
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
        // The intersection numbers are canonicalized differently depending on the
        // dimension of the CY, so they can only be checked once it is known.
        process_int_nums(intnums.clone(), is_threefold)
            .map_err(|e| IoError::invalid_field("intnums", e.to_string()))?;

        Ok(Input {
            name,
            generators,
            grading_vector,
            q,
            nefpart,
            intnums,
            find_gv,
            is_threefold,
            max_deg,
            min_points,
            target_points,
            prec,
            n_threads,
            pool_size,
        })
    }

    /// Run the computation specified by this input.
    ///
    /// The results are sorted by degree, and then by curve class and reference surface
    /// index. The last two are only used to break ties, since the order in which
    /// curve classes of equal degree are computed is not deterministic.
    pub fn compute(&self) -> Vec<InvariantResult> {
        let mut results = self.compute_unsorted();
        results.sort_unstable_by(|a, b| {
            let (curve_a, curve_b) = (&a.0 .0, &b.0 .0);
            self.degree(curve_a)
                .cmp(&self.degree(curve_b))
                .then_with(|| curve_a.as_slice().cmp(curve_b.as_slice()))
                .then_with(|| a.0 .1.cmp(&b.0 .1))
        });
        results
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
    fn compute_unsorted(&self) -> Vec<InvariantResult> {
        compute_gvgw_strings(
            self.generators.clone(),
            self.grading_vector.clone(),
            self.q.clone(),
            self.nefpart.clone(),
            self.intnums.clone(),
            self.find_gv,
            self.is_threefold,
            self.max_deg,
            self.min_points,
            self.target_points.clone(),
            self.n_threads,
            self.pool_size,
            self.prec,
        )
    }

    /// Write the results of a computation as a YAML document.
    ///
    /// GV invariants are written as integers, while GW invariants are written as
    /// quoted strings, since they are either fractions or floating-point numbers with
    /// a precision that YAML cannot represent.
    pub fn write_results(
        &self,
        writer: &mut impl Write,
        results: &[InvariantResult],
    ) -> std::io::Result<()> {
        let key = if self.find_gv { "gv" } else { "gw" };
        writeln!(writer, "---")?;
        if let Some(name) = &self.name {
            writeln!(writer, "name: {}", quote_string(name))?;
        }
        writeln!(writer, "invariants: {key}")?;
        writeln!(writer, "is_threefold: {}", self.is_threefold)?;
        if results.is_empty() {
            return writeln!(writer, "results: []");
        }
        writeln!(writer, "results:")?;
        for ((curve, idx), value) in results {
            write!(writer, "  - {{curve_class: [")?;
            for (i, c) in curve.iter().enumerate() {
                if i > 0 {
                    write!(writer, ", ")?;
                }
                write!(writer, "{c}")?;
            }
            write!(writer, "]")?;
            if !self.is_threefold {
                write!(writer, ", surface_index: {idx}")?;
            }
            if self.find_gv {
                writeln!(writer, ", {key}: {value}}}")?;
            } else {
                writeln!(writer, ", {key}: '{value}'}}")?;
            }
        }
        Ok(())
    }
}

/// Deduce whether the CY is a threefold from the shape of the GLSM charge matrix and
/// the nef partition.
fn infer_is_threefold(q: &DMatrix<i32>, nefpart: &[DVector<i32>]) -> bool {
    let ambient_dim = q.nrows() as i64 - q.ncols() as i64;
    let cy_codim = if nefpart.is_empty() {
        1
    } else {
        nefpart.len() as i64
    };
    ambient_dim - cy_codim == 3
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

/// Read a value that must be a non-negative integer.
fn as_usize(data: &Yaml, field: &str) -> Result<usize, IoError> {
    usize::try_from(as_i64(data, field)?)
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
            format!("expected vectors with {h11} entries to match the generators"),
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
        assert!(input.find_gv);
        assert!(input.is_threefold);
        assert_eq!(input.min_points, Some(20));
        assert_eq!(input.max_deg, None);
        assert_eq!(input.target_points, None);
        assert_eq!(input.prec, None);
        assert_eq!(input.pool_size, DEFAULT_POOL_SIZE);
    }

    #[test]
    fn test_parse_multiple_documents() {
        let data = format!("---\n{THREEFOLD}---\n{THREEFOLD}invariants: gw\n");
        let inputs = Input::load_all(&data).unwrap();
        assert_eq!(inputs.len(), 2);
        assert!(inputs[0].find_gv);
        assert!(!inputs[1].find_gv);
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

        assert_eq!(
            load_one(&format!("{THREEFOLD}target_points: [1, 2]\n")),
            Err(IoError::ConflictingTruncations)
        );
        assert!(matches!(
            load_one(&format!("{base}target_points: [[1, 2, 3]]\n")),
            Err(IoError::InvalidField { field, .. }) if field == "target_points"
        ));
    }

    #[test]
    fn test_infer_is_threefold() {
        assert!(load_one(THREEFOLD).unwrap().is_threefold);

        // Eight rays, two Kähler parameters, and a codimension-two CY, i.e. a fourfold.
        let data = THREEFOLD.replace(
            "q: [[1, 1, 1, 0, 1, 2], [0, 0, -1, 1, 1, -1]]",
            "q: [[1, 1, 1, 0, 1, 2, 1, 0], [0, 0, -1, 1, 1, -1, 0, 1]]",
        );
        let data = format!("{data}nefpart: [[0, 1], [2, 3, 4, 5, 6, 7]]\n");
        let input = load_one(&data).unwrap();
        assert_eq!(input.nefpart.len(), 2);
        assert!(!input.is_threefold);

        assert!(
            load_one(&format!("{data}is_threefold: true\n"))
                .unwrap()
                .is_threefold
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
            load_one(&format!("{THREEFOLD}prec: -1\n")),
            Err(IoError::InvalidField { field, .. }) if field == "prec"
        ));
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
        assert!(Input::load_all("generators: [[1]], bad").is_err());
    }

    #[test]
    fn test_write_results() {
        let mut input = load_one(THREEFOLD).unwrap();
        let results = vec![
            ((dvector![1, 0], 0), "3".to_owned()),
            ((dvector![0, 1], 2), "-6".to_owned()),
        ];

        let mut buf = Vec::new();
        input.write_results(&mut buf, &results).unwrap();
        assert_eq!(
            String::from_utf8(buf).unwrap(),
            "---\nname: test\ninvariants: gv\nis_threefold: true\nresults:\n  \
             - {curve_class: [1, 0], gv: 3}\n  - {curve_class: [0, 1], gv: -6}\n"
        );

        input.name = None;
        input.find_gv = false;
        input.is_threefold = false;
        let mut buf = Vec::new();
        input.write_results(&mut buf, &results).unwrap();
        assert_eq!(
            String::from_utf8(buf).unwrap(),
            "---\ninvariants: gw\nis_threefold: false\nresults:\n  \
             - {curve_class: [1, 0], surface_index: 0, gw: '3'}\n  \
             - {curve_class: [0, 1], surface_index: 2, gw: '-6'}\n"
        );

        let mut buf = Vec::new();
        input.write_results(&mut buf, &[]).unwrap();
        assert_eq!(
            String::from_utf8(buf).unwrap(),
            "---\ninvariants: gw\nis_threefold: false\nresults: []\n"
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
        let results = input.compute();
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
        let results = input.compute();
        assert!(!results.is_empty());

        let mut buf = Vec::new();
        input.write_results(&mut buf, &results).unwrap();
        let output = YamlLoader::load_from_str(&String::from_utf8(buf).unwrap()).unwrap();
        assert_eq!(output.len(), 1);
        assert_eq!(output[0]["name"].as_str(), Some("test"));
        assert_eq!(output[0]["invariants"].as_str(), Some("gv"));
        assert_eq!(output[0]["is_threefold"].as_bool(), Some(true));
        let parsed = output[0]["results"].as_vec().unwrap();
        assert_eq!(parsed.len(), results.len());
        for (entry, ((curve, _), value)) in parsed.iter().zip(results.iter()) {
            assert_eq!(
                as_int_vec(&entry["curve_class"], "").unwrap(),
                curve.as_slice()
            );
            assert_eq!(entry["gv"].as_i64().unwrap().to_string(), *value);
        }
    }
}
