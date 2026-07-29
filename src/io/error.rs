//! A module for errors that can happen when reading and writing data.

use core::fmt;

/// An error enum for input and output errors.
#[derive(Debug, Clone, PartialEq)]
pub enum IoError {
    /// The input could not be parsed as YAML.
    InvalidYaml(String),
    /// The input contains no YAML document.
    EmptyInput,
    /// An input document is not a mapping of field names to values.
    NotAMapping,
    /// A required field is missing.
    MissingField(&'static str),
    /// A field that is not part of the input schema was found.
    UnknownField(String),
    /// More than one of the fields that select the curve classes to compute was given.
    ConflictingTruncations,
    /// A field has an invalid type, shape, or value.
    InvalidField {
        /// The name of the offending field.
        field: String,
        /// The reason why the value is not acceptable.
        reason: String,
    },
}

impl IoError {
    /// Construct an [`IoError::InvalidField`] error.
    pub fn invalid_field(field: impl Into<String>, reason: impl Into<String>) -> Self {
        IoError::InvalidField {
            field: field.into(),
            reason: reason.into(),
        }
    }
}

/// Implement Display trait for IoError.
impl fmt::Display for IoError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            IoError::InvalidYaml(msg) => write!(f, "the input is not valid YAML: {msg}"),
            IoError::EmptyInput => write!(f, "the input contains no YAML document"),
            IoError::NotAMapping => write!(
                f,
                "each input document must be a mapping of field names to values"
            ),
            IoError::MissingField(field) => write!(f, "the required field \"{field}\" is missing"),
            IoError::UnknownField(field) => write!(f, "unknown field \"{field}\""),
            IoError::ConflictingTruncations => write!(
                f,
                "at most one of the fields \"max_deg\", \"min_points\", and \
                 \"target_points\" may be given"
            ),
            IoError::InvalidField { field, reason } => {
                write!(f, "invalid value for the field \"{field}\": {reason}")
            }
        }
    }
}

impl std::error::Error for IoError {}
