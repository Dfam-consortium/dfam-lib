//! Error type shared across the crate.

use thiserror::Error;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Error)]
pub enum Error {
    /// A substitution matrix could not be parsed, or is internally inconsistent.
    #[error("matrix error: {0}")]
    Matrix(String),

    /// An alignment's edit script does not agree with its coordinates, or the
    /// supplied sequences are too short for the span it claims.
    #[error("alignment error: {0}")]
    Alignment(String),

    /// A multiple alignment's rows disagree in width, or a row is malformed.
    #[error("MSA error: {0}")]
    Msa(String),

    /// Rescoring or divergence calculation could not proceed.
    #[error("scoring error: {0}")]
    Scoring(String),

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}
