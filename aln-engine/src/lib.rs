//! Aligner and search-engine abstractions, plus parallel drivers.
//!
//! # Two traits, not one
//!
//! The suite uses two genuinely different kinds of aligner, and collapsing them
//! into a single trait produces a parameter set where most fields are
//! meaningless for any given implementor:
//!
//! * [`PairwiseAligner`] — two sequences in, one alignment out.  Full dynamic
//!   programming, no seeding, no database.  This is parasail, Farrar's striped
//!   SSE2 code, Monardo, and Phil Green's SWAT.  It is the direct analogue of
//!   GIRI's `SWAligner`, and it is the *only* one `acons`/`autocons` need.
//! * [`SearchEngine`] — a query and a subject database in, many HSPs out, with
//!   seeding, masking, complexity adjustment and score cutoffs.  This is
//!   rmblast, crossmatch and HMMER, and it is the analogue of RepeatMasker's
//!   `SearchEngineI`.
//!
//! # No global state
//!
//! GIRI configures alignment through static members — `MultipleAlignment::
//! setAligner`, `PairwiseAlignment::setScoreMatrix`, `ThreadedAligner::
//! setMultithreaded`.  Everything here is owned by the aligner instance or
//! passed as a parameter, so two differently-configured aligners can run
//! concurrently in one process.

pub mod driver;
pub mod engine;
pub mod params;
pub mod traits;

pub use driver::{align_pairs, all_vs_all, one_to_many};
pub use engine::{SearchEngine, SearchParams, SeqSource};
pub use params::{AlignMode, AlignParams, AlignerCaps};
pub use traits::{DynAligner, PairwiseAligner};

use thiserror::Error;

pub type Result<T> = std::result::Result<T, EngineError>;

#[derive(Debug, Error)]
pub enum EngineError {
    #[error(transparent)]
    Core(#[from] aln_core::Error),

    /// The aligner cannot honour the requested parameters — for example a
    /// score-only backend asked to produce a traceback.
    #[error("unsupported by {aligner}: {what}")]
    Unsupported { aligner: &'static str, what: String },

    /// The backend failed at run time (allocation, FFI error, overflow).
    #[error("{aligner} failed: {message}")]
    Backend { aligner: &'static str, message: String },

    #[error("invalid parameters: {0}")]
    Params(String),

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}

impl EngineError {
    pub fn unsupported(aligner: &'static str, what: impl Into<String>) -> Self {
        EngineError::Unsupported { aligner, what: what.into() }
    }

    pub fn backend(aligner: &'static str, message: impl Into<String>) -> Self {
        EngineError::Backend { aligner, message: message.into() }
    }
}
