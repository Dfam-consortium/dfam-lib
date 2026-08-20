//! Core types for the Dfam / GIRI alignment stack.
//!
//! This crate holds everything that is *not* an aligner: sequences, scoring
//! matrices, the pairwise alignment representation, divergence and rescoring
//! statistics, and the multiple alignment.  Aligner backends
//! (`aln-parasail`, `aln-farrar`, `aln-rmblast`, `aln-reference`) depend on
//! this crate and on `aln-engine`; nothing here depends on them.
//!
//! # Conventions, in one place
//!
//! These are the three things to internalise before using the crate; each is
//! documented at length in its own module.
//!
//! 1. **Coordinates are 0-based, half-open, forward-strand** ([`align`]).
//!    1-based closed coordinates exist only at I/O boundaries.
//! 2. **Gap `-` and padding `' '` are different characters** ([`seq`]).  A gap
//!    means "present but deleted"; padding means "not present".
//! 3. **Matrices are asymmetric and indexed `matrix[subject][query]`**
//!    ([`matrix`]), where subject is the consensus and query is genomic.
//!
//! # Layout
//!
//! | module | contents |
//! |--------|----------|
//! | [`seq`] | strand, IUPAC tables, gap/padding conventions, [`Sequence`] |
//! | [`matrix`] | [`SubstMatrix`]: parsing, frequencies, Karlin-Altschul lambda |
//! | [`align`] | [`Alignment`], [`EditScript`] — one pairwise result |
//! | [`stats`] | Kimura / K2P-Gap divergence, rescoring, complexity adjustment |
//! | [`msa`] | [`MultiAlign`], and assembly of one from pairwise alignments |
//! | [`result`] | [`SearchResult`] — an alignment plus reporting annotation |
//! | [`fmt`] | RepeatMasker-compatible `.out` and crossmatch writers |

pub mod align;
pub mod consensus;
pub mod error;
pub mod fmt;
pub mod giri;
pub mod io;
pub mod matrix;
pub mod msa;
pub mod result;
pub mod seq;
pub mod stats;
/// UCSC `.2bit` random-access reader.
///
/// Sequence storage rather than alignment, but every tool that extends a
/// consensus into genomic flanks needs it, and a second copy in a second
/// project is how two readers drift apart.
pub mod twobit;

/// Reading crossmatch-style pairwise output.
///
/// The counterpart to [`fmt::to_crossmatch`], which has been able to write this
/// format all along while nothing here could read it back. Parsing lived in
/// dfam-curator, so any other consumer of the format had to depend on a
/// curation tool to get at it.
pub mod crossmatch;

pub use align::{Alignment, EditOp, EditScript, IdentityCounts};
pub use consensus::{build_consensus_from_sequences, ConsensusParams};
pub use error::{Error, Result};
pub use fmt::{
    annotation_line, to_caf, to_cigar_record, to_crossmatch, to_out_line, AlignmentMode,
};
pub use matrix::SubstMatrix;
pub use msa::{assemble_msa, InsertionPolicy, MsaMember, MultiAlign, SequenceRow};
pub use result::SearchResult;
pub use seq::{Sequence, Strand};
pub use stats::{
    classify_subst, k2p, k2p_gap_divergence, kimura_divergence, kimura_stats, mean_kimura,
    rescore, Divergence, KimuraStats,
    RescoreParams, RescoreResult, SubstClass,
};
