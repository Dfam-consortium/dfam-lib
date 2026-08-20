//! The seeded database-search abstraction — the analogue of RepeatMasker's
//! `SearchEngineI`.
//!
//! Where [`PairwiseAligner`](crate::PairwiseAligner) does full dynamic
//! programming on two sequences, a [`SearchEngine`] seeds, extends, filters and
//! reports many HSPs across a database.  rmblast, crossmatch and HMMER go here.
//!
//! `SearchEngineI` exposes its configuration as ~17 getter/setter pairs on the
//! object.  [`SearchParams`] collects the same knobs into one owned struct, so
//! an engine can be shared across threads without the setters racing.

use std::path::PathBuf;

use aln_core::{Alignment, Sequence, SubstMatrix};

use crate::Result;

/// How raw scores are reported — `SearchEngineI`'s `scoreMode`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ScoreMode {
    /// `SearchEngineI::basicScoreMode`.
    #[default]
    Basic,
    /// `SearchEngineI::complexityAdjustedScoreMode` — Phil Green's adjustment.
    /// See [`aln_core::stats::complexity_adjust`].
    ComplexityAdjusted,
}

/// Where an engine reads sequences from.
///
/// Search engines differ in what they can consume: rmblast wants a prepared
/// BLAST database, crossmatch wants FASTA on disk, and an in-process aligner
/// wants sequences in memory.  Making the source explicit lets an engine reject
/// what it cannot use up front instead of failing mid-search.
#[derive(Debug, Clone)]
pub enum SeqSource {
    /// FASTA file on disk.
    Fasta(PathBuf),
    /// UCSC 2bit file.
    TwoBit(PathBuf),
    /// Prepared BLAST database prefix (the path passed to `makeblastdb -out`).
    BlastDb(PathBuf),
    /// Sequences already in memory.
    Memory(Vec<Sequence>),
}

impl SeqSource {
    /// Number of sequences, when that is knowable without reading the source.
    pub fn len_hint(&self) -> Option<usize> {
        match self {
            SeqSource::Memory(v) => Some(v.len()),
            _ => None,
        }
    }

    pub fn path(&self) -> Option<&std::path::Path> {
        match self {
            SeqSource::Fasta(p) | SeqSource::TwoBit(p) | SeqSource::BlastDb(p) => Some(p),
            SeqSource::Memory(_) => None,
        }
    }
}

/// Search configuration — the fields of `SearchEngineI`, owned rather than set.
#[derive(Debug, Clone)]
pub struct SearchParams {
    /// Substitution matrix.  Engines that shell out write it to a temp file.
    pub matrix: Option<SubstMatrix>,

    /// Gap-open penalty (`setGapInit`).  Signed, as RepeatMasker writes it.
    pub gap_init: i32,
    /// Insertion gap-extension penalty (`setInsGapExt`).
    pub ins_gap_ext: i32,
    /// Deletion gap-extension penalty (`setDelGapExt`).
    pub del_gap_ext: i32,

    /// Minimum word/seed length (`setMinMatch`).
    pub min_match: u32,
    /// Minimum reportable score (`setMinScore`).
    pub min_score: i32,
    /// Banded-alignment half-width (`setBandwidth`); `None` for unbanded.
    pub bandwidth: Option<u32>,
    /// Overlap allowed between reported hits, as a percentage (`setMaskLevel`).
    pub mask_level: u32,
    /// Raw word size passed straight through (`setWordRaw`).
    pub word_raw: Option<u32>,

    pub score_mode: ScoreMode,
    /// Whether to produce alignment strings, not just coordinates
    /// (`setGenerateAlignments`).  Turning this off is markedly faster.
    pub generate_alignments: bool,

    /// Worker threads (`setCores`).  `None` lets the engine decide.
    pub cores: Option<usize>,
    /// Scratch directory for engines that shell out (`setTempDir`).
    pub temp_dir: Option<PathBuf>,
    /// Path to the external binary (`setPathToEngine`), where applicable.
    pub path_to_engine: Option<PathBuf>,
}

impl Default for SearchParams {
    fn default() -> Self {
        SearchParams {
            matrix: None,
            gap_init: -25,
            ins_gap_ext: -5,
            del_gap_ext: -5,
            min_match: 7,
            min_score: 150,
            bandwidth: None,
            mask_level: 80,
            word_raw: None,
            score_mode: ScoreMode::Basic,
            generate_alignments: true,
            cores: None,
            temp_dir: None,
            path_to_engine: None,
        }
    }
}

/// A seeded database search.
///
/// Implementors own their configuration; [`SearchParams`] is passed at
/// construction so `search` can take `&self` and run concurrently.
pub trait SearchEngine: Send + Sync {
    /// Engine name, for diagnostics and for `.out` provenance lines.
    fn name(&self) -> &'static str;

    /// Version string of the underlying engine — `SearchEngineI::getVersion`.
    fn version(&self) -> String;

    /// The configuration this engine was built with.
    fn params(&self) -> &SearchParams;

    /// Source kinds this engine can consume.  Checked before `search` so an
    /// unusable source fails immediately with a clear message.
    fn accepts(&self, source: &SeqSource) -> bool;

    /// Run the search and return every HSP that cleared the cutoffs.
    ///
    /// Under this crate's conventions the **query** is the genomic sequence and
    /// the **subject** is the consensus library — the same assignment
    /// RepeatMasker's matrices were built for.
    fn search(&self, query: &SeqSource, subject: &SeqSource) -> Result<Vec<Alignment>>;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn memory_sources_report_their_length() {
        let s = SeqSource::Memory(vec![Sequence::new("a", b"ACGT".to_vec())]);
        assert_eq!(s.len_hint(), Some(1));
        assert!(s.path().is_none());
    }

    #[test]
    fn file_sources_expose_a_path_but_no_length() {
        let s = SeqSource::Fasta(PathBuf::from("/tmp/x.fa"));
        assert_eq!(s.len_hint(), None);
        assert_eq!(s.path().unwrap().to_str(), Some("/tmp/x.fa"));
    }

    #[test]
    fn defaults_match_the_repeatmasker_baseline() {
        let p = SearchParams::default();
        assert_eq!(p.gap_init, -25);
        assert_eq!(p.min_score, 150);
        assert_eq!(p.mask_level, 80);
        assert_eq!(p.score_mode, ScoreMode::Basic);
        assert!(p.generate_alignments);
    }
}
