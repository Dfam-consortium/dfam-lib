//! rmblast backend — a [`SearchEngine`] over the Rust port of rmblastn.
//!
//! This is the *seeded database search* side of the stack, the analogue of
//! RepeatMasker's `NCBIBlastSearchEngine`.  Where [`PairwiseAligner`] backends
//! (`aln-parasail`, `aln-reference`) do full dynamic programming on two
//! sequences, this seeds, extends, filters by mask level, and reports many HSPs
//! across a database.
//!
//! [`PairwiseAligner`]: aln_engine::PairwiseAligner
//!
//! # Coordinates line up already
//!
//! `rmblast-lib`'s `Hsp` uses 0-based half-open offsets with subject coordinates
//! always on the plus strand — the same convention [`aln_core::Alignment`]
//! adopted.  Conversion is therefore a field copy, not a re-basing, which is
//! most of why this backend is thin.  The one thing that does need care is the
//! scoring matrix; see [`matrix`].
//!
//! # Sentinels
//!
//! rmblast works in BLASTNA encoding with a sentinel byte at each end of every
//! sequence (`encode_iupac`).  Offsets it reports exclude the sentinels, so no
//! adjustment is needed on the way out — `offsets_exclude_sentinels` pins that.

pub mod batch;
pub mod matrix;
pub mod pairwise;

use aln_core::align::{Alignment, EditOp, EditScript};
use aln_core::{Sequence, Strand};
use aln_engine::engine::ScoreMode;
use aln_engine::{EngineError, Result, SearchEngine, SearchParams, SeqSource};
use rayon::prelude::*;

use rmblast_lib::encoding::encode_iupac;
use rmblast_lib::hits::{EditOp as RmEditOp, Strand as RmStrand};
use rmblast_lib::matrix::ScoreMatrix;
use rmblast_lib::options::{MtMode, SearchParams as RmParams, SeedMode};
use rmblast_lib::output::AlignResult;
use rmblast_lib::search::{apply_mask_level, build_query_lookup, search_with_query_lookup};
use rmblast_lib::seq::SubjectDb;

pub use pairwise::{PreparedSubject, RmblastPairwise};

pub(crate) const NAME: &str = "rmblast";

// ── rmblast-specific options ──────────────────────────────────────────────────

/// Knobs rmblast has that [`SearchParams`] does not.
///
/// Defaults reproduce the invocation `dfam-curator`'s `run_rmblastn` uses:
/// DUST off, complexity adjustment on, and X-drops derived from the minimum
/// score.
#[derive(Debug, Clone)]
pub struct RmblastOptions {
    /// Ungapped X-drop.  `None` derives `min_score * 2`.
    pub xdrop_ungap: Option<i32>,
    /// Gapped X-drop.  `None` derives `min_score / 2`.
    pub xdrop_gap: Option<i32>,
    /// Final gapped X-drop.  `None` derives `min_score`.
    pub xdrop_gap_final: Option<i32>,
    /// Low-complexity query masking.  RepeatMasker runs `-dust no`.
    pub dust: bool,
    pub seed_mode: SeedMode,
    pub mt_mode: MtMode,
    /// Karlin-Altschul ungapped cutoff; `None` uses rmblast's fallback.
    pub ungapped_cutoff: Option<i32>,
}

impl Default for RmblastOptions {
    fn default() -> Self {
        RmblastOptions {
            xdrop_ungap: None,
            xdrop_gap: None,
            xdrop_gap_final: None,
            dust: false,
            seed_mode: SeedMode::Combined,
            mt_mode: MtMode::SplitByDb,
            ungapped_cutoff: None,
        }
    }
}

// ── Engine ────────────────────────────────────────────────────────────────────

/// A configured rmblast search.
pub struct RmblastEngine {
    params: SearchParams,
    rm: RmParams,
    matrix: ScoreMatrix,
    /// Sized from `SearchParams::cores`, used to spread a multi-query search.
    /// `None` for a single core, where a pool would only add overhead.
    pool: Option<rayon::ThreadPool>,
}

impl RmblastEngine {
    /// Build from generic [`SearchParams`], converting the matrix it carries.
    ///
    /// Fails if no matrix is set: rmblast has no usable default for nucleotide
    /// searches with RepeatMasker-style scoring.
    pub fn new(params: SearchParams, opts: RmblastOptions) -> Result<Self> {
        let subst = params.matrix.as_ref().ok_or_else(|| {
            EngineError::Params("rmblast needs a scoring matrix; SearchParams::matrix is None".into())
        })?;
        let m = matrix::to_rmblast(subst)?;
        Self::assemble(params, opts, m)
    }

    /// Build using a matrix file read by rmblast directly.
    ///
    /// Preferred when reproducing a real `rmblastn` run: the file is parsed by
    /// rmblast's own reader with no conversion in between.  Note the file must
    /// be in **NCBI** layout (`Matrices/ncbi/nt/...`), not crossmatch layout —
    /// the two are transposes; see [`matrix`].
    pub fn with_matrix_file(
        params: SearchParams,
        opts: RmblastOptions,
        path: &str,
    ) -> Result<Self> {
        let m = ScoreMatrix::from_file(path).map_err(|e| {
            EngineError::backend(NAME, format!("cannot read matrix {path}: {e}"))
        })?;
        Self::assemble(params, opts, m)
    }

    fn assemble(params: SearchParams, opts: RmblastOptions, matrix: ScoreMatrix) -> Result<Self> {
        if params.ins_gap_ext != params.del_gap_ext {
            return Err(EngineError::unsupported(
                NAME,
                format!(
                    "rmblast has a single gap-extension penalty, but insertion ({}) and \
                     deletion ({}) differ",
                    params.ins_gap_ext, params.del_gap_ext
                ),
            ));
        }
        let min_score = params.min_score.max(0);
        let gap_extend = params.ins_gap_ext.unsigned_abs() as i32;
        // NCBI and crossmatch cost a length-k gap differently:
        //
        //   crossmatch / aln_core::stats::rescore :  open + (k-1) * extend
        //   NCBI / rmblast                        :  open +  k    * extend
        //
        // so NCBI's open must be reduced by one extension to describe the same
        // scoring system.  This is why RepeatMasker passes `-gapopen 20` for a
        // crossmatch `gap_init` of -25 with `gap_ext` of -5, and why
        // `dfam-curator`'s BlastParams carries the same 20/5 pair.
        let gap_open = (params.gap_init.unsigned_abs() as i32)
            .checked_sub(gap_extend)
            .filter(|v| *v >= 0)
            .ok_or_else(|| {
                EngineError::Params(format!(
                    "gap_init ({}) must be at least as large in magnitude as the gap \
                     extension ({}); NCBI's open cost excludes the first position",
                    params.gap_init, params.ins_gap_ext
                ))
            })?;
        let rm = RmParams {
            gap_open,
            gap_extend,
            matrix_name: matrix.name.clone(),
            word_size: params.word_raw.unwrap_or(params.min_match) as usize,
            xdrop_ungap: opts.xdrop_ungap.unwrap_or(min_score * 2),
            xdrop_gap: opts.xdrop_gap.unwrap_or(min_score / 2),
            xdrop_gap_final: opts.xdrop_gap_final.unwrap_or(min_score),
            min_raw_gapped_score: min_score,
            ungapped_cutoff: opts.ungapped_cutoff,
            complexity_adjust: params.score_mode == ScoreMode::ComplexityAdjusted,
            dust: opts.dust,
            mask_level: params.mask_level,
            num_threads: params.cores.unwrap_or(1).max(1),
            mt_mode: opts.mt_mode,
            seed_mode: opts.seed_mode,
        };
        let pool = (rm.num_threads > 1)
            .then(|| {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(rm.num_threads)
                    .build()
                    .map_err(|e| EngineError::backend(NAME, format!("thread pool: {e}")))
            })
            .transpose()?;
        Ok(RmblastEngine { params, rm, matrix, pool })
    }

    /// The rmblast-side parameters, for diagnostics and for reproducing a run on
    /// the `rmblastn` command line.
    pub fn rmblast_params(&self) -> &RmParams {
        &self.rm
    }

    /// Search and return rmblast's own per-HSP statistics alongside each
    /// alignment.
    ///
    /// `AlignStats` carries Kimura divergence, CpG counts and the percent
    /// substitution/gap figures that RepeatMasker's `.out` writer needs.  They
    /// are computed during the search, so taking them here is free compared with
    /// recomputing via [`aln_core::stats`].
    pub fn search_with_stats(
        &self,
        query: &SeqSource,
        subject: &SeqSource,
    ) -> Result<Vec<(Alignment, rmblast_lib::stats::AlignStats)>> {
        let results = self.run(query, subject)?;
        results
            .into_iter()
            .map(|r| {
                let stats = r.stats.clone();
                to_alignment(&r).map(|a| (a, stats))
            })
            .collect()
    }

    /// Search with cross-subject masking disabled; see [`Self::run_masked`].
    pub(crate) fn search_unmasked(
        &self,
        query: &SeqSource,
        subject: &SeqSource,
    ) -> Result<Vec<Alignment>> {
        self.run_masked(query, subject, 101)?
            .iter()
            .map(to_alignment)
            .collect()
    }

    fn run(&self, query: &SeqSource, subject: &SeqSource) -> Result<Vec<AlignResult>> {
        self.run_masked(query, subject, self.rm.mask_level)
    }

    /// Search with an explicit mask level rather than the configured one.
    ///
    /// `apply_mask_level` ranks a query's HSPs across *all* subjects at once
    /// (NCBI's `Blast_HSPResultsApplyMasklevel`).  The batched entry points in
    /// [`crate::batch`] need that suppressed — in an all-against-all search a
    /// sequence's perfect self-hit would otherwise mask every real cross-hit —
    /// so they pass 101 here and filter per pair themselves.
    fn run_masked(
        &self,
        query: &SeqSource,
        subject: &SeqSource,
        mask_level: u32,
    ) -> Result<Vec<AlignResult>> {
        let queries = self.load_queries(query)?;
        let mut out = Vec::new();

        match subject {
            SeqSource::Memory(subjects) => {
                let names: Vec<String> = subjects.iter().map(|s| s.name.clone()).collect();
                let encoded: Vec<Vec<u8>> =
                    subjects.iter().map(|s| encode_iupac(&s.seq)).collect();
                // Queries are independent — each builds its own lookup table and
                // scans every subject — so this is the axis to spread work on.
                let one = |q: &Sequence| {
                    let eq = encode_iupac(&q.seq);
                    let (lookup, _) = build_query_lookup(&eq, &self.rm);
                    let mut per_query: Vec<AlignResult> = Vec::new();
                    for (s, name) in encoded.iter().zip(&names) {
                        per_query.extend(search_with_query_lookup(
                            // `prepared`: rmblast-lib gained this 10th parameter
                            // on 2026-08-15 for reusing a precomputed subject
                            // across searches. `None` reproduces the previous
                            // behaviour exactly. If that parameter exists so this
                            // all-vs-all loop can reuse subjects, wiring it here
                            // is the optimisation — this is a compile fix, not a
                            // considered answer to that question.
                            &lookup, &eq, &q.name, s, name, &self.rm, &self.matrix, 0, &[],
                            None,
                        ));
                    }
                    apply_mask_level(&mut per_query, mask_level, &names);
                    per_query
                };
                match &self.pool {
                    Some(pool) if queries.len() > 1 => {
                        let per: Vec<Vec<AlignResult>> =
                            pool.install(|| queries.par_iter().map(one).collect());
                        out.extend(per.into_iter().flatten());
                    }
                    _ => out.extend(queries.iter().flat_map(one)),
                }
            }
            SeqSource::Fasta(p) | SeqSource::TwoBit(p) => {
                let path = p.to_str().ok_or_else(|| {
                    EngineError::backend(NAME, format!("non-UTF8 database path {}", p.display()))
                })?;
                let db = SubjectDb::open(path).map_err(|e| {
                    EngineError::backend(NAME, format!("cannot open database {path}: {e}"))
                })?;
                let names: Vec<String> =
                    db.sequences().iter().map(|s| s.name.clone()).collect();
                for q in &queries {
                    let eq = encode_iupac(&q.seq);
                    let (lookup, _) = build_query_lookup(&eq, &self.rm);
                    let mut per_query: Vec<AlignResult> = Vec::new();
                    for name in &names {
                        let seq = db.get_full_sequence_blastna(name).map_err(|e| {
                            EngineError::backend(NAME, format!("reading {name}: {e}"))
                        })?;
                        let n_mask = db.get_n_mask(name);
                        per_query.extend(search_with_query_lookup(
                            &lookup, &eq, &q.name, &seq, name, &self.rm, &self.matrix, 0, &n_mask,
                            None,  // see the note on the in-memory path above
                        ));
                    }
                    apply_mask_level(&mut per_query, mask_level, &names);
                    out.extend(per_query);
                }
            }
            SeqSource::BlastDb(p) => {
                return Err(EngineError::unsupported(
                    NAME,
                    format!(
                        "this backend is the in-process rmblast port, which reads FASTA \
                         and 2bit directly; a prepared BLAST database ({}) would need the \
                         external rmblastn binary",
                        p.display()
                    ),
                ))
            }
        }
        Ok(out)
    }

    fn load_queries(&self, source: &SeqSource) -> Result<Vec<Sequence>> {
        match source {
            SeqSource::Memory(v) => Ok(v.clone()),
            SeqSource::Fasta(p) => {
                let f = std::fs::File::open(p)?;
                let mut reader =
                    rmblast_lib::seq::fasta::FastaReader::new(std::io::BufReader::new(f));
                let mut out = Vec::new();
                while let Some(rec) = reader.next_record().map_err(|e| {
                    EngineError::backend(NAME, format!("reading {}: {e}", p.display()))
                })? {
                    out.push(Sequence::new(rec.id.clone(), rec.bases().to_vec()));
                }
                Ok(out)
            }
            other => Err(EngineError::unsupported(
                NAME,
                format!("queries cannot be read from {other:?}"),
            )),
        }
    }
}

impl SearchEngine for RmblastEngine {
    fn name(&self) -> &'static str {
        NAME
    }

    fn version(&self) -> String {
        format!("rmblast-lib {} (in-process)", env!("CARGO_PKG_VERSION"))
    }

    fn params(&self) -> &SearchParams {
        &self.params
    }

    fn accepts(&self, source: &SeqSource) -> bool {
        matches!(
            source,
            SeqSource::Memory(_) | SeqSource::Fasta(_) | SeqSource::TwoBit(_)
        )
    }

    fn search(&self, query: &SeqSource, subject: &SeqSource) -> Result<Vec<Alignment>> {
        self.run(query, subject)?.iter().map(to_alignment).collect()
    }
}

// ── Conversion ────────────────────────────────────────────────────────────────

/// Convert one rmblast HSP into an [`Alignment`].
///
/// A field copy: both sides use 0-based half-open offsets with plus-strand
/// subject coordinates, and both name their edit operations from the same
/// perspective (`GapInQuery` means the query row shows `-`).
pub fn to_alignment(r: &AlignResult) -> Result<Alignment> {
    let h = &r.hsp;
    let mut edits = EditScript::new();
    for &(op, n) in &h.edit_script.ops {
        let mapped = match op {
            RmEditOp::Sub => EditOp::Sub,
            RmEditOp::GapInQuery => EditOp::GapInQuery,
            RmEditOp::GapInSubject => EditOp::GapInSubject,
        };
        edits.push(mapped, n);
    }

    let mut a = Alignment::new(
        r.query_id.clone(),
        r.subject_id.clone(),
        h.q_start as usize,
        h.s_start as usize,
        match h.strand {
            RmStrand::Plus => Strand::Plus,
            RmStrand::Minus => Strand::Minus,
        },
        h.score,
        edits,
    );
    a.query_len = Some(h.q_len as usize);
    a.subj_len = Some(h.s_len as usize);

    // rmblast reports spans independently of the edit script, so disagreement
    // means a convention has drifted — catch it here rather than downstream.
    if a.query_end != h.q_end as usize || a.subj_end != h.s_end as usize {
        return Err(EngineError::backend(
            NAME,
            format!(
                "rmblast HSP is internally inconsistent: reported query {}..{} / subject \
                 {}..{}, but the edit script implies query {}..{} / subject {}..{}",
                h.q_start, h.q_end, h.s_start, h.s_end,
                a.query_start, a.query_end, a.subj_start, a.subj_end
            ),
        ));
    }
    a.validate()
        .map_err(|e| EngineError::backend(NAME, format!("rmblast HSP failed validation: {e}")))?;
    Ok(a)
}
