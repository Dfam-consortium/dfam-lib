//! Parallel alignment drivers — the replacement for GIRI's `ThreadedAligner`.
//!
//! `ThreadedAligner` hand-rolls a pthread pool with two mutexes and a global
//! `isMultithreaded` flag, and its `alignManyToMany` shares a single pair of
//! iterators across workers under lock.  These functions do the same work with
//! rayon: no global state, no hand-written synchronisation, and thread count set
//! by the caller's rayon pool rather than a static.
//!
//! Ordering is deterministic — results come back in input order regardless of
//! thread scheduling, which `ThreadedAligner` does not guarantee.

use aln_core::{Alignment, Sequence};
use rayon::prelude::*;

use crate::traits::PairwiseAligner;
use crate::Result;

/// Align every query against one subject, preparing the subject once.
///
/// Empty queries are skipped, as GIRI does.  Alignments that fail to clear
/// [`AlignParams::min_score`](crate::AlignParams::min_score) are dropped by the
/// backend and simply do not appear.
///
/// Results are in query order.
pub fn one_to_many<A: PairwiseAligner>(
    aligner: &A,
    subject: &Sequence,
    queries: &[Sequence],
) -> Result<Vec<Alignment>> {
    let prepared = aligner.prepare_subject(subject)?;
    let results: Vec<Result<Option<Alignment>>> = queries
        .par_iter()
        .map(|q| {
            if q.is_empty() {
                return Ok(None);
            }
            aligner.align_prepared(&prepared, q)
        })
        .collect();
    collect(results)
}

/// Align every sequence against every other, skipping self-comparisons.
///
/// GIRI skips self-alignment by pointer identity (`&ref == &*sptr`), which
/// silently fails to skip duplicated sequences that happen to be distinct
/// objects.  This skips by index, so `seqs[i]` is never aligned to itself, and
/// genuine duplicates are still compared.
///
/// This is the reference-selection pass in `autocons`.  Cost is `O(n^2)`
/// alignments; each subject is prepared once and reused across the inner loop.
///
/// Results are ordered by subject, then by query.
pub fn all_vs_all<A: PairwiseAligner>(
    aligner: &A,
    seqs: &[Sequence],
) -> Result<Vec<Alignment>> {
    let per_subject: Vec<Result<Vec<Alignment>>> = seqs
        .par_iter()
        .enumerate()
        .map(|(si, subject)| {
            if subject.is_empty() {
                return Ok(Vec::new());
            }
            let prepared = aligner.prepare_subject(subject)?;
            let mut out = Vec::new();
            for (qi, query) in seqs.iter().enumerate() {
                if qi == si || query.is_empty() {
                    continue;
                }
                if let Some(a) = aligner.align_prepared(&prepared, query)? {
                    out.push(a);
                }
            }
            Ok(out)
        })
        .collect();

    let mut all = Vec::new();
    for r in per_subject {
        all.extend(r?);
    }
    Ok(all)
}

/// Align explicit `(query, subject)` pairs — GIRI's `align_pairs`.
///
/// Each pair prepares its own subject, so this is the right driver only when
/// subjects genuinely differ; use [`one_to_many`] when they do not.
///
/// Results are in pair order.
pub fn align_pairs<A: PairwiseAligner>(
    aligner: &A,
    pairs: &[(Sequence, Sequence)],
) -> Result<Vec<Alignment>> {
    let results: Vec<Result<Option<Alignment>>> = pairs
        .par_iter()
        .map(|(query, subject)| {
            if query.is_empty() || subject.is_empty() {
                return Ok(None);
            }
            aligner.align(query, subject)
        })
        .collect();
    collect(results)
}

/// Sum of alignment scores — GIRI's drivers return this alongside the results.
pub fn total_score(alignments: &[Alignment]) -> i64 {
    alignments.iter().map(|a| a.score as i64).sum()
}

/// Flatten `Vec<Result<Option<_>>>`, propagating the first error.
fn collect(results: Vec<Result<Option<Alignment>>>) -> Result<Vec<Alignment>> {
    let mut out = Vec::with_capacity(results.len());
    for r in results {
        if let Some(a) = r? {
            out.push(a);
        }
    }
    Ok(out)
}
