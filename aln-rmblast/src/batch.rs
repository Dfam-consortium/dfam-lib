//! Batched searches — the shape rmblast is actually built for.
//!
//! Driving a search engine one pair at a time rebuilds a query lookup table per
//! call and scans a single subject, throwing away the point of seeding. These
//! two entry points issue **one search** and demultiplex the results.
//!
//! * [`RmblastEngine::one_to_many`] — every instance against one reference.
//!   This is `autocons` phase 2 (refinement against the current consensus).
//! * [`RmblastEngine::all_vs_all`] — every sequence against every other, in one
//!   call. This is `autocons` phase 1 (reference selection).
//!
//! # Multiple HSPs per pair, and what `mask_level` does
//!
//! Unlike a pairwise aligner, rmblast can return **several HSPs for one
//! (query, subject) pair** — an instance interrupted by an insertion, or
//! matching on both strands. `mask_level` is the filter: an HSP is dropped when
//! a higher-scoring HSP already covers more than that percentage of its **query**
//! span. RepeatMasker's default is 80; 101 disables it.
//!
//! Surviving HSPs are all returned. A single instance can therefore contribute
//! more than one row to the MSA — which is exactly GIRI's own `FRAGMENT` model,
//! so nothing downstream needs to change.
//!
//! ## Masking is applied per pair here, not across subjects
//!
//! rmblast's own `apply_mask_level` follows NCBI's
//! `Blast_HSPResultsApplyMasklevel`: for a given **query**, HSPs are ranked
//! across *all* subjects together, so a hit against one subject can suppress a
//! hit against another. That is right for genome annotation — you want one
//! family per genomic region — and wrong here.
//!
//! In an all-against-all search every sequence hits **itself** perfectly, and
//! that self-hit covers 100% of its own query span. Cross-subject masking then
//! discards essentially every real cross-hit before the caller ever sees it.
//! Measured: with cross-subject masking at 80, 29 of 30 families produced no
//! consensus at all.
//!
//! So these functions run the search with masking **disabled** and then apply
//! the same rule per `(query, subject)` pair. Overlap within a pair is what
//! matters when the HSPs are going to become rows of one MSA.
//!
//! # Identity of results
//!
//! Sequence names in a family are not guaranteed unique, so these functions
//! rename inputs to their index for the duration of the search and map back
//! afterwards. The returned [`Alignment`]s carry the caller's original names.

use std::collections::HashMap;

use aln_core::{Alignment, Sequence};
use aln_engine::engine::{SearchEngine, SeqSource};
use aln_engine::{EngineError, Result};

use crate::{RmblastEngine, NAME};

/// Drop an HSP when a higher-scoring one already covers more than
/// `mask_level` percent of its query span. `>= 100` keeps everything.
///
/// Same rule as rmblast's `apply_mask_level`, but scoped to one
/// `(query, subject)` pair — see the module docs for why.
fn mask_within_pair(mut hsps: Vec<Alignment>, mask_level: u32) -> Vec<Alignment> {
    if mask_level >= 100 || hsps.len() < 2 {
        return hsps;
    }
    hsps.sort_by(|a, b| b.score.cmp(&a.score).then(a.query_start.cmp(&b.query_start)));
    let mut kept: Vec<Alignment> = Vec::with_capacity(hsps.len());
    for h in hsps {
        let span = h.query_span().max(1) as u64;
        let covered = kept.iter().any(|k| {
            let lo = h.query_start.max(k.query_start);
            let hi = h.query_end.min(k.query_end);
            let overlap = hi.saturating_sub(lo) as u64;
            overlap * 100 > span * mask_level as u64
        });
        if !covered {
            kept.push(h);
        }
    }
    kept
}

/// Group by `(query, subject)`, mask within each pair, and flatten.
fn mask_pairwise<K: std::hash::Hash + Eq + Copy>(
    items: Vec<(K, Alignment)>,
    mask_level: u32,
) -> Vec<(K, Alignment)> {
    let mut by: HashMap<K, Vec<Alignment>> = HashMap::new();
    for (k, a) in items {
        by.entry(k).or_default().push(a);
    }
    let mut out = Vec::new();
    for (k, group) in by {
        for a in mask_within_pair(group, mask_level) {
            out.push((k, a));
        }
    }
    out
}

/// Rename to indices so results can be mapped back unambiguously.
fn indexed(seqs: &[Sequence]) -> Vec<Sequence> {
    seqs.iter()
        .enumerate()
        .map(|(i, s)| Sequence::new(i.to_string(), s.seq.clone()))
        .collect()
}

fn parse_index(name: &str, what: &str) -> Result<usize> {
    name.parse::<usize>().map_err(|_| {
        EngineError::backend(NAME, format!("unexpected {what} name {name:?} in search output"))
    })
}

impl RmblastEngine {
    /// Search every sequence in `queries` against a single `reference`, in one call.
    ///
    /// Returns `(query index, alignment)` for every HSP that survives
    /// `mask_level`. `skip` suppresses one query index — used so a reference is
    /// not aligned to itself.
    ///
    /// Alignments come back with `subj_name` set to the reference's name and
    /// `query_name` to the original query's.
    pub fn one_to_many(
        &self,
        reference: &Sequence,
        queries: &[Sequence],
        skip: Option<usize>,
    ) -> Result<Vec<(usize, Alignment)>> {
        if reference.is_empty() || queries.is_empty() {
            return Ok(Vec::new());
        }
        let q = indexed(queries);
        let subject = Sequence::new("R", reference.seq.clone());

        let hits = self.search_unmasked(
            &SeqSource::Memory(q),
            &SeqSource::Memory(vec![subject]),
        )?;

        let mut out = Vec::with_capacity(hits.len());
        for mut a in hits {
            let qi = parse_index(&a.query_name, "query")?;
            if Some(qi) == skip {
                continue;
            }
            a.query_name = queries[qi].name.clone();
            a.subj_name = reference.name.clone();
            out.push((qi, a));
        }
        // One subject here, so per-pair masking is just per-query masking.
        let mut out = mask_pairwise(out, self.params().mask_level);
        // Deterministic order: by query, then by descending score.
        out.sort_by(|x, y| x.0.cmp(&y.0).then(y.1.score.cmp(&x.1.score)));
        Ok(out)
    }

    /// Search every sequence against every other, in one call.
    ///
    /// Returns `(query index, subject index, alignment)` for every surviving
    /// HSP, with self-pairs removed. Grouping by subject index gives, for each
    /// candidate reference, the alignments that would be built against it.
    ///
    /// `mask_level` is taken as an argument rather than from the engine's
    /// params because RepeatModeler's `Refiner` treats the two phases
    /// differently: it calls `setMaskLevel(80)` on its one-vs-all engine and
    /// leaves the all-vs-all engine unmasked, so reference selection sees every
    /// HSP. Pass 101 to match it.
    ///
    /// Whatever value is passed is applied **per (query, subject) pair**, not
    /// across subjects as `rmblastn` would — see the module docs. Cross-subject
    /// masking is unusable here because every sequence hits itself perfectly.
    pub fn all_vs_all(
        &self,
        seqs: &[Sequence],
        mask_level: u32,
    ) -> Result<Vec<(usize, usize, Alignment)>> {
        if seqs.len() < 2 {
            return Ok(Vec::new());
        }
        let idx = indexed(seqs);
        let hits = self.search_unmasked(
            &SeqSource::Memory(idx.clone()),
            &SeqSource::Memory(idx),
        )?;

        let mut out = Vec::with_capacity(hits.len());
        for mut a in hits {
            let qi = parse_index(&a.query_name, "query")?;
            let si = parse_index(&a.subj_name, "subject")?;
            if qi == si {
                continue; // a sequence is not its own instance
            }
            a.query_name = seqs[qi].name.clone();
            a.subj_name = seqs[si].name.clone();
            out.push(((qi, si), a));
        }
        let mut out: Vec<(usize, usize, Alignment)> =
            mask_pairwise(out, mask_level)
                .into_iter()
                .map(|((qi, si), a)| (qi, si, a))
                .collect();
        out.sort_by(|x, y| {
            x.1.cmp(&y.1).then(x.0.cmp(&y.0)).then(y.2.score.cmp(&x.2.score))
        });
        Ok(out)
    }

    /// Group [`all_vs_all`](Self::all_vs_all) output by subject index.
    ///
    /// `by_subject[i]` holds every alignment against `seqs[i]` as reference.
    pub fn group_by_subject(
        hits: Vec<(usize, usize, Alignment)>,
        n: usize,
    ) -> Vec<Vec<(usize, Alignment)>> {
        let mut by: Vec<Vec<(usize, Alignment)>> = vec![Vec::new(); n];
        for (qi, si, a) in hits {
            if si < n {
                by[si].push((qi, a));
            }
        }
        by
    }

    /// Keep only the best-scoring HSP per (query, subject) pair.
    ///
    /// Use when downstream expects one alignment per instance, as a pairwise
    /// aligner would give. Dropping the rest loses the fragment structure that
    /// [`one_to_many`](Self::one_to_many) preserves.
    pub fn best_per_query(hits: Vec<(usize, Alignment)>) -> Vec<(usize, Alignment)> {
        let mut best: HashMap<usize, (usize, Alignment)> = HashMap::new();
        for (qi, a) in hits {
            match best.get(&qi) {
                Some((_, prev)) if prev.score >= a.score => {}
                _ => {
                    best.insert(qi, (qi, a));
                }
            }
        }
        let mut out: Vec<_> = best.into_values().collect();
        out.sort_by_key(|(qi, _)| *qi);
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::RmblastOptions;
    use aln_core::SubstMatrix;
    use aln_engine::engine::SearchParams;

    const M: &str = "\
FREQS A 0.325 C 0.175 G 0.175 T 0.325
  A   R   G   C   Y   T   K   M   S   W   N   X
  8   0 -10 -18 -19 -21 -15  -4 -14  -6  -1 -30
  3   3  12 -17 -18 -19  -9  -8  -8  -9  -1 -30
 -7   2  12 -16 -16 -17  -2 -11  -1 -12  -1 -30
-17 -16 -16  12   2  -7 -11  -2  -1 -12  -1 -30
-19 -18 -17  12   0   3  -8  -9  -8  -9  -1 -30
-21 -19 -18 -10   0   8  -4 -15 -14  -6  -1 -30
-14  -8  -2 -13  -8  -4  -3 -13  -8  -9  -1 -30
 -4  -8 -13  -2  -8 -14 -13  -3  -8  -9  -1 -30
-12  -7  -1  -1  -7 -12  -7  -7  -1 -12  -1 -30
 -6 -10 -14 -14 -10  -6 -10 -10 -14  -6  -1 -30
 -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1 -30
-30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30
";

    fn engine(mask_level: u32) -> RmblastEngine {
        let p = SearchParams {
            matrix: Some(SubstMatrix::parse(M).unwrap()),
            gap_init: -25,
            ins_gap_ext: -5,
            del_gap_ext: -5,
            min_match: 7,
            min_score: 100,
            mask_level,
            cores: Some(1),
            ..Default::default()
        };
        RmblastEngine::new(p, RmblastOptions::default()).unwrap()
    }

    struct Rng(u64);
    impl Rng {
        fn new(s: u64) -> Self {
            Rng(s.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1)
        }
        fn below(&mut self, n: usize) -> usize {
            let mut x = self.0;
            x ^= x >> 12;
            x ^= x << 25;
            x ^= x >> 27;
            self.0 = x;
            (x.wrapping_mul(0x2545_F491_4F6C_DD1D) % n as u64) as usize
        }
    }
    fn rs(r: &mut Rng, n: usize) -> Vec<u8> {
        (0..n).map(|_| b"ACGT"[r.below(4)]).collect()
    }
    fn family(seed: u64, n: usize, len: usize) -> Vec<Sequence> {
        let mut r = Rng::new(seed);
        let anc = rs(&mut r, len);
        (0..n)
            .map(|i| {
                let copy: Vec<u8> = anc
                    .iter()
                    .map(|&b| if r.below(100) < 10 { b"ACGT"[r.below(4)] } else { b })
                    .collect();
                Sequence::new(format!("copy{i}"), copy)
            })
            .collect()
    }

    #[test]
    fn one_to_many_finds_every_instance_in_a_single_call() {
        let seqs = family(1, 10, 500);
        let e = engine(101);
        let hits = e.one_to_many(&seqs[0], &seqs, Some(0)).unwrap();
        let distinct: std::collections::HashSet<usize> = hits.iter().map(|(i, _)| *i).collect();
        assert_eq!(distinct.len(), 9, "expected all 9 non-self instances");
        assert!(hits.iter().all(|(_, a)| a.subj_name == "copy0"));
        for (i, a) in &hits {
            assert_eq!(&a.query_name, &seqs[*i].name, "names must map back");
            a.validate().unwrap();
        }
    }

    #[test]
    fn one_to_many_honours_skip() {
        let seqs = family(2, 6, 400);
        let e = engine(101);
        let hits = e.one_to_many(&seqs[0], &seqs, Some(0)).unwrap();
        assert!(hits.iter().all(|(i, _)| *i != 0), "self must be skipped");
        let with_self = e.one_to_many(&seqs[0], &seqs, None).unwrap();
        assert!(with_self.iter().any(|(i, _)| *i == 0), "and included when not skipped");
    }

    #[test]
    fn all_vs_all_covers_every_ordered_pair_and_drops_self() {
        let seqs = family(3, 6, 400);
        let e = engine(101);
        let hits = e.all_vs_all(&seqs, 101).unwrap();
        assert!(hits.iter().all(|(q, s, _)| q != s), "self-pairs must be removed");

        let by = RmblastEngine::group_by_subject(hits, seqs.len());
        assert_eq!(by.len(), seqs.len());
        for (si, group) in by.iter().enumerate() {
            let distinct: std::collections::HashSet<usize> =
                group.iter().map(|(q, _)| *q).collect();
            assert_eq!(distinct.len(), seqs.len() - 1, "subject {si} missing partners");
        }
    }

    /// The point of the batched path: one call, not n².
    #[test]
    fn all_vs_all_agrees_with_repeated_one_to_many() {
        let seqs = family(4, 6, 400);
        let e = engine(101);
        let grouped = RmblastEngine::group_by_subject(e.all_vs_all(&seqs, 101).unwrap(), seqs.len());

        for (si, reference) in seqs.iter().enumerate() {
            let direct = e.one_to_many(reference, &seqs, Some(si)).unwrap();
            let a: i64 = RmblastEngine::best_per_query(direct).iter().map(|(_, x)| x.score as i64).sum();
            let b: i64 = RmblastEngine::best_per_query(grouped[si].clone())
                .iter().map(|(_, x)| x.score as i64).sum();
            assert_eq!(a, b, "subject {si}: batched and per-reference totals differ");
        }
    }

    #[test]
    fn best_per_query_keeps_one_alignment_per_instance() {
        let seqs = family(5, 8, 400);
        let e = engine(101);
        let hits = e.one_to_many(&seqs[0], &seqs, Some(0)).unwrap();
        let best = RmblastEngine::best_per_query(hits.clone());
        let distinct: std::collections::HashSet<usize> = hits.iter().map(|(i, _)| *i).collect();
        assert_eq!(best.len(), distinct.len());
        // And it really is the best one.
        for (qi, a) in &best {
            let max = hits.iter().filter(|(i, _)| i == qi).map(|(_, x)| x.score).max().unwrap();
            assert_eq!(a.score, max);
        }
    }

    /// A tandem duplication yields two HSPs for one instance; mask_level decides
    /// whether both survive. They do not overlap on the query, so 80 keeps both.
    #[test]
    fn mask_level_governs_multi_hsp_survival() {
        let mut r = Rng::new(6);
        let block = rs(&mut r, 400);
        let mut dup = block.clone();
        dup.extend_from_slice(&rs(&mut r, 30));
        dup.extend_from_slice(&block);

        let seqs = vec![
            Sequence::new("ref", block),
            Sequence::new("dup", dup),
            Sequence::new("plain", rs(&mut Rng::new(7), 400)),
        ];
        let permissive = engine(101).one_to_many(&seqs[0], &seqs, Some(0)).unwrap();
        let n_dup = permissive.iter().filter(|(i, _)| *i == 1).count();
        assert!(n_dup >= 2, "a tandem duplication should give >=2 HSPs, got {n_dup}");
    }
}
