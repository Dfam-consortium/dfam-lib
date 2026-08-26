//! An alignment plus the annotation needed to report it.
//!
//! [`Alignment`] deliberately carries only the geometric and scoring facts,
//! because `autocons` holds `O(n^2)` of them during its all-against-all pass.
//! Reporting needs more: percent divergence, library classification, a cluster
//! id. [`SearchResult`] is that fuller record — the analogue of RepeatMasker's
//! `SearchResult.pm` object, minus the search-engine plumbing.
//!
//! Nothing here recomputes anything. Build one from an [`Alignment`] and the
//! [`RescoreResult`] you already have.

use crate::align::Alignment;
use crate::stats::RescoreResult;

/// An alignment together with everything the output formats print.
#[derive(Debug, Clone)]
pub struct SearchResult {
    pub alignment: Alignment,

    /// Percent substitutions — the `.out` "perc div." column.
    pub pct_diverge: f64,
    /// Percent of the query deleted relative to the subject.
    pub pct_delete: f64,
    /// Percent inserted into the query relative to the subject.
    pub pct_insert: f64,

    /// Repeat class/family, e.g. `LINE/L1`.  The `.out` column after the name.
    pub subj_class: Option<String>,
    /// Cluster/annotation id linking fragments of one insertion.
    pub id: Option<u32>,
    /// Taxonomic lineage identifier, when the library carries one.
    pub lineage_id: Option<String>,
    /// `*` when a higher-scoring match overlaps this one.
    pub overlap: Option<char>,

    /// Matrix name, printed in the crossmatch alignment footer.
    pub matrix_name: Option<String>,
    /// CpG-modified Kimura divergence (percent).
    pub kimura_divergence: Option<f64>,
    /// Unadjusted Kimura divergence (percent).
    pub raw_kimura_divergence: Option<f64>,
    pub cpg_sites: Option<u32>,
}

impl SearchResult {
    /// Wrap an alignment with no annotation.
    pub fn new(alignment: Alignment) -> Self {
        SearchResult {
            alignment,
            pct_diverge: 0.0,
            pct_delete: 0.0,
            pct_insert: 0.0,
            subj_class: None,
            id: None,
            lineage_id: None,
            overlap: None,
            matrix_name: None,
            kimura_divergence: None,
            raw_kimura_divergence: None,
            cpg_sites: None,
        }
    }

    /// Build from an alignment and a rescoring pass.
    ///
    /// `pct_diverge` is taken as substitutions over aligned columns, matching
    /// the `.out` "perc div." column — *not* the Kimura divergence, which is
    /// reported separately.
    pub fn from_rescore(alignment: Alignment, r: &RescoreResult) -> Self {
        let d = &r.divergence;
        let denom = d.well_characterized as f64;
        let pct_diverge = if denom > 0.0 {
            ((d.transitions + d.transversions as f64) / denom) * 100.0
        } else {
            0.0
        };
        SearchResult {
            alignment,
            pct_diverge,
            pct_delete: r.pct_delete,
            pct_insert: r.pct_insert,
            subj_class: None,
            id: None,
            lineage_id: None,
            overlap: None,
            matrix_name: None,
            kimura_divergence: d.value,
            raw_kimura_divergence: None,
            cpg_sites: Some(d.cpg_sites),
        }
    }

    /// Bases of the query beyond the alignment's end, or 0 if unknown.
    pub fn query_left(&self) -> u64 {
        self.alignment.query_remaining().unwrap_or(0) as u64
    }

    /// Bases of the subject beyond the alignment's end, or 0 if unknown.
    pub fn subj_left(&self) -> u64 {
        self.alignment.subj_remaining().unwrap_or(0) as u64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::{EditOp, EditScript};
    use crate::Strand;

    fn aln() -> Alignment {
        let mut e = EditScript::new();
        e.push(EditOp::Sub, 10);
        let mut a = Alignment::new("q", "s", 5, 2, Strand::Plus, 100, e);
        a.query_len = Some(40);
        a.subj_len = Some(30);
        a
    }

    #[test]
    fn remaining_bases_come_from_the_alignment() {
        let r = SearchResult::new(aln());
        assert_eq!(r.query_left(), 40 - 15);
        assert_eq!(r.subj_left(), 30 - 12);
    }

    #[test]
    fn unknown_lengths_report_zero_remaining() {
        let mut a = aln();
        a.query_len = None;
        a.subj_len = None;
        let r = SearchResult::new(a);
        assert_eq!(r.query_left(), 0);
        assert_eq!(r.subj_left(), 0);
    }
}
