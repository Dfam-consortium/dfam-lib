//! Parity tests for [`aln_core::stats::rescore`] against RepeatMasker's
//! `SearchResult.pm::rescoreAlignment`.
//!
//! The expected values below were produced by running the Perl directly against
//! `Matrices/crossmatch/14p35g.matrix` with `gapOpenPenalty => -25` and
//! `gapExtPenalty => -5`.  To regenerate:
//!
//! ```perl
//! use lib "/usr/local/RepeatMasker";
//! use SearchResult; use Matrix;
//! my $mat = Matrix->new(fileName =>
//!     "/usr/local/RepeatMasker/Matrices/crossmatch/14p35g.matrix");
//! my $r = SearchResult->new(
//!     queryName => "q", subjName => "s",
//!     queryStart => 1, queryEnd => $ungapped_q,
//!     subjStart  => 1, subjEnd  => $ungapped_s,
//!     orientation => "", score => 0,
//!     queryString => $q, subjString => $s);
//! my ($score, $kimura, $cpg, $pi, $pd, $ps, $xf, $wcb, $ti, $tv) =
//!     $r->rescoreAlignment(scoreMatrix => $mat,
//!                          gapOpenPenalty => -25, gapExtPenalty => -5,
//!                          # plus scoreCpGMod / divCpGMod / complexityAdjust
//!                         );
//! ```
//!
//! These are self-contained — the matrix is inlined, so they do not need a
//! RepeatMasker installation to run.

use aln_core::stats::{rescore, RescoreParams};
use aln_core::SubstMatrix;

const M14P35G: &str = "\
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

/// One row of Perl output.
struct Case {
    query: &'static str,
    subject: &'static str,
    score_cpg_mod: bool,
    div_cpg_mod: bool,
    complexity_adjust: bool,
    /// `None` where the Perl returned its `100.00` saturation initialiser.
    kimura: Option<f64>,
    score: i32,
    cpg_sites: u32,
    pct_insert: f64,
    pct_delete: f64,
    well_characterized: u32,
    transitions: f64,
    transversions: u32,
}

const CASES: &[Case] = &[
    // Perfect match.
    Case {
        query: "ACGTACGTAC", subject: "ACGTACGTAC",
        score_cpg_mod: false, div_cpg_mod: false, complexity_adjust: false,
        score: 100, kimura: Some(0.0), cpg_sites: 2,
        pct_insert: 0.0, pct_delete: 0.0,
        well_characterized: 10, transitions: 0.0, transversions: 0,
    },
    // One transversion (A -> T at column 4).
    Case {
        query: "ACGTTCGTAC", subject: "ACGTACGTAC",
        score_cpg_mod: false, div_cpg_mod: false, complexity_adjust: false,
        score: 71, kimura: Some(10.846615), cpg_sites: 2,
        pct_insert: 0.0, pct_delete: 0.0,
        well_characterized: 10, transitions: 0.0, transversions: 1,
    },
    // Deletion: a 3 bp gap in the query.
    Case {
        query: "ACG---TAC", subject: "ACGTTGTAC",
        score_cpg_mod: false, div_cpg_mod: false, complexity_adjust: false,
        score: 25, kimura: Some(0.0), cpg_sites: 1,
        pct_insert: 0.0, pct_delete: 50.0,
        well_characterized: 6, transitions: 0.0, transversions: 0,
    },
    // Insertion: the same gap on the subject side.
    Case {
        query: "ACGTTGTAC", subject: "ACG---TAC",
        score_cpg_mod: false, div_cpg_mod: false, complexity_adjust: false,
        score: 25, kimura: Some(0.0), cpg_sites: 1,
        pct_insert: 50.0, pct_delete: 0.0,
        well_characterized: 6, transitions: 0.0, transversions: 0,
    },
    // CpG-rich, no modification: the G -> A transition counts in full.
    Case {
        query: "CATTACGCGTTACG", subject: "CGTTACGCGTTACG",
        score_cpg_mod: false, div_cpg_mod: false, complexity_adjust: false,
        score: 125, kimura: Some(7.707534), cpg_sites: 4,
        pct_insert: 0.0, pct_delete: 0.0,
        well_characterized: 14, transitions: 1.0, transversions: 0,
    },
    // Both CpG modifications: score rises and divergence collapses.
    Case {
        query: "CATTACGCGTTACG", subject: "CGTTACGCGTTACG",
        score_cpg_mod: true, div_cpg_mod: true, complexity_adjust: false,
        score: 144, kimura: Some(0.719437), cpg_sites: 4,
        pct_insert: 0.0, pct_delete: 0.0,
        well_characterized: 14, transitions: 0.1, transversions: 0,
    },
    // divCpGMod alone: divergence changes, score does not.
    Case {
        query: "CATTACGCGTTACG", subject: "CGTTACGCGTTACG",
        score_cpg_mod: false, div_cpg_mod: true, complexity_adjust: false,
        score: 125, kimura: Some(0.719437), cpg_sites: 4,
        pct_insert: 0.0, pct_delete: 0.0,
        well_characterized: 14, transitions: 0.1, transversions: 0,
    },
    // Complexity adjustment flattens a homopolymer to zero.
    Case {
        query: "AAAAAAAAAAAAAAAAAAAA", subject: "AAAAAAAAAAAAAAAAAAAA",
        score_cpg_mod: false, div_cpg_mod: false, complexity_adjust: true,
        score: 0, kimura: Some(0.0), cpg_sites: 0,
        pct_insert: 0.0, pct_delete: 0.0,
        well_characterized: 20, transitions: 0.0, transversions: 0,
    },
    // Complexity adjustment barely touches typical composition.
    Case {
        query: "ATATATATACGCGATATATAT", subject: "ATATATATACGCGATATATAT",
        score_cpg_mod: false, div_cpg_mod: false, complexity_adjust: true,
        score: 174, kimura: Some(0.0), cpg_sites: 2,
        pct_insert: 0.0, pct_delete: 0.0,
        well_characterized: 21, transitions: 0.0, transversions: 0,
    },
    // Two separate single-base deletions — two gap opens, no extensions.
    Case {
        query: "ACGT-ACG-TACG", subject: "ACGTAACGTTACG",
        score_cpg_mod: false, div_cpg_mod: false, complexity_adjust: false,
        score: 62, kimura: Some(0.0), cpg_sites: 3,
        pct_insert: 0.0, pct_delete: 18.181818,
        well_characterized: 11, transitions: 0.0, transversions: 0,
    },
    // Saturated: every other column is a transition, so K2P has no value.
    Case {
        query: "CGCGCGCGCGCG", subject: "CACACACACACA",
        score_cpg_mod: false, div_cpg_mod: true, complexity_adjust: false,
        score: 12, kimura: None, cpg_sites: 0,
        pct_insert: 0.0, pct_delete: 0.0,
        well_characterized: 12, transitions: 6.0, transversions: 0,
    },
];

#[test]
fn rescore_matches_repeatmasker() {
    let m = SubstMatrix::parse(M14P35G).unwrap();

    for (i, c) in CASES.iter().enumerate() {
        let params = RescoreParams {
            gap_open: -25,
            ins_gap_extend: -5,
            del_gap_extend: -5,
            score_cpg_mod: c.score_cpg_mod,
            div_cpg_mod: c.div_cpg_mod,
            complexity_adjust: c.complexity_adjust,
            ..RescoreParams::new(&m)
        };
        let r = rescore(c.query.as_bytes(), c.subject.as_bytes(), &params)
            .unwrap_or_else(|e| panic!("case {i} ({}/{}): {e}", c.query, c.subject));

        let label = format!(
            "case {i}: {} / {} (scoreCpG={} divCpG={} cadj={})",
            c.query, c.subject, c.score_cpg_mod, c.div_cpg_mod, c.complexity_adjust
        );

        assert_eq!(r.score, c.score, "{label}: score");
        assert_eq!(r.divergence.cpg_sites, c.cpg_sites, "{label}: CpG sites");
        assert_eq!(
            r.divergence.well_characterized, c.well_characterized,
            "{label}: well-characterised bases"
        );
        assert_eq!(
            r.divergence.transversions, c.transversions,
            "{label}: transversions"
        );
        assert!(
            (r.divergence.transitions - c.transitions).abs() < 1e-9,
            "{label}: transitions {} != {}",
            r.divergence.transitions,
            c.transitions
        );
        assert!(
            (r.pct_insert - c.pct_insert).abs() < 1e-6,
            "{label}: pctIns {} != {}",
            r.pct_insert,
            c.pct_insert
        );
        assert!(
            (r.pct_delete - c.pct_delete).abs() < 1e-6,
            "{label}: pctDel {} != {}",
            r.pct_delete,
            c.pct_delete
        );

        match (r.divergence.value, c.kimura) {
            (Some(got), Some(want)) => assert!(
                (got - want).abs() < 1e-5,
                "{label}: kimura {got} != {want}"
            ),
            (None, None) => {
                // The Perl reports its 100.00 initialiser here.
                assert_eq!(r.divergence.or_repeatmasker_default(100.0), 100.0, "{label}");
            }
            (got, want) => panic!("{label}: kimura {got:?} != {want:?}"),
        }
    }
}

/// The cumulative position scores must always end at the raw score, including
/// on the CpG-corrected path where earlier entries are retroactively patched.
#[test]
fn position_scores_end_at_the_raw_score() {
    let m = SubstMatrix::parse(M14P35G).unwrap();
    for c in CASES {
        let params = RescoreParams {
            gap_open: -25,
            ins_gap_extend: -5,
            del_gap_extend: -5,
            score_cpg_mod: c.score_cpg_mod,
            div_cpg_mod: c.div_cpg_mod,
            complexity_adjust: false,
            ..RescoreParams::new(&m)
        };
        let r = rescore(c.query.as_bytes(), c.subject.as_bytes(), &params).unwrap();
        assert_eq!(
            *r.position_scores.last().unwrap(),
            r.raw_score,
            "{} / {}",
            c.query,
            c.subject
        );
        assert_eq!(r.position_scores.len(), c.query.len());
    }
}
