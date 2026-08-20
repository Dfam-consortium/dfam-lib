//! Parity tests for [`aln_core::stats::kimura_divergence`] and
//! [`aln_core::stats::k2p_gap_divergence`] against RepeatMasker's
//! `SearchResult.pm::calcKimuraDivergence` / `calcK2PGapDivergence`.
//!
//! These are separate code paths from `rescoreAlignment`'s internal divergence
//! calculation — see `rescore_perl_parity.rs` for that one — and they differ
//! from it in observable ways, which is exactly why they are pinned separately.
//!
//! Expected values regenerated with:
//!
//! ```perl
//! use lib "/usr/local/RepeatMasker"; use SearchResult;
//! my $r = SearchResult->new(
//!     queryName => "q", subjName => "s",
//!     queryStart => 1, queryEnd => $ungapped_q,
//!     subjStart  => 1, subjEnd  => $ungapped_s,
//!     orientation => "", score => 0,
//!     queryString => $q, subjString => $s);
//! my ($kimura, $ti, $tv, $wcb, $cpg)          = $r->calcKimuraDivergence(divCpGMod => $d);
//! my ($k2pgap, $ti, $tv, $wcb, $cpg, $gaplen) = $r->calcK2PGapDivergence(divCpGMod => $d);
//! ```

use aln_core::stats::{Masking, k2p_gap_divergence, kimura_divergence};

struct Case {
    query: &'static str,
    subject: &'static str,
    div_cpg_mod: bool,

    /// `calcKimuraDivergence`: `None` where the Perl fell back to `100.00`.
    kimura: Option<f64>,
    /// `calcK2PGapDivergence` — note this stays defined under saturation.
    k2p_gap: Option<f64>,

    transitions: f64,
    transversions: u32,
    well_characterized: u32,
    /// CpG sites, which both routines count only when `div_cpg_mod` is on.
    cpg_sites: u32,
    gap_len: u32,
}

const CASES: &[Case] = &[
    Case {
        query: "ACGTACGTAC", subject: "ACGTACGTAC", div_cpg_mod: false,
        kimura: Some(0.0), k2p_gap: Some(0.0),
        transitions: 0.0, transversions: 0, well_characterized: 10,
        cpg_sites: 0, gap_len: 0,
    },
    // Same alignment, CpG counting switched on: two sites appear.
    Case {
        query: "ACGTACGTAC", subject: "ACGTACGTAC", div_cpg_mod: true,
        kimura: Some(0.0), k2p_gap: Some(0.0),
        transitions: 0.0, transversions: 0, well_characterized: 10,
        cpg_sites: 2, gap_len: 0,
    },
    Case {
        query: "CATTACGCGTTACG", subject: "CGTTACGCGTTACG", div_cpg_mod: false,
        kimura: Some(7.70753399), k2p_gap: Some(7.70753399),
        transitions: 1.0, transversions: 0, well_characterized: 14,
        cpg_sites: 0, gap_len: 0,
    },
    Case {
        query: "CATTACGCGTTACG", subject: "CGTTACGCGTTACG", div_cpg_mod: true,
        kimura: Some(0.71943687), k2p_gap: Some(0.71943687),
        transitions: 0.1, transversions: 0, well_characterized: 14,
        cpg_sites: 4, gap_len: 0,
    },
    // Query gaps: Kimura ignores them entirely, K2P-Gap charges for them.
    Case {
        query: "ACGT-ACG-TACG", subject: "ACGTAACGTTACG", div_cpg_mod: false,
        kimura: Some(0.0), k2p_gap: Some(6.02386456),
        transitions: 0.0, transversions: 0, well_characterized: 11,
        cpg_sites: 0, gap_len: 2,
    },
    Case {
        query: "ACGT-ACG-TACG", subject: "ACGTAACGTTACG", div_cpg_mod: true,
        kimura: Some(0.0), k2p_gap: Some(6.02386456),
        transitions: 0.0, transversions: 0, well_characterized: 11,
        cpg_sites: 3, gap_len: 2,
    },
    // Subject gaps: skipped outright by Kimura, counted by K2P-Gap.
    Case {
        query: "ACGTTGTAC", subject: "ACG---TAC", div_cpg_mod: false,
        kimura: Some(0.0), k2p_gap: Some(13.94647196),
        transitions: 0.0, transversions: 0, well_characterized: 6,
        cpg_sites: 0, gap_len: 3,
    },
    Case {
        query: "ACGTTGTAC", subject: "ACG---TAC", div_cpg_mod: true,
        kimura: Some(0.0), k2p_gap: Some(13.94647196),
        transitions: 0.0, transversions: 0, well_characterized: 6,
        cpg_sites: 1, gap_len: 3,
    },
    // Three transitions, two of which sit outside any CpG context, so the
    // discount changes nothing.
    Case {
        query: "TACGATACGATG", subject: "CACGGTACGGTG", div_cpg_mod: true,
        kimura: Some(34.65735903), k2p_gap: Some(34.65735903),
        transitions: 3.0, transversions: 0, well_characterized: 12,
        cpg_sites: 2, gap_len: 0,
    },
    Case {
        query: "TACGATACGATG", subject: "CACGGTACGGTG", div_cpg_mod: false,
        kimura: Some(34.65735903), k2p_gap: Some(34.65735903),
        transitions: 3.0, transversions: 0, well_characterized: 12,
        cpg_sites: 0, gap_len: 0,
    },
    // Saturated by transitions.  Kimura gives up; K2P-Gap substitutes a literal
    // 1 for the inner log term and returns -50.
    Case {
        query: "CGCGCGCGCGCG", subject: "CACACACACACA", div_cpg_mod: true,
        kimura: None, k2p_gap: Some(-50.0),
        transitions: 6.0, transversions: 0, well_characterized: 12,
        cpg_sites: 0, gap_len: 0,
    },
    // Saturated by transversions — the (1 - 2q) term goes negative.
    Case {
        query: "AAAACCCCGGGG", subject: "TTTTGGGGCCCC", div_cpg_mod: false,
        kimura: None, k2p_gap: Some(-50.0),
        transitions: 0.0, transversions: 12, well_characterized: 12,
        cpg_sites: 0, gap_len: 0,
    },
    Case {
        query: "ATGCATGCATGC", subject: "ATGCATGCATGC", div_cpg_mod: true,
        kimura: Some(0.0), k2p_gap: Some(0.0),
        transitions: 0.0, transversions: 0, well_characterized: 12,
        cpg_sites: 0, gap_len: 0,
    },
];

fn close(got: Option<f64>, want: Option<f64>, label: &str) {
    match (got, want) {
        (Some(g), Some(w)) => assert!((g - w).abs() < 1e-6, "{label}: {g} != {w}"),
        (None, None) => {}
        _ => panic!("{label}: {got:?} != {want:?}"),
    }
}

#[test]
fn kimura_divergence_matches_repeatmasker() {
    for (i, c) in CASES.iter().enumerate() {
        let d = kimura_divergence(c.query.as_bytes(), c.subject.as_bytes(), c.div_cpg_mod, Masking::Ignore)
            .unwrap();
        let label = format!("case {i}: {} / {} divCpG={}", c.query, c.subject, c.div_cpg_mod);

        close(d.value, c.kimura, &format!("{label}: kimura"));
        assert!(
            (d.transitions - c.transitions).abs() < 1e-9,
            "{label}: transitions {} != {}",
            d.transitions,
            c.transitions
        );
        assert_eq!(d.transversions, c.transversions, "{label}: transversions");
        assert_eq!(
            d.well_characterized, c.well_characterized,
            "{label}: well-characterised"
        );
        assert_eq!(d.cpg_sites, c.cpg_sites, "{label}: CpG sites");
    }
}

#[test]
fn k2p_gap_divergence_matches_repeatmasker() {
    for (i, c) in CASES.iter().enumerate() {
        let d = k2p_gap_divergence(c.query.as_bytes(), c.subject.as_bytes(), c.div_cpg_mod, Masking::Ignore)
            .unwrap();
        let label = format!("case {i}: {} / {} divCpG={}", c.query, c.subject, c.div_cpg_mod);

        close(d.value, c.k2p_gap, &format!("{label}: k2p-gap"));
        assert_eq!(d.gap_len, c.gap_len, "{label}: gap length");
        assert_eq!(
            d.well_characterized, c.well_characterized,
            "{label}: well-characterised"
        );
        assert_eq!(d.cpg_sites, c.cpg_sites, "{label}: CpG sites");
    }
}

/// The two routines diverge precisely where gaps are involved: `calcKimura`
/// ignores gap columns, `calcK2PGap` penalises them.  Without gaps they agree.
#[test]
fn the_two_routines_agree_only_when_there_are_no_gaps() {
    for c in CASES {
        let k = kimura_divergence(c.query.as_bytes(), c.subject.as_bytes(), c.div_cpg_mod, Masking::Ignore)
            .unwrap();
        let g = k2p_gap_divergence(c.query.as_bytes(), c.subject.as_bytes(), c.div_cpg_mod, Masking::Ignore)
            .unwrap();
        if g.gap_len == 0 {
            if let (Some(kv), Some(gv)) = (k.value, g.value) {
                assert!(
                    (kv - gv).abs() < 1e-6,
                    "{} / {}: ungapped but {kv} != {gv}",
                    c.query,
                    c.subject
                );
            }
        } else {
            assert_ne!(k.value, g.value, "{} / {}", c.query, c.subject);
        }
    }
}
