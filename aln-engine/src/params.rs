//! Alignment parameters and backend capability reporting.

use aln_core::SubstMatrix;

/// Which dynamic-programming recurrence to run.
///
/// The names map onto parasail's function families (`sw_`, `nw_`, `sg_`) and
/// onto the classical algorithms.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AlignMode {
    /// Smith-Waterman: best local alignment, free ends on both sequences.
    /// GIRI's `SWAligner` and `acons` use this.
    #[default]
    Local,
    /// Needleman-Wunsch: both sequences aligned end to end.
    Global,
    /// Semi-global: end gaps are free on the flagged sides.
    ///
    /// `free_query_ends = false, free_subject_ends = true` is the usual
    /// "fit a consensus into a genomic region" shape.
    SemiGlobal {
        free_query_ends: bool,
        free_subject_ends: bool,
    },
}

/// Everything a [`PairwiseAligner`](crate::PairwiseAligner) needs to score.
///
/// Gap penalties are **positive magnitudes** here — the aligner subtracts them.
/// This differs from [`SubstMatrix::gap_open`], which reports the signed value
/// as written in a GIRI `GAP` line; [`AlignParams::from_matrix`] flips the sign.
#[derive(Debug, Clone)]
pub struct AlignParams {
    pub mode: AlignMode,

    /// Cost of opening a gap, including its first position.
    pub gap_open: u32,
    /// Cost of each position beyond the first.
    pub gap_extend: u32,

    /// Discard alignments scoring below this.  GIRI keeps anything `> 0`.
    pub min_score: i32,

    /// Ask the backend for a traceback.  Turning this off lets score-only
    /// kernels run, which are markedly faster in parasail.
    pub traceback: bool,

    /// Band half-width for banded backends; `None` means unbanded.
    pub bandwidth: Option<u32>,
}

impl Default for AlignParams {
    fn default() -> Self {
        AlignParams {
            mode: AlignMode::Local,
            gap_open: 25,
            gap_extend: 5,
            min_score: 1,
            traceback: true,
            bandwidth: None,
        }
    }
}

impl AlignParams {
    /// Take gap penalties from a matrix's `GAP` line when it carries one,
    /// converting the signed GIRI values to positive magnitudes.
    pub fn from_matrix(matrix: &SubstMatrix) -> Self {
        let mut p = AlignParams::default();
        if let Some(go) = matrix.gap_open() {
            p.gap_open = go.unsigned_abs();
        }
        if let Some(ge) = matrix.gap_extend() {
            p.gap_extend = ge.unsigned_abs();
        }
        p
    }

    /// Reject combinations no backend can honour.
    pub fn validate(&self) -> crate::Result<()> {
        if self.gap_extend > self.gap_open {
            return Err(crate::EngineError::Params(format!(
                "gap_extend ({}) exceeds gap_open ({}); affine gaps would get \
                 cheaper as they lengthen",
                self.gap_extend, self.gap_open
            )));
        }
        Ok(())
    }
}

/// What a backend can actually do.
///
/// Checked by the drivers before dispatch so a mismatch surfaces as a clear
/// error rather than a silently wrong alignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignerCaps {
    /// Backend name, for diagnostics.
    pub name: &'static str,
    /// Modes the backend implements.
    pub modes: &'static [AlignMode],
    /// Whether it can return an edit script at all.  Several of parasail's
    /// fastest striped kernels are score-only.
    pub traceback: bool,
    /// Whether it honours [`AlignParams::bandwidth`].
    pub banded: bool,
    /// SIMD instruction set in use, if any — `"sse2"`, `"avx2"`, `"scalar"`.
    pub simd: &'static str,
}

impl AlignerCaps {
    pub fn supports(&self, params: &AlignParams) -> Option<String> {
        if params.traceback && !self.traceback {
            return Some(format!("{} is score-only; no traceback available", self.name));
        }
        if params.bandwidth.is_some() && !self.banded {
            return Some(format!("{} does not implement banding", self.name));
        }
        if !self.modes.contains(&params.mode) {
            return Some(format!("{} does not implement {:?}", self.name, params.mode));
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const CAPS: AlignerCaps = AlignerCaps {
        name: "test",
        modes: &[AlignMode::Local],
        traceback: false,
        banded: false,
        simd: "scalar",
    };

    #[test]
    fn caps_reject_traceback_from_a_score_only_backend() {
        let p = AlignParams { traceback: true, ..Default::default() };
        assert!(CAPS.supports(&p).is_some());
        let p = AlignParams { traceback: false, ..Default::default() };
        assert!(CAPS.supports(&p).is_none());
    }

    #[test]
    fn caps_reject_an_unimplemented_mode() {
        let p = AlignParams { traceback: false, mode: AlignMode::Global, ..Default::default() };
        assert!(CAPS.supports(&p).is_some());
    }

    #[test]
    fn validate_rejects_gap_extend_above_gap_open() {
        let p = AlignParams { gap_open: 5, gap_extend: 10, ..Default::default() };
        assert!(p.validate().is_err());
    }

    #[test]
    fn gap_penalties_from_a_giri_matrix_become_positive_magnitudes() {
        let m = SubstMatrix::parse("GAP -25 -5\n  A   C\n  1  -1\n -1   1\n").unwrap();
        let p = AlignParams::from_matrix(&m);
        assert_eq!(p.gap_open, 25);
        assert_eq!(p.gap_extend, 5);
    }
}
