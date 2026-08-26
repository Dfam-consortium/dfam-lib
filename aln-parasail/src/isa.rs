//! Runtime SIMD selection.
//!
//! parasail ships its own cpuid/xgetbv dispatcher, but it is not vendored here:
//! `is_x86_feature_detected!` does the same job, is already correct about
//! OS-level AVX state (the part hand-rolled cpuid code usually gets wrong), and
//! keeps three more C files out of the build.
//!
//! Each kernel is compiled into its own static library with its own `-m` flag
//! (see `build.rs`), so calling an AVX2 kernel on a host without AVX2 would be
//! an illegal instruction — this module is what prevents that.

use crate::ffi;

/// SIMD instruction set a kernel was compiled for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Isa {
    Avx2,
    Sse41,
    Sse2,
}

impl Isa {
    /// The widest instruction set this CPU supports.
    ///
    /// x86-64 guarantees SSE2, so this never fails on a supported target.
    pub fn detect() -> Isa {
        #[cfg(any(target_arch = "x86_64", target_arch = "x86"))]
        {
            if std::arch::is_x86_feature_detected!("avx2") {
                return Isa::Avx2;
            }
            if std::arch::is_x86_feature_detected!("sse4.1") {
                return Isa::Sse41;
            }
        }
        Isa::Sse2
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Isa::Avx2 => "avx2",
            Isa::Sse41 => "sse41",
            Isa::Sse2 => "sse2",
        }
    }

    /// SIMD register width in bits.  Profile layout depends on this, so a
    /// profile built for one ISA cannot be handed to a kernel of another.
    pub fn register_bits(self) -> usize {
        match self {
            Isa::Avx2 => 256,
            Isa::Sse41 | Isa::Sse2 => 128,
        }
    }

    /// The matching profile constructor.
    ///
    /// The `_sat` constructors fill the 8-, 16- and 32-bit layouts at once, so
    /// one profile serves the whole saturation-fallback chain.
    pub fn profile_create(self) -> ffi::ProfileCreateFn {
        match self {
            Isa::Avx2 => ffi::parasail_profile_create_avx_256_sat,
            Isa::Sse41 | Isa::Sse2 => ffi::parasail_profile_create_sse_128_sat,
        }
    }
}

/// Lane width to attempt.  Narrower lanes are faster but overflow sooner;
/// parasail reports overflow via `parasail_result_is_saturated`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Width {
    W8,
    W16,
    W32,
}

/// Widths in the order the fallback chain tries them.
pub const FALLBACK_CHAIN: [Width; 3] = [Width::W8, Width::W16, Width::W32];

/// Which recurrence to run.
///
/// The semi-global variants are named from **parasail's** point of view, where
/// `q` is `s1` and `d` is `s2`.  Because this crate profiles the *subject* as
/// parasail's `s1`, the mapping to [`AlignMode`](aln_engine::AlignMode) is
/// inverted — see `ParasailAligner::kernel`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kernel {
    /// Smith-Waterman, local.
    Sw,
    /// Needleman-Wunsch, global.
    Nw,
    /// Semi-global, all four ends free.
    Sg,
    /// Semi-global, both ends of `s1` free.
    SgQx,
    /// Semi-global, both ends of `s2` free.
    SgDx,
}

macro_rules! dispatch {
    ($kernel:expr, $isa:expr, $width:expr, $v:ident, $( $k:ident => $prefix:ident ),+ $(,)?) => {
        match $kernel {
            $(
                Kernel::$k => match ($isa, $width) {
                    (Isa::Sse2,  Width::W8)  => paste_fn!($v, $prefix, sse2_128_8),
                    (Isa::Sse2,  Width::W16) => paste_fn!($v, $prefix, sse2_128_16),
                    (Isa::Sse2,  Width::W32) => paste_fn!($v, $prefix, sse2_128_32),
                    (Isa::Sse41, Width::W8)  => paste_fn!($v, $prefix, sse41_128_8),
                    (Isa::Sse41, Width::W16) => paste_fn!($v, $prefix, sse41_128_16),
                    (Isa::Sse41, Width::W32) => paste_fn!($v, $prefix, sse41_128_32),
                    (Isa::Avx2,  Width::W8)  => paste_fn!($v, $prefix, avx2_256_8),
                    (Isa::Avx2,  Width::W16) => paste_fn!($v, $prefix, avx2_256_16),
                    (Isa::Avx2,  Width::W32) => paste_fn!($v, $prefix, avx2_256_32),
                },
            )+
        }
    };
}

/// Concatenating identifiers needs either `paste` or explicit spelling.  The
/// table is small and fixed, so it is spelled out rather than taking a
/// proc-macro dependency.
macro_rules! paste_fn {
    ($v:ident, sw, $suffix:ident) => { paste_sw!($v, $suffix) };
    ($v:ident, nw, $suffix:ident) => { paste_nw!($v, $suffix) };
    ($v:ident, sg, $suffix:ident) => { paste_sg!($v, $suffix) };
    ($v:ident, sg_qx, $suffix:ident) => { paste_sg_qx!($v, $suffix) };
    ($v:ident, sg_dx, $suffix:ident) => { paste_sg_dx!($v, $suffix) };
}

macro_rules! paste_sw {
    (trace, sse2_128_8) => { ffi::parasail_sw_trace_striped_profile_sse2_128_8 };
    (score, sse2_128_8) => { ffi::parasail_sw_striped_profile_sse2_128_8 };
    (trace, sse2_128_16) => { ffi::parasail_sw_trace_striped_profile_sse2_128_16 };
    (score, sse2_128_16) => { ffi::parasail_sw_striped_profile_sse2_128_16 };
    (trace, sse2_128_32) => { ffi::parasail_sw_trace_striped_profile_sse2_128_32 };
    (score, sse2_128_32) => { ffi::parasail_sw_striped_profile_sse2_128_32 };
    (trace, sse41_128_8) => { ffi::parasail_sw_trace_striped_profile_sse41_128_8 };
    (score, sse41_128_8) => { ffi::parasail_sw_striped_profile_sse41_128_8 };
    (trace, sse41_128_16) => { ffi::parasail_sw_trace_striped_profile_sse41_128_16 };
    (score, sse41_128_16) => { ffi::parasail_sw_striped_profile_sse41_128_16 };
    (trace, sse41_128_32) => { ffi::parasail_sw_trace_striped_profile_sse41_128_32 };
    (score, sse41_128_32) => { ffi::parasail_sw_striped_profile_sse41_128_32 };
    (trace, avx2_256_8) => { ffi::parasail_sw_trace_striped_profile_avx2_256_8 };
    (score, avx2_256_8) => { ffi::parasail_sw_striped_profile_avx2_256_8 };
    (trace, avx2_256_16) => { ffi::parasail_sw_trace_striped_profile_avx2_256_16 };
    (score, avx2_256_16) => { ffi::parasail_sw_striped_profile_avx2_256_16 };
    (trace, avx2_256_32) => { ffi::parasail_sw_trace_striped_profile_avx2_256_32 };
    (score, avx2_256_32) => { ffi::parasail_sw_striped_profile_avx2_256_32 };
}

macro_rules! paste_nw {
    (trace, sse2_128_8) => { ffi::parasail_nw_trace_striped_profile_sse2_128_8 };
    (score, sse2_128_8) => { ffi::parasail_nw_striped_profile_sse2_128_8 };
    (trace, sse2_128_16) => { ffi::parasail_nw_trace_striped_profile_sse2_128_16 };
    (score, sse2_128_16) => { ffi::parasail_nw_striped_profile_sse2_128_16 };
    (trace, sse2_128_32) => { ffi::parasail_nw_trace_striped_profile_sse2_128_32 };
    (score, sse2_128_32) => { ffi::parasail_nw_striped_profile_sse2_128_32 };
    (trace, sse41_128_8) => { ffi::parasail_nw_trace_striped_profile_sse41_128_8 };
    (score, sse41_128_8) => { ffi::parasail_nw_striped_profile_sse41_128_8 };
    (trace, sse41_128_16) => { ffi::parasail_nw_trace_striped_profile_sse41_128_16 };
    (score, sse41_128_16) => { ffi::parasail_nw_striped_profile_sse41_128_16 };
    (trace, sse41_128_32) => { ffi::parasail_nw_trace_striped_profile_sse41_128_32 };
    (score, sse41_128_32) => { ffi::parasail_nw_striped_profile_sse41_128_32 };
    (trace, avx2_256_8) => { ffi::parasail_nw_trace_striped_profile_avx2_256_8 };
    (score, avx2_256_8) => { ffi::parasail_nw_striped_profile_avx2_256_8 };
    (trace, avx2_256_16) => { ffi::parasail_nw_trace_striped_profile_avx2_256_16 };
    (score, avx2_256_16) => { ffi::parasail_nw_striped_profile_avx2_256_16 };
    (trace, avx2_256_32) => { ffi::parasail_nw_trace_striped_profile_avx2_256_32 };
    (score, avx2_256_32) => { ffi::parasail_nw_striped_profile_avx2_256_32 };
}

macro_rules! paste_sg {
    (trace, sse2_128_8) => { ffi::parasail_sg_trace_striped_profile_sse2_128_8 };
    (score, sse2_128_8) => { ffi::parasail_sg_striped_profile_sse2_128_8 };
    (trace, sse2_128_16) => { ffi::parasail_sg_trace_striped_profile_sse2_128_16 };
    (score, sse2_128_16) => { ffi::parasail_sg_striped_profile_sse2_128_16 };
    (trace, sse2_128_32) => { ffi::parasail_sg_trace_striped_profile_sse2_128_32 };
    (score, sse2_128_32) => { ffi::parasail_sg_striped_profile_sse2_128_32 };
    (trace, sse41_128_8) => { ffi::parasail_sg_trace_striped_profile_sse41_128_8 };
    (score, sse41_128_8) => { ffi::parasail_sg_striped_profile_sse41_128_8 };
    (trace, sse41_128_16) => { ffi::parasail_sg_trace_striped_profile_sse41_128_16 };
    (score, sse41_128_16) => { ffi::parasail_sg_striped_profile_sse41_128_16 };
    (trace, sse41_128_32) => { ffi::parasail_sg_trace_striped_profile_sse41_128_32 };
    (score, sse41_128_32) => { ffi::parasail_sg_striped_profile_sse41_128_32 };
    (trace, avx2_256_8) => { ffi::parasail_sg_trace_striped_profile_avx2_256_8 };
    (score, avx2_256_8) => { ffi::parasail_sg_striped_profile_avx2_256_8 };
    (trace, avx2_256_16) => { ffi::parasail_sg_trace_striped_profile_avx2_256_16 };
    (score, avx2_256_16) => { ffi::parasail_sg_striped_profile_avx2_256_16 };
    (trace, avx2_256_32) => { ffi::parasail_sg_trace_striped_profile_avx2_256_32 };
    (score, avx2_256_32) => { ffi::parasail_sg_striped_profile_avx2_256_32 };
}

macro_rules! paste_sg_qx {
    (trace, sse2_128_8) => { ffi::parasail_sg_qx_trace_striped_profile_sse2_128_8 };
    (score, sse2_128_8) => { ffi::parasail_sg_qx_striped_profile_sse2_128_8 };
    (trace, sse2_128_16) => { ffi::parasail_sg_qx_trace_striped_profile_sse2_128_16 };
    (score, sse2_128_16) => { ffi::parasail_sg_qx_striped_profile_sse2_128_16 };
    (trace, sse2_128_32) => { ffi::parasail_sg_qx_trace_striped_profile_sse2_128_32 };
    (score, sse2_128_32) => { ffi::parasail_sg_qx_striped_profile_sse2_128_32 };
    (trace, sse41_128_8) => { ffi::parasail_sg_qx_trace_striped_profile_sse41_128_8 };
    (score, sse41_128_8) => { ffi::parasail_sg_qx_striped_profile_sse41_128_8 };
    (trace, sse41_128_16) => { ffi::parasail_sg_qx_trace_striped_profile_sse41_128_16 };
    (score, sse41_128_16) => { ffi::parasail_sg_qx_striped_profile_sse41_128_16 };
    (trace, sse41_128_32) => { ffi::parasail_sg_qx_trace_striped_profile_sse41_128_32 };
    (score, sse41_128_32) => { ffi::parasail_sg_qx_striped_profile_sse41_128_32 };
    (trace, avx2_256_8) => { ffi::parasail_sg_qx_trace_striped_profile_avx2_256_8 };
    (score, avx2_256_8) => { ffi::parasail_sg_qx_striped_profile_avx2_256_8 };
    (trace, avx2_256_16) => { ffi::parasail_sg_qx_trace_striped_profile_avx2_256_16 };
    (score, avx2_256_16) => { ffi::parasail_sg_qx_striped_profile_avx2_256_16 };
    (trace, avx2_256_32) => { ffi::parasail_sg_qx_trace_striped_profile_avx2_256_32 };
    (score, avx2_256_32) => { ffi::parasail_sg_qx_striped_profile_avx2_256_32 };
}

macro_rules! paste_sg_dx {
    (trace, sse2_128_8) => { ffi::parasail_sg_dx_trace_striped_profile_sse2_128_8 };
    (score, sse2_128_8) => { ffi::parasail_sg_dx_striped_profile_sse2_128_8 };
    (trace, sse2_128_16) => { ffi::parasail_sg_dx_trace_striped_profile_sse2_128_16 };
    (score, sse2_128_16) => { ffi::parasail_sg_dx_striped_profile_sse2_128_16 };
    (trace, sse2_128_32) => { ffi::parasail_sg_dx_trace_striped_profile_sse2_128_32 };
    (score, sse2_128_32) => { ffi::parasail_sg_dx_striped_profile_sse2_128_32 };
    (trace, sse41_128_8) => { ffi::parasail_sg_dx_trace_striped_profile_sse41_128_8 };
    (score, sse41_128_8) => { ffi::parasail_sg_dx_striped_profile_sse41_128_8 };
    (trace, sse41_128_16) => { ffi::parasail_sg_dx_trace_striped_profile_sse41_128_16 };
    (score, sse41_128_16) => { ffi::parasail_sg_dx_striped_profile_sse41_128_16 };
    (trace, sse41_128_32) => { ffi::parasail_sg_dx_trace_striped_profile_sse41_128_32 };
    (score, sse41_128_32) => { ffi::parasail_sg_dx_striped_profile_sse41_128_32 };
    (trace, avx2_256_8) => { ffi::parasail_sg_dx_trace_striped_profile_avx2_256_8 };
    (score, avx2_256_8) => { ffi::parasail_sg_dx_striped_profile_avx2_256_8 };
    (trace, avx2_256_16) => { ffi::parasail_sg_dx_trace_striped_profile_avx2_256_16 };
    (score, avx2_256_16) => { ffi::parasail_sg_dx_striped_profile_avx2_256_16 };
    (trace, avx2_256_32) => { ffi::parasail_sg_dx_trace_striped_profile_avx2_256_32 };
    (score, avx2_256_32) => { ffi::parasail_sg_dx_striped_profile_avx2_256_32 };
}

/// Resolve `(kernel, ISA, lane width)` to a traceback entry point.
pub fn resolve(kernel: Kernel, isa: Isa, width: Width) -> ffi::TraceProfileFn {
    dispatch!(kernel, isa, width, trace,
        Sw => sw,
        Nw => nw,
        Sg => sg,
        SgQx => sg_qx,
        SgDx => sg_dx,
    )
}

/// Resolve to the score-only entry point for the same kernel.
///
/// These skip the O(mn) traceback matrix, so `parasail_result_get_cigar` is not
/// valid on their results — only the score and end positions are populated.
pub fn resolve_score_only(kernel: Kernel, isa: Isa, width: Width) -> ffi::TraceProfileFn {
    dispatch!(kernel, isa, width, score,
        Sw => sw,
        Nw => nw,
        Sg => sg,
        SgQx => sg_qx,
        SgDx => sg_dx,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detect_returns_something_this_cpu_supports() {
        let isa = Isa::detect();
        match isa {
            Isa::Avx2 => assert!(std::arch::is_x86_feature_detected!("avx2")),
            Isa::Sse41 => assert!(std::arch::is_x86_feature_detected!("sse4.1")),
            Isa::Sse2 => {}
        }
    }

    #[test]
    fn every_combination_resolves_to_a_distinct_symbol() {
        let mut seen = std::collections::HashSet::new();
        for kernel in [Kernel::Sw, Kernel::Nw, Kernel::Sg, Kernel::SgQx, Kernel::SgDx] {
            for isa in [Isa::Sse2, Isa::Sse41, Isa::Avx2] {
                for width in FALLBACK_CHAIN {
                    let f = resolve(kernel, isa, width) as usize;
                    assert!(
                        seen.insert(f),
                        "duplicate symbol for {kernel:?}/{isa:?}/{width:?} — \
                         a macro arm is wired to the wrong kernel"
                    );
                }
            }
        }
        assert_eq!(seen.len(), 45);
    }

    #[test]
    fn register_width_tracks_the_isa() {
        assert_eq!(Isa::Avx2.register_bits(), 256);
        assert_eq!(Isa::Sse2.register_bits(), 128);
    }
}
