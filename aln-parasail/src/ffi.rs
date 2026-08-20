//! Raw declarations for the vendored parasail subset.
//!
//! Written by hand rather than generated with bindgen: the surface is about
//! twenty functions plus two structs, and hand-writing keeps `bindgen` (and
//! therefore libclang) off the build-dependency list.
//!
//! Only the striped **traceback** kernels are vendored.  parasail's own
//! `_sat` wrappers and cpuid dispatcher are not — saturation fallback and ISA
//! selection happen in Rust (see [`crate::isa`] and
//! [`ParasailAligner::align_prepared`](crate::ParasailAligner::align_prepared)).

#![allow(non_camel_case_types)]

use std::os::raw::{c_char, c_int, c_void};

/// `parasail_matrix_t`.  Field order must match `parasail.h` exactly.
#[repr(C)]
pub struct parasail_matrix_t {
    pub name: *const c_char,
    pub matrix: *const c_int,
    pub mapper: *const c_int,
    /// Side length of the square matrix — note this is `alphabet.len() + 1`;
    /// parasail appends a `'*'` catch-all row and column.
    pub size: c_int,
    pub max: c_int,
    pub min: c_int,
    pub user_matrix: *mut c_int,
    pub type_: c_int,
    pub length: c_int,
    pub alphabet: *const c_char,
    pub query: *const c_char,
}

#[repr(C)]
pub struct parasail_profile_data {
    pub score: *mut c_void,
    pub matches: *mut c_void,
    pub similar: *mut c_void,
}

/// `parasail_profile_t`.
///
/// Note `s1` is **borrowed**, not copied — `parasail_profile_new` stores the
/// caller's pointer.  Whatever owns those bytes must outlive the profile.
#[repr(C)]
pub struct parasail_profile_t {
    pub s1: *const c_char,
    pub s1_len: c_int,
    pub matrix: *const parasail_matrix_t,
    pub profile8: parasail_profile_data,
    pub profile16: parasail_profile_data,
    pub profile32: parasail_profile_data,
    pub profile64: parasail_profile_data,
    pub free: Option<unsafe extern "C" fn(*mut c_void)>,
    pub stop: c_int,
}

/// Opaque — only reached through the accessor functions below.
#[repr(C)]
pub struct parasail_result_t {
    _private: [u8; 0],
}

#[repr(C)]
pub struct parasail_cigar_t {
    /// Packed `(length << 4) | op_code`; decode with [`parasail_cigar_decode_op`]
    /// and [`parasail_cigar_decode_len`].
    pub seq: *mut u32,
    pub len: c_int,
    /// 0-based start offset in `s1`.
    pub beg_query: c_int,
    /// 0-based start offset in `s2`.
    pub beg_ref: c_int,
}

/// A striped kernel taking a prepared profile — traceback or score-only;
/// both families share this signature.
///
/// `(profile, s2, s2_len, gap_open, gap_extend) -> result`.  Gap penalties are
/// positive magnitudes.
pub type TraceProfileFn = unsafe extern "C" fn(
    *const parasail_profile_t,
    *const c_char,
    c_int,
    c_int,
    c_int,
) -> *mut parasail_result_t;

/// A profile constructor for one SIMD register width.
pub type ProfileCreateFn =
    unsafe extern "C" fn(*const c_char, c_int, *const parasail_matrix_t) -> *mut parasail_profile_t;

extern "C" {
    // ── Matrices ──────────────────────────────────────────────────────────
    pub fn parasail_matrix_create(
        alphabet: *const c_char,
        match_: c_int,
        mismatch: c_int,
    ) -> *mut parasail_matrix_t;
    pub fn parasail_matrix_set_value(
        matrix: *mut parasail_matrix_t,
        row: c_int,
        col: c_int,
        value: c_int,
    );
    pub fn parasail_matrix_free(matrix: *mut parasail_matrix_t);

    // ── Profiles ──────────────────────────────────────────────────────────
    // The `_sat` constructors populate the 8-, 16- and 32-bit layouts in one
    // object, which is what lets a single profile serve the whole saturation
    // fallback chain.
    pub fn parasail_profile_create_sse_128_sat(
        s1: *const c_char,
        s1_len: c_int,
        matrix: *const parasail_matrix_t,
    ) -> *mut parasail_profile_t;
    pub fn parasail_profile_create_avx_256_sat(
        s1: *const c_char,
        s1_len: c_int,
        matrix: *const parasail_matrix_t,
    ) -> *mut parasail_profile_t;
    pub fn parasail_profile_free(profile: *mut parasail_profile_t);

    // ── Results ───────────────────────────────────────────────────────────
    pub fn parasail_result_get_score(result: *const parasail_result_t) -> c_int;
    /// Last aligned offset in `s1` (inclusive).
    pub fn parasail_result_get_end_query(result: *const parasail_result_t) -> c_int;
    /// Last aligned offset in `s2` (inclusive).
    pub fn parasail_result_get_end_ref(result: *const parasail_result_t) -> c_int;
    /// Non-zero when the lane width overflowed and the result is unusable.
    pub fn parasail_result_is_saturated(result: *const parasail_result_t) -> c_int;
    pub fn parasail_result_free(result: *mut parasail_result_t);

    // ── CIGAR ─────────────────────────────────────────────────────────────
    pub fn parasail_result_get_cigar(
        result: *mut parasail_result_t,
        seq_a: *const c_char,
        len_a: c_int,
        seq_b: *const c_char,
        len_b: c_int,
        matrix: *const parasail_matrix_t,
    ) -> *mut parasail_cigar_t;
    pub fn parasail_cigar_decode_op(cigar_int: u32) -> c_char;
    pub fn parasail_cigar_decode_len(cigar_int: u32) -> u32;
    pub fn parasail_cigar_free(cigar: *mut parasail_cigar_t);
}

/// Declare one kernel.
///
/// Naming is `parasail_<alg>_trace_striped_profile_<isa>_<width>`.  The `sg`
/// source additionally emits every semi-global variant via `sg_helper.h`, so
/// `sg_qx` and `sg_dx` come from the same object file as `sg`.
macro_rules! decl_kernel {
    ($name:ident) => {
        extern "C" {
            pub fn $name(
                profile: *const parasail_profile_t,
                s2: *const c_char,
                s2_len: c_int,
                open: c_int,
                gap: c_int,
            ) -> *mut parasail_result_t;
        }
    };
}

macro_rules! decl_family {
    ($($name:ident),+ $(,)?) => {
        $( decl_kernel!($name); )+
    };
}

decl_family!(
    // Smith-Waterman — local.
    parasail_sw_trace_striped_profile_sse2_128_8,
    parasail_sw_trace_striped_profile_sse2_128_16,
    parasail_sw_trace_striped_profile_sse2_128_32,
    parasail_sw_trace_striped_profile_sse41_128_8,
    parasail_sw_trace_striped_profile_sse41_128_16,
    parasail_sw_trace_striped_profile_sse41_128_32,
    parasail_sw_trace_striped_profile_avx2_256_8,
    parasail_sw_trace_striped_profile_avx2_256_16,
    parasail_sw_trace_striped_profile_avx2_256_32,
    // Needleman-Wunsch — global.
    parasail_nw_trace_striped_profile_sse2_128_8,
    parasail_nw_trace_striped_profile_sse2_128_16,
    parasail_nw_trace_striped_profile_sse2_128_32,
    parasail_nw_trace_striped_profile_sse41_128_8,
    parasail_nw_trace_striped_profile_sse41_128_16,
    parasail_nw_trace_striped_profile_sse41_128_32,
    parasail_nw_trace_striped_profile_avx2_256_8,
    parasail_nw_trace_striped_profile_avx2_256_16,
    parasail_nw_trace_striped_profile_avx2_256_32,
    // Semi-global, all four ends free.
    parasail_sg_trace_striped_profile_sse2_128_8,
    parasail_sg_trace_striped_profile_sse2_128_16,
    parasail_sg_trace_striped_profile_sse2_128_32,
    parasail_sg_trace_striped_profile_sse41_128_8,
    parasail_sg_trace_striped_profile_sse41_128_16,
    parasail_sg_trace_striped_profile_sse41_128_32,
    parasail_sg_trace_striped_profile_avx2_256_8,
    parasail_sg_trace_striped_profile_avx2_256_16,
    parasail_sg_trace_striped_profile_avx2_256_32,
    // Semi-global, both ends of s1 free (`qx`).
    parasail_sg_qx_trace_striped_profile_sse2_128_8,
    parasail_sg_qx_trace_striped_profile_sse2_128_16,
    parasail_sg_qx_trace_striped_profile_sse2_128_32,
    parasail_sg_qx_trace_striped_profile_sse41_128_8,
    parasail_sg_qx_trace_striped_profile_sse41_128_16,
    parasail_sg_qx_trace_striped_profile_sse41_128_32,
    parasail_sg_qx_trace_striped_profile_avx2_256_8,
    parasail_sg_qx_trace_striped_profile_avx2_256_16,
    parasail_sg_qx_trace_striped_profile_avx2_256_32,
    // Semi-global, both ends of s2 free (`dx`).
    parasail_sg_dx_trace_striped_profile_sse2_128_8,
    parasail_sg_dx_trace_striped_profile_sse2_128_16,
    parasail_sg_dx_trace_striped_profile_sse2_128_32,
    parasail_sg_dx_trace_striped_profile_sse41_128_8,
    parasail_sg_dx_trace_striped_profile_sse41_128_16,
    parasail_sg_dx_trace_striped_profile_sse41_128_32,
    parasail_sg_dx_trace_striped_profile_avx2_256_8,
    parasail_sg_dx_trace_striped_profile_avx2_256_16,
    parasail_sg_dx_trace_striped_profile_avx2_256_32,
);

// The score-only counterparts.  Same signature; they skip the O(mn) traceback
// matrix and keep only O(m) column state, so they are what a score prepass
// should call.  `parasail_result_get_cigar` is NOT valid on their results —
// only the score and end positions are populated.
decl_family!(
    // Smith-Waterman — local.
    parasail_sw_striped_profile_sse2_128_8,
    parasail_sw_striped_profile_sse2_128_16,
    parasail_sw_striped_profile_sse2_128_32,
    parasail_sw_striped_profile_sse41_128_8,
    parasail_sw_striped_profile_sse41_128_16,
    parasail_sw_striped_profile_sse41_128_32,
    parasail_sw_striped_profile_avx2_256_8,
    parasail_sw_striped_profile_avx2_256_16,
    parasail_sw_striped_profile_avx2_256_32,
    // Needleman-Wunsch — global.
    parasail_nw_striped_profile_sse2_128_8,
    parasail_nw_striped_profile_sse2_128_16,
    parasail_nw_striped_profile_sse2_128_32,
    parasail_nw_striped_profile_sse41_128_8,
    parasail_nw_striped_profile_sse41_128_16,
    parasail_nw_striped_profile_sse41_128_32,
    parasail_nw_striped_profile_avx2_256_8,
    parasail_nw_striped_profile_avx2_256_16,
    parasail_nw_striped_profile_avx2_256_32,
    // semi-global, both ends free.
    parasail_sg_striped_profile_sse2_128_8,
    parasail_sg_striped_profile_sse2_128_16,
    parasail_sg_striped_profile_sse2_128_32,
    parasail_sg_striped_profile_sse41_128_8,
    parasail_sg_striped_profile_sse41_128_16,
    parasail_sg_striped_profile_sse41_128_32,
    parasail_sg_striped_profile_avx2_256_8,
    parasail_sg_striped_profile_avx2_256_16,
    parasail_sg_striped_profile_avx2_256_32,
    // semi-global, both ends of s1 free.
    parasail_sg_qx_striped_profile_sse2_128_8,
    parasail_sg_qx_striped_profile_sse2_128_16,
    parasail_sg_qx_striped_profile_sse2_128_32,
    parasail_sg_qx_striped_profile_sse41_128_8,
    parasail_sg_qx_striped_profile_sse41_128_16,
    parasail_sg_qx_striped_profile_sse41_128_32,
    parasail_sg_qx_striped_profile_avx2_256_8,
    parasail_sg_qx_striped_profile_avx2_256_16,
    parasail_sg_qx_striped_profile_avx2_256_32,
    // semi-global, both ends of s2 free.
    parasail_sg_dx_striped_profile_sse2_128_8,
    parasail_sg_dx_striped_profile_sse2_128_16,
    parasail_sg_dx_striped_profile_sse2_128_32,
    parasail_sg_dx_striped_profile_sse41_128_8,
    parasail_sg_dx_striped_profile_sse41_128_16,
    parasail_sg_dx_striped_profile_sse41_128_32,
    parasail_sg_dx_striped_profile_avx2_256_8,
    parasail_sg_dx_striped_profile_avx2_256_16,
    parasail_sg_dx_striped_profile_avx2_256_32,
);

/// parasail matrix type tag for a plain square substitution matrix.
pub const PARASAIL_MATRIX_TYPE_SQUARE: c_int = 0;
