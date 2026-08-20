//! parasail SIMD backend.
//!
//! Wraps the striped-traceback kernels of [parasail] 2.6.2 (Daily, 2016) behind
//! [`PairwiseAligner`].  Only a curated subset is vendored — see `build.rs`.
//!
//! [parasail]: https://github.com/jeffdaily/parasail
//!
//! # Two inversions to keep straight
//!
//! parasail names its arguments `s1` (profiled) and `s2` (streamed), and treats
//! `s1` as the query.  This crate profiles the **subject**, because the hot loop
//! in `acons` is one consensus against many instances and the profile must be
//! reusable across them.  So:
//!
//! ```text
//!   parasail s1  =  aln-core subject   (profiled, "query" in parasail's docs)
//!   parasail s2  =  aln-core query     (streamed, "database" in parasail's docs)
//! ```
//!
//! Two consequences, both of which silently produce wrong answers if missed and
//! are therefore pinned by tests:
//!
//! 1. **The scoring matrix is transposed.**  parasail looks scores up as
//!    `matrix[s2][s1]`; `aln-core` defines them as `matrix[subject][query]`.
//!    With the mapping above those are transposes of one another, so
//!    [`ParasailMatrix::from_subst`] writes `pm[a][b] = subst[b][a]`.
//! 2. **The CIGAR opcodes are swapped** relative to SAM's convention.  In
//!    parasail's output `'I'` consumes `s1` alone and `'D'` consumes `s2`
//!    alone, so here `'I'` is a gap in the *query* and `'D'` a gap in the
//!    *subject* — the reverse of [`EditScript::from_cigar`], which must not be
//!    used on parasail output.
//!
//! # Saturation and ISA selection
//!
//! parasail's own `_sat` wrappers and cpuid dispatcher are not vendored (they
//! would drag in nearly the whole library). Instead the 8/16/32-bit fallback
//! chain runs in [`ParasailAligner::align_prepared`] and the SIMD level is
//! chosen by `is_x86_feature_detected!` in [`isa`].

pub mod ffi;
pub mod isa;

use std::os::raw::{c_char, c_int};
use std::ptr::NonNull;

use aln_core::align::{Alignment, EditOp, EditScript};
use aln_core::{Sequence, Strand, SubstMatrix};
use aln_engine::{AlignMode, AlignParams, AlignerCaps, EngineError, PairwiseAligner, Result};

use isa::{Isa, Kernel, FALLBACK_CHAIN};

const NAME: &str = "parasail";

// ── Matrix ────────────────────────────────────────────────────────────────────

/// An owned parasail matrix, transposed from an [`SubstMatrix`].
///
/// parasail allocates an `(n+1) x (n+1)` grid: one row and column past the
/// alphabet act as a `'*'` catch-all that every unmapped byte falls into.  That
/// catch-all is filled from the matrix's `N` row and column, which gives
/// unknown input the same treatment `aln-reference` gives it.
pub struct ParasailMatrix {
    raw: NonNull<ffi::parasail_matrix_t>,
    /// Kept alive because parasail stores the pointer, not a copy.
    _alphabet: std::ffi::CString,
}

// Read-only once built; the kernels take it as `const`.
unsafe impl Send for ParasailMatrix {}
unsafe impl Sync for ParasailMatrix {}

impl ParasailMatrix {
    /// Build the transposed parasail form of `subst`.
    ///
    /// Fails if the alphabet has no `N`, since there would be nothing sensible
    /// to put in the catch-all row.
    pub fn from_subst(subst: &SubstMatrix) -> Result<Self> {
        let n = subst.size();
        let n_idx = subst.index_of(b'N').ok_or_else(|| {
            EngineError::unsupported(
                NAME,
                format!(
                    "matrix alphabet {:?} has no 'N'; parasail needs one to fill \
                     its '*' catch-all row",
                    String::from_utf8_lossy(subst.alphabet())
                ),
            )
        })?;

        let alphabet = std::ffi::CString::new(subst.alphabet().to_vec())
            .map_err(|e| EngineError::backend(NAME, format!("alphabet has an interior NUL: {e}")))?;

        // SAFETY: alphabet is a valid NUL-terminated C string. The match /
        // mismatch values are placeholders; every cell is overwritten below.
        let raw = unsafe { ffi::parasail_matrix_create(alphabet.as_ptr(), 1, -1) };
        let raw = NonNull::new(raw)
            .ok_or_else(|| EngineError::backend(NAME, "parasail_matrix_create returned NULL"))?;

        // SAFETY: raw is a live matrix created just above.
        let size = unsafe { raw.as_ref().size } as usize;
        if size != n + 1 {
            unsafe { ffi::parasail_matrix_free(raw.as_ptr()) };
            return Err(EngineError::backend(
                NAME,
                format!("expected a {}x{} parasail matrix, got {size}x{size}", n + 1, n + 1),
            ));
        }
        // SAFETY: created with an explicit alphabet, so it is a user matrix and
        // parasail_matrix_set_value will accept writes.
        if unsafe { raw.as_ref().user_matrix }.is_null() {
            unsafe { ffi::parasail_matrix_free(raw.as_ptr()) };
            return Err(EngineError::backend(
                NAME,
                "parasail matrix is not writable (user_matrix is NULL)",
            ));
        }

        let star = n; // index of the '*' catch-all
        let set = |row: usize, col: usize, v: i32| {
            // SAFETY: row and col are both < size, checked above.
            unsafe { ffi::parasail_matrix_set_value(raw.as_ptr(), row as c_int, col as c_int, v) };
        };

        // pm[a][b] must equal subst.score(subject = b, query = a).
        for a in 0..n {
            for b in 0..n {
                set(a, b, subst.score_idx(b, a));
            }
            // Unknown subject byte: score it as if the subject were N.
            set(a, star, subst.score_idx(n_idx, a));
            // Unknown query byte: score it as if the query were N.
            set(star, a, subst.score_idx(a, n_idx));
        }
        set(star, star, subst.score_idx(n_idx, n_idx));

        Ok(ParasailMatrix { raw, _alphabet: alphabet })
    }

    fn as_ptr(&self) -> *const ffi::parasail_matrix_t {
        self.raw.as_ptr()
    }
}

impl Drop for ParasailMatrix {
    fn drop(&mut self) {
        // SAFETY: raw was produced by parasail_matrix_create and is freed once.
        unsafe { ffi::parasail_matrix_free(self.raw.as_ptr()) };
    }
}

// ── Profile ───────────────────────────────────────────────────────────────────

/// A subject prepared for repeated alignment.
///
/// Holds the subject bytes itself: `parasail_profile_new` stores the caller's
/// pointer rather than copying, so the buffer has to outlive the profile.  The
/// `Box<[u8]>` keeps the allocation at a stable address even when this struct
/// moves.
pub struct ParasailProfile {
    name: String,
    seq: Box<[u8]>,
    raw: NonNull<ffi::parasail_profile_t>,
    isa: Isa,
}

// The kernels take the profile as `const` and write only to per-call scratch,
// so concurrent alignment against one profile is sound.
unsafe impl Send for ParasailProfile {}
unsafe impl Sync for ParasailProfile {}

impl ParasailProfile {
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    pub fn name(&self) -> &str {
        &self.name
    }
}

impl Drop for ParasailProfile {
    fn drop(&mut self) {
        // SAFETY: raw came from a parasail profile constructor and is freed
        // once. `seq` outlives this call because it is dropped afterwards.
        unsafe { ffi::parasail_profile_free(self.raw.as_ptr()) };
    }
}

// ── Aligner ───────────────────────────────────────────────────────────────────

/// The parasail-backed aligner.
pub struct ParasailAligner {
    subst: SubstMatrix,
    matrix: ParasailMatrix,
    params: AlignParams,
    isa: Isa,
}

impl ParasailAligner {
    /// Build an aligner, detecting the widest SIMD level this CPU supports.
    pub fn new(subst: SubstMatrix, params: AlignParams) -> Result<Self> {
        Self::with_isa(subst, params, Isa::detect())
    }

    /// Build against a specific SIMD level.
    ///
    /// Intended for tests that want to prove every ISA agrees; production code
    /// should use [`new`](Self::new).  Requesting an ISA the CPU lacks is
    /// rejected here rather than faulting inside the kernel.
    pub fn with_isa(subst: SubstMatrix, params: AlignParams, isa: Isa) -> Result<Self> {
        params.validate()?;
        let available = Isa::detect();
        let ok = match isa {
            Isa::Sse2 => true,
            Isa::Sse41 => matches!(available, Isa::Sse41 | Isa::Avx2),
            Isa::Avx2 => matches!(available, Isa::Avx2),
        };
        if !ok {
            return Err(EngineError::unsupported(
                NAME,
                format!("this CPU does not support {}", isa.as_str()),
            ));
        }
        let matrix = ParasailMatrix::from_subst(&subst)?;
        Ok(ParasailAligner { subst, matrix, params, isa })
    }

    pub fn subst_matrix(&self) -> &SubstMatrix {
        &self.subst
    }

    pub fn isa(&self) -> Isa {
        self.isa
    }

    /// Map an [`AlignMode`] to a parasail kernel.
    ///
    /// **The semi-global mapping is inverted**, because parasail's `q` refers to
    /// `s1` — which is this crate's *subject*, not its query.  So free query
    /// ends select `sg_dx` (both ends of `s2` free) and free subject ends select
    /// `sg_qx`.  `sg_qx_selects_free_subject_ends` pins this.
    fn kernel(&self) -> Kernel {
        match self.params.mode {
            AlignMode::Local => Kernel::Sw,
            AlignMode::Global => Kernel::Nw,
            AlignMode::SemiGlobal { free_query_ends, free_subject_ends } => {
                match (free_query_ends, free_subject_ends) {
                    (true, true) => Kernel::Sg,
                    (true, false) => Kernel::SgDx,
                    (false, true) => Kernel::SgQx,
                    (false, false) => Kernel::Nw,
                }
            }
        }
    }

    fn free_ends(&self) -> (bool, bool) {
        match self.params.mode {
            AlignMode::Local => (true, true),
            AlignMode::Global => (false, false),
            AlignMode::SemiGlobal { free_query_ends, free_subject_ends } => {
                (free_query_ends, free_subject_ends)
            }
        }
    }
}

impl PairwiseAligner for ParasailAligner {
    type Profile = ParasailProfile;

    fn name(&self) -> &'static str {
        NAME
    }

    fn caps(&self) -> AlignerCaps {
        AlignerCaps {
            name: "parasail",
            modes: &[
                AlignMode::Local,
                AlignMode::Global,
                AlignMode::SemiGlobal { free_query_ends: true, free_subject_ends: true },
                AlignMode::SemiGlobal { free_query_ends: true, free_subject_ends: false },
                AlignMode::SemiGlobal { free_query_ends: false, free_subject_ends: true },
                AlignMode::SemiGlobal { free_query_ends: false, free_subject_ends: false },
            ],
            traceback: true,
            banded: false,
            // Reported as the compiled-in maximum; the instance's actual level
            // is available from `ParasailAligner::isa`.
            simd: "sse2/sse41/avx2",
        }
    }

    fn params(&self) -> &AlignParams {
        &self.params
    }

    fn prepare_subject(&self, subject: &Sequence) -> Result<ParasailProfile> {
        let seq: Box<[u8]> = subject.seq.clone().into_boxed_slice();
        if seq.is_empty() {
            return Err(EngineError::backend(NAME, "cannot profile an empty subject"));
        }
        if seq.len() > c_int::MAX as usize {
            return Err(EngineError::backend(
                NAME,
                format!("subject of {} bases exceeds parasail's int length", seq.len()),
            ));
        }

        let create = self.isa.profile_create();
        // SAFETY: `seq` is a live, non-empty buffer owned by the returned
        // struct, so it outlives the profile that borrows it. The matrix
        // outlives the profile because callers hold the aligner across use.
        let raw = unsafe {
            create(seq.as_ptr() as *const c_char, seq.len() as c_int, self.matrix.as_ptr())
        };
        let raw = NonNull::new(raw).ok_or_else(|| {
            EngineError::backend(NAME, "parasail profile constructor returned NULL")
        })?;

        Ok(ParasailProfile { name: subject.name.clone(), seq, raw, isa: self.isa })
    }

    fn align_prepared(
        &self,
        subject: &ParasailProfile,
        query: &Sequence,
    ) -> Result<Option<Alignment>> {
        let Some(result) = self.run_kernel(subject, query, Trace::Yes)? else {
            return Ok(None);
        };
        // From here on `result` must be freed on every path.
        let outcome = self.build_alignment(result, subject, query);
        // SAFETY: result is live and freed exactly once.
        unsafe { ffi::parasail_result_free(result.as_ptr()) };
        outcome
    }

    /// Uses parasail's score-only kernels, which never allocate the O(mn)
    /// traceback matrix.
    fn score_prepared(
        &self,
        subject: &ParasailProfile,
        query: &Sequence,
    ) -> Result<Option<i32>> {
        let Some(result) = self.run_kernel(subject, query, Trace::No)? else {
            return Ok(None);
        };
        // SAFETY: result is live; freed immediately after reading the score.
        let score = unsafe { ffi::parasail_result_get_score(result.as_ptr()) };
        unsafe { ffi::parasail_result_free(result.as_ptr()) };
        Ok((score >= self.params.min_score).then_some(score))
    }
}

/// Which kernel family to run.  parasail keeps them separate; the score-only
/// ones skip the traceback matrix entirely.
#[derive(Copy, Clone, PartialEq, Eq)]
enum Trace {
    Yes,
    No,
}

impl ParasailAligner {
    /// Run the kernel, widening lanes until one does not saturate.
    ///
    /// Returns the raw result, which the caller must free.  `Ok(None)` means
    /// there was nothing to align.
    fn run_kernel(
        &self,
        subject: &ParasailProfile,
        query: &Sequence,
        trace: Trace,
    ) -> Result<Option<NonNull<ffi::parasail_result_t>>> {
        if subject.isa != self.isa {
            return Err(EngineError::backend(
                NAME,
                format!(
                    "profile was built for {} but this aligner runs {} — profile \
                     layout depends on register width",
                    subject.isa.as_str(),
                    self.isa.as_str()
                ),
            ));
        }
        if query.is_empty() || subject.is_empty() {
            return Ok(None);
        }
        if query.len() > c_int::MAX as usize {
            return Err(EngineError::backend(
                NAME,
                format!("query of {} bases exceeds parasail's int length", query.len()),
            ));
        }

        let kernel = self.kernel();
        let open = self.params.gap_open as c_int;
        let ext = self.params.gap_extend as c_int;

        // Try progressively wider lanes until one does not overflow.
        let mut result = None;
        for width in FALLBACK_CHAIN {
            let f = match trace {
                Trace::Yes => isa::resolve(kernel, self.isa, width),
                Trace::No => isa::resolve_score_only(kernel, self.isa, width),
            };
            // SAFETY: profile and query are live for the call; the ISA matches
            // the kernel, checked above and at construction.
            let r = unsafe {
                f(
                    subject.raw.as_ptr(),
                    query.seq.as_ptr() as *const c_char,
                    query.len() as c_int,
                    open,
                    ext,
                )
            };
            let Some(r) = NonNull::new(r) else {
                return Err(EngineError::backend(
                    NAME,
                    format!("{:?} kernel at width {width:?} returned NULL", kernel),
                ));
            };
            // SAFETY: r is a live result.
            if unsafe { ffi::parasail_result_is_saturated(r.as_ptr()) } != 0 {
                unsafe { ffi::parasail_result_free(r.as_ptr()) };
                continue;
            }
            result = Some(r);
            break;
        }

        let Some(result) = result else {
            return Err(EngineError::backend(
                NAME,
                "alignment saturated even at 32-bit lanes",
            ));
        };
        Ok(Some(result))
    }
}

impl ParasailAligner {
    /// Convert a parasail result into an [`Alignment`].
    ///
    /// Does not free `result`; the caller owns it.
    fn build_alignment(
        &self,
        result: NonNull<ffi::parasail_result_t>,
        subject: &ParasailProfile,
        query: &Sequence,
    ) -> Result<Option<Alignment>> {
        // SAFETY: result is live for the duration of this function.
        let score = unsafe { ffi::parasail_result_get_score(result.as_ptr()) };
        if score < self.params.min_score {
            return Ok(None);
        }

        // SAFETY: both sequences are live and their lengths were bounds-checked.
        let cigar = unsafe {
            ffi::parasail_result_get_cigar(
                result.as_ptr(),
                subject.seq.as_ptr() as *const c_char,
                subject.len() as c_int,
                query.seq.as_ptr() as *const c_char,
                query.len() as c_int,
                self.matrix.as_ptr(),
            )
        };
        let Some(cigar) = NonNull::new(cigar) else {
            return Err(EngineError::backend(NAME, "parasail_result_get_cigar returned NULL"));
        };

        let built = self.decode_cigar(cigar, subject, query, score);
        // SAFETY: cigar is live and freed exactly once.
        unsafe { ffi::parasail_cigar_free(cigar.as_ptr()) };
        built
    }

    fn decode_cigar(
        &self,
        cigar: NonNull<ffi::parasail_cigar_t>,
        subject: &ParasailProfile,
        query: &Sequence,
        score: i32,
    ) -> Result<Option<Alignment>> {
        // SAFETY: cigar is live.
        let c = unsafe { cigar.as_ref() };
        if c.len < 0 || c.seq.is_null() {
            return Err(EngineError::backend(NAME, "parasail returned a malformed CIGAR"));
        }
        // SAFETY: c.seq points to c.len packed operations.
        let ops = unsafe { std::slice::from_raw_parts(c.seq, c.len as usize) };

        let mut script = EditScript::new();
        for &packed in ops {
            // SAFETY: decoding a value parasail itself produced.
            let op_char = unsafe { ffi::parasail_cigar_decode_op(packed) } as u8;
            let run = unsafe { ffi::parasail_cigar_decode_len(packed) };
            let op = match op_char {
                b'M' | b'=' | b'X' => EditOp::Sub,
                // parasail's 'I' consumes s1 alone. s1 is our subject, so the
                // query is what has the gap.  This is the reverse of SAM.
                b'I' => EditOp::GapInQuery,
                // 'D' consumes s2 alone — our query — so the subject is gapped.
                b'D' => EditOp::GapInSubject,
                other => {
                    return Err(EngineError::backend(
                        NAME,
                        format!("unexpected CIGAR operator {:?}", other as char),
                    ))
                }
            };
            script.push(op, run);
        }

        if script.is_empty() {
            return Ok(None);
        }

        // beg_query indexes s1 (subject); beg_ref indexes s2 (query).
        let mut subj_start = c.beg_query.max(0) as usize;
        let mut query_start = c.beg_ref.max(0) as usize;

        let (free_q, free_s) = self.free_ends();
        trim_free_ends(&mut script, &mut query_start, &mut subj_start, free_q, free_s);

        if script.is_empty() {
            return Ok(None);
        }

        let mut aln = Alignment::new(
            query.name.clone(),
            subject.name.clone(),
            query_start,
            subj_start,
            Strand::Plus,
            score,
            script,
        );
        aln.query_len = Some(query.len());
        aln.subj_len = Some(subject.len());

        // The FFI boundary is exactly where a coordinate slip would go unnoticed,
        // so the invariant is checked rather than assumed.
        aln.validate().map_err(|e| {
            EngineError::backend(NAME, format!("parasail produced an inconsistent alignment: {e}"))
        })?;

        if self.needs_traceback_check() {
            self.verify_traceback(&aln, subject, query)?;
        }
        Ok(Some(aln))
    }

    /// Whether this mode's traceback needs checking against its own score.
    ///
    /// Measured over 400 randomised pairs per mode (see
    /// `tests/traceback_consistency.rs`): `Local`, `Global` and the two
    /// one-sided semi-global modes were self-consistent in every case, while
    /// `SemiGlobal { free_query_ends: true, free_subject_ends: true }` — which
    /// maps to parasail's plain `sg` — disagreed in about 2.5%.  In every such
    /// case the *score* still matched the reference aligner exactly; it is only
    /// the reconstructed path that is wrong, typically ending short of both
    /// sequences, which is not a legal semi-global endpoint.
    ///
    /// The check is therefore scoped to that one mode rather than paid for on
    /// every alignment.
    fn needs_traceback_check(&self) -> bool {
        matches!(
            self.params.mode,
            AlignMode::SemiGlobal { free_query_ends: true, free_subject_ends: true }
        )
    }

    /// Re-score the reported traceback and fail if it does not reproduce the
    /// reported score.
    ///
    /// Returning an error is deliberate: a silently wrong alignment is far worse
    /// than a loud one, and the caller can either fall back to `aln-reference`
    /// or use a one-sided semi-global mode, both of which are unaffected.
    fn verify_traceback(
        &self,
        aln: &Alignment,
        subject: &ParasailProfile,
        query: &Sequence,
    ) -> Result<()> {
        let (gq, gs) = aln.gapped(&query.seq, &subject.seq).map_err(|e| {
            EngineError::backend(NAME, format!("cannot expand parasail traceback: {e}"))
        })?;
        let rp = aln_core::stats::RescoreParams {
            gap_open: -(self.params.gap_open as i32),
            ins_gap_extend: -(self.params.gap_extend as i32),
            del_gap_extend: -(self.params.gap_extend as i32),
            ..aln_core::stats::RescoreParams::new(&self.subst)
        };
        let rescored = aln_core::stats::rescore(&gq, &gs, &rp)
            .map_err(|e| EngineError::backend(NAME, format!("cannot rescore traceback: {e}")))?;
        if rescored.score != aln.score {
            return Err(EngineError::backend(
                NAME,
                format!(
                    "parasail's semi-global traceback is inconsistent with its own score \
                     (reported {}, path rescores to {}; cigar {}). The score is reliable — \
                     only the reconstructed path is not. Use a one-sided semi-global mode \
                     or aln-reference if you need the alignment itself.",
                    aln.score,
                    rescored.score,
                    aln.edits.to_cigar()
                ),
            ));
        }
        Ok(())
    }
}

/// Strip leading and trailing gap runs that sit on a free end.
///
/// For semi-global and global modes parasail's CIGAR spans both sequences in
/// full, including the unaligned tails.  A run on a *free* end is not part of
/// the alignment and is removed, advancing the corresponding start offset;
/// a run on a penalised end is real and is kept.  Local alignments arrive
/// already trimmed, so this is a no-op for them.
fn trim_free_ends(
    script: &mut EditScript,
    query_start: &mut usize,
    subj_start: &mut usize,
    free_query_ends: bool,
    free_subject_ends: bool,
) {
    // A leading GapInQuery consumes subject bases — trimmable iff the subject's
    // ends are free.  A leading GapInSubject consumes query bases.
    while let Some(&(op, run)) = script.ops.first() {
        match op {
            EditOp::GapInQuery if free_subject_ends => *subj_start += run as usize,
            EditOp::GapInSubject if free_query_ends => *query_start += run as usize,
            _ => break,
        }
        script.ops.remove(0);
    }
    while let Some(&(op, _)) = script.ops.last() {
        let trimmable = match op {
            EditOp::GapInQuery => free_subject_ends,
            EditOp::GapInSubject => free_query_ends,
            EditOp::Sub => false,
        };
        if !trimmable {
            break;
        }
        script.ops.pop();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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

    fn matrix() -> SubstMatrix {
        SubstMatrix::parse(M14P35G).unwrap()
    }

    fn aligner(mode: AlignMode) -> ParasailAligner {
        let p = AlignParams { mode, gap_open: 25, gap_extend: 5, min_score: 1, ..Default::default() };
        ParasailAligner::new(matrix(), p).unwrap()
    }

    fn seq(name: &str, s: &[u8]) -> Sequence {
        Sequence::new(name, s.to_vec())
    }

    #[test]
    fn matrix_is_transposed_relative_to_parasail() {
        // 14p35g scores G/A as -7 and A/G as -10, so a transposition slip is
        // detectable from the score alone.
        let a = aligner(AlignMode::Local);
        // subject G, query A  =>  aln-core score_idx(G, A) = -7.
        // Embed it between matches so the local alignment spans all three.
        let aln = a
            .align(&seq("q", b"CAC"), &seq("s", b"CGC"))
            .unwrap()
            .unwrap();
        assert_eq!(aln.score, 12 + -7 + 12, "matrix appears transposed the wrong way");
    }

    #[test]
    fn identical_sequences_score_the_diagonal_sum() {
        let a = aligner(AlignMode::Local);
        let aln = a.align(&seq("q", b"ACGT"), &seq("s", b"ACGT")).unwrap().unwrap();
        assert_eq!(aln.score, 40); // 8 + 12 + 12 + 8
        assert_eq!(aln.edits.to_cigar(), "4M");
    }

    #[test]
    fn local_alignment_reports_offsets_into_both_sequences() {
        let a = aligner(AlignMode::Local);
        let query = seq("q", b"TTTTACGTACGTTTTT");
        let subject = seq("s", b"ACGTACGT");
        let aln = a.align(&query, &subject).unwrap().unwrap();
        assert_eq!((aln.query_start, aln.query_end), (4, 12));
        assert_eq!((aln.subj_start, aln.subj_end), (0, 8));
        aln.validate().unwrap();
    }

    #[test]
    fn cigar_opcodes_are_not_read_as_sam() {
        // The query carries an extra base, so the subject must be the gapped
        // side.  Reading parasail's 'I'/'D' with SAM semantics would put the
        // gap on the query and this assertion would fail.
        let a = aligner(AlignMode::Global);
        let query = seq("q", b"ACGTTACGT");
        let subject = seq("s", b"ACGTACGT");
        let aln = a.align(&query, &subject).unwrap().unwrap();
        assert_eq!(aln.edits.query_consumed(), 9);
        assert_eq!(aln.edits.subject_consumed(), 8);
        assert!(
            aln.edits.ops.iter().any(|&(op, _)| op == EditOp::GapInSubject),
            "expected the subject to carry the gap, got {}",
            aln.edits.to_cigar()
        );
    }

    #[test]
    fn global_mode_spans_both_sequences() {
        let a = aligner(AlignMode::Global);
        let query = seq("q", b"ACGTAC");
        let subject = seq("s", b"ACGTGC");
        let aln = a.align(&query, &subject).unwrap().unwrap();
        assert_eq!((aln.query_start, aln.query_end), (0, 6));
        assert_eq!((aln.subj_start, aln.subj_end), (0, 6));
    }

    #[test]
    fn sg_qx_selects_free_subject_ends() {
        // Free subject ends: the whole query must be used, the subject's flanks
        // are free.  If the sg_qx/sg_dx mapping were inverted this would trim
        // the query instead and the spans would come out swapped.
        let a = aligner(AlignMode::SemiGlobal {
            free_query_ends: false,
            free_subject_ends: true,
        });
        let query = seq("q", b"ACGTACGT");
        let subject = seq("s", b"TTTTTTACGTACGTTTTTTT");
        let aln = a.align(&query, &subject).unwrap().unwrap();
        assert_eq!((aln.query_start, aln.query_end), (0, 8), "query should be fully used");
        assert_eq!((aln.subj_start, aln.subj_end), (6, 14), "subject flanks should be free");
    }

    #[test]
    fn sg_dx_selects_free_query_ends() {
        let a = aligner(AlignMode::SemiGlobal {
            free_query_ends: true,
            free_subject_ends: false,
        });
        let query = seq("q", b"TTTTTTACGTACGTTTTTTT");
        let subject = seq("s", b"ACGTACGT");
        let aln = a.align(&query, &subject).unwrap().unwrap();
        assert_eq!((aln.subj_start, aln.subj_end), (0, 8), "subject should be fully used");
        assert_eq!((aln.query_start, aln.query_end), (6, 14), "query flanks should be free");
    }

    #[test]
    fn unknown_bytes_land_in_the_catch_all_row() {
        let a = aligner(AlignMode::Local);
        // '@' is outside the alphabet; the '*' row is filled from N, whose
        // score against G is -1.
        let aln = a.align(&seq("q", b"AC@T"), &seq("s", b"ACGT")).unwrap().unwrap();
        assert_eq!(aln.score, 8 + 12 - 1 + 8);
    }

    #[test]
    fn a_matrix_without_n_is_rejected() {
        let m = SubstMatrix::parse("  A   C\n  1  -1\n -1   1\n").unwrap();
        assert!(ParasailAligner::new(m, AlignParams::default()).is_err());
    }

    #[test]
    fn min_score_suppresses_weak_alignments() {
        let p = AlignParams { mode: AlignMode::Local, min_score: 1000, ..Default::default() };
        let a = ParasailAligner::new(matrix(), p).unwrap();
        assert!(a.align(&seq("q", b"ACGT"), &seq("s", b"ACGT")).unwrap().is_none());
    }

    #[test]
    fn empty_input_yields_no_alignment() {
        let a = aligner(AlignMode::Local);
        assert!(a.align(&seq("q", b""), &seq("s", b"ACGT")).unwrap().is_none());
    }

    #[test]
    fn a_profile_can_be_reused_across_queries() {
        let a = aligner(AlignMode::Local);
        let subject = seq("s", b"ACGTACGTACGT");
        let profile = a.prepare_subject(&subject).unwrap();
        for q in [b"ACGTACGT".as_slice(), b"GTACGTAC".as_slice(), b"ACGT".as_slice()] {
            let query = seq("q", q);
            let reused = a.align_prepared(&profile, &query).unwrap();
            let fresh = a.align(&query, &subject).unwrap();
            assert_eq!(reused.map(|x| x.score), fresh.map(|x| x.score));
        }
    }

    #[test]
    fn score_only_kernels_agree_with_the_traceback_kernels() {
        // The whole point of the score-only path is that it is the same
        // computation without the traceback matrix.  If the two disagree, the
        // prepass would reject alignments the real pass would have kept.
        let mut rng = 0x5eed_u64;
        let mut next = || {
            rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; rng
        };
        let mut rand_seq = |n: usize| {
            Sequence::new("s", (0..n).map(|_| b"ACGT"[(next() % 4) as usize]).collect::<Vec<u8>>())
        };
        for mode in [
            AlignMode::Local,
            AlignMode::Global,
            AlignMode::SemiGlobal { free_query_ends: true, free_subject_ends: true },
            AlignMode::SemiGlobal { free_query_ends: true, free_subject_ends: false },
            AlignMode::SemiGlobal { free_query_ends: false, free_subject_ends: true },
        ] {
            let p = AlignParams { mode, gap_open: 25, gap_extend: 5, min_score: i32::MIN, ..Default::default() };
            let al = ParasailAligner::new(matrix(), p).unwrap();
            for n in [17usize, 64, 250] {
                let a = rand_seq(n);
                let b = rand_seq(n + 11);
                let full = al.align(&b, &a).unwrap().map(|r| r.score);
                let quick = al.score(&b, &a).unwrap();
                assert_eq!(full, quick, "{mode:?} at len {n}: traceback kernel and \
                                         score-only kernel disagree");
            }
        }
    }

    #[test]
    fn a_long_alignment_survives_the_saturation_fallback() {
        // Well past what 8-bit lanes can accumulate (score ~ 10 per base),
        // so this only succeeds if the fallback to wider lanes works.
        let a = aligner(AlignMode::Local);
        let s: Vec<u8> = b"ACGTACGTAC".iter().cycle().take(4000).copied().collect();
        let aln = a.align(&seq("q", &s), &seq("s", &s)).unwrap().unwrap();
        assert!(aln.score > 30_000, "score was {}", aln.score);
        assert_eq!(aln.edits.query_consumed(), 4000);
    }

    #[test]
    fn every_supported_isa_agrees() {
        let query = seq("q", b"ACGTTTACGTACGATCGATCGAAA");
        let subject = seq("s", b"ACGTACGTACGATCGATCGTAA");

        let mut scores = Vec::new();
        for isa in [Isa::Sse2, Isa::Sse41, Isa::Avx2] {
            let p = AlignParams {
                mode: AlignMode::Local,
                gap_open: 25,
                gap_extend: 5,
                min_score: 1,
                ..Default::default()
            };
            match ParasailAligner::with_isa(matrix(), p, isa) {
                Ok(a) => {
                    let aln = a.align(&query, &subject).unwrap().unwrap();
                    scores.push((isa, aln.score));
                }
                // The CPU lacks this level; nothing to compare.
                Err(_) => continue,
            }
        }
        assert!(!scores.is_empty(), "no ISA was available");
        let first = scores[0].1;
        for (isa, score) in &scores {
            assert_eq!(*score, first, "{isa:?} disagreed");
        }
    }
}
