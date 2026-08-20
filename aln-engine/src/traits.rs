//! The pairwise-aligner trait and its object-safe wrapper.

use aln_core::{Alignment, Sequence};

use crate::params::{AlignParams, AlignerCaps};
use crate::Result;

/// A full dynamic-programming pairwise aligner.
///
/// # Why there is a separate `prepare_subject` step
///
/// Striped SIMD aligners build a query profile — the substitution scores for one
/// sequence, laid out in SIMD lane order — and reuse it across every opposing
/// sequence.  GIRI models this with `SWAligner::getProfile`; parasail exposes it
/// as `parasail_profile_create_sat` plus the `_profile_` function family.
/// Rebuilding it per pair throws away most of the speedup in a one-to-many pass,
/// which is exactly the shape `acons` runs.
///
/// The **subject** is the prepared side, because the suite's hot loop is one
/// reference/consensus against many instances, and the reference is the subject
/// under this crate's conventions (see [`aln_core::matrix`]).
///
/// Backends with no profile concept can set `type Profile = Sequence` and make
/// `prepare_subject` a clone.
pub trait PairwiseAligner: Send + Sync {
    /// Backend-specific prepared form of the subject.
    type Profile: Send + Sync;

    /// Backend name, for diagnostics.
    fn name(&self) -> &'static str;

    /// What this backend supports.
    fn caps(&self) -> AlignerCaps;

    /// The parameters this aligner was configured with.
    fn params(&self) -> &AlignParams;

    /// Build the reusable prepared form of the subject.
    fn prepare_subject(&self, subject: &Sequence) -> Result<Self::Profile>;

    /// Align one query against a prepared subject.
    ///
    /// Returns `Ok(None)` when no alignment clears
    /// [`AlignParams::min_score`] — a normal outcome, not an error.
    fn align_prepared(&self, subject: &Self::Profile, query: &Sequence)
        -> Result<Option<Alignment>>;

    /// Align a single pair.  The default prepares the subject and discards it;
    /// override only if a backend has a cheaper one-shot path.
    fn align(&self, query: &Sequence, subject: &Sequence) -> Result<Option<Alignment>> {
        let prepared = self.prepare_subject(subject)?;
        self.align_prepared(&prepared, query)
    }

    /// Score a query against a prepared subject without building a traceback.
    ///
    /// Returns `Ok(None)` under the same conditions as
    /// [`align_prepared`](Self::align_prepared).
    ///
    /// The default runs the full alignment and takes its score, so it saves
    /// nothing.  Backends with separate score-only kernels should override it:
    /// parasail's skip the O(mn) traceback matrix entirely, which is the whole
    /// reason to ask for a score alone.  Use it for a cheap prepass when most
    /// pairs will be rejected by [`AlignParams::min_score`] anyway.
    fn score_prepared(&self, subject: &Self::Profile, query: &Sequence) -> Result<Option<i32>> {
        Ok(self.align_prepared(subject, query)?.map(|a| a.score))
    }

    /// Score a single pair without building a traceback.
    fn score(&self, query: &Sequence, subject: &Sequence) -> Result<Option<i32>> {
        let prepared = self.prepare_subject(subject)?;
        self.score_prepared(&prepared, query)
    }
}

/// Object-safe view of a [`PairwiseAligner`], for runtime backend selection.
///
/// `PairwiseAligner` has an associated type and so cannot be made into a trait
/// object.  This trait erases it, at the cost of not letting the caller hold a
/// profile across calls — [`align_one_to_many`](DynAligner::align_one_to_many)
/// therefore exists so the profile can still be reused internally.
///
/// Every `PairwiseAligner` implements this automatically:
///
/// ```ignore
/// let aligner: Box<dyn DynAligner> = match cli.backend {
///     Backend::Parasail => Box::new(ParasailAligner::new(matrix, params)?),
///     Backend::Farrar   => Box::new(FarrarAligner::new(matrix, params)?),
/// };
/// ```
pub trait DynAligner: Send + Sync {
    fn name(&self) -> &'static str;
    fn caps(&self) -> AlignerCaps;
    fn params(&self) -> &AlignParams;

    fn align(&self, query: &Sequence, subject: &Sequence) -> Result<Option<Alignment>>;

    /// Align many queries against one subject, preparing the subject once.
    ///
    /// Sequential — see [`crate::driver::one_to_many`] for the parallel form.
    fn align_one_to_many(
        &self,
        subject: &Sequence,
        queries: &[Sequence],
    ) -> Result<Vec<Alignment>>;
}

impl<A> DynAligner for A
where
    A: PairwiseAligner,
{
    fn name(&self) -> &'static str {
        PairwiseAligner::name(self)
    }

    fn caps(&self) -> AlignerCaps {
        PairwiseAligner::caps(self)
    }

    fn params(&self) -> &AlignParams {
        PairwiseAligner::params(self)
    }

    fn align(&self, query: &Sequence, subject: &Sequence) -> Result<Option<Alignment>> {
        PairwiseAligner::align(self, query, subject)
    }

    fn align_one_to_many(
        &self,
        subject: &Sequence,
        queries: &[Sequence],
    ) -> Result<Vec<Alignment>> {
        let prepared = self.prepare_subject(subject)?;
        let mut out = Vec::new();
        for q in queries {
            if q.is_empty() {
                continue;
            }
            if let Some(a) = self.align_prepared(&prepared, q)? {
                out.push(a);
            }
        }
        Ok(out)
    }
}
