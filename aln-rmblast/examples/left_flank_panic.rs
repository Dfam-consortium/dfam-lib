//! Minimal reproduction of an rmblast-lib bug, for reporting upstream.
//!
//! # Symptom
//!
//! `attempt to subtract with overflow` at
//! `rmblast-lib/src/search/gapped.rs:184` whenever a gapped alignment needs a
//! non-trivial **left extension** — i.e. the hit does not start at the very
//! beginning of the subject. Five bases of subject left-flank is enough.
//!
//! # Cause
//!
//! In the `REVERSE` (left-extension) pass the `b` pointer is initialised
//! unconditionally, before the loop whose bound would make it safe:
//!
//! ```ignore
//! let mut b_cur: *const u8 = if REVERSE {
//!     unsafe { b.as_ptr().add(n - 1 - first_b_index) }   // line 184
//! } else {
//!     std::hint::black_box(unsafe { b.as_ptr().add(first_b_index) })
//! };
//!
//! let mut b_idx = first_b_index;
//! while b_idx < inner_end {   // inner_end = b_size.min(n)
//! ```
//!
//! When `first_b_index >= n` the loop body never runs, but `n - 1 -
//! first_b_index` has already underflowed.
//!
//! # Impact
//!
//! * **Debug builds**: panics. Any debug-built consumer cannot align a hit that
//!   starts partway into a subject — the ordinary RepeatMasker shape.
//! * **Release builds**: the subtraction wraps to a huge `usize` and `.add()`
//!   forms an out-of-bounds pointer. It is never dereferenced (the loop is
//!   skipped), so results are correct, but constructing the pointer is itself
//!   undefined behaviour under Rust's rules and is not guaranteed to stay
//!   benign across compiler versions.
//!
//! # Suggested fix
//!
//! Compute `b_cur` only when the loop will execute, or clamp:
//!
//! ```ignore
//! let mut b_cur: *const u8 = if first_b_index >= inner_end {
//!     b.as_ptr()          // never read; the loop below does not run
//! } else if REVERSE {
//!     unsafe { b.as_ptr().add(n - 1 - first_b_index) }
//! } else {
//!     std::hint::black_box(unsafe { b.as_ptr().add(first_b_index) })
//! };
//! ```
//!
//! # Running
//!
//! ```sh
//! cargo run -p aln-rmblast --example left_flank_panic            # panics at lflank 5
//! cargo run -p aln-rmblast --example left_flank_panic --release  # all OK
//! ```

use aln_core::{Sequence, SubstMatrix};
use aln_engine::engine::{SearchEngine, SearchParams, SeqSource};
use aln_rmblast::{RmblastEngine, RmblastOptions};

const MATRIX: &str = "\
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

fn random_seq(r: &mut Rng, n: usize) -> Vec<u8> {
    (0..n).map(|_| b"ACGT"[r.below(4)]).collect()
}

fn engine() -> RmblastEngine {
    let params = SearchParams {
        matrix: Some(SubstMatrix::parse(MATRIX).unwrap()),
        gap_init: -25,
        ins_gap_ext: -5,
        del_gap_ext: -5,
        min_match: 7,
        min_score: 100,
        mask_level: 101,
        ..Default::default()
    };
    RmblastEngine::new(params, RmblastOptions::default()).unwrap()
}

fn main() {
    std::panic::set_hook(Box::new(|_| {}));
    let core = random_seq(&mut Rng::new(2), 400);

    println!(
        "400 bp query vs the same 400 bp preceded by N bases of left flank\n\
         (build: {})\n",
        if cfg!(debug_assertions) { "debug" } else { "release" }
    );

    for left_flank in 0..=8usize {
        let mut subject = random_seq(&mut Rng::new(77), left_flank);
        subject.extend_from_slice(&core);

        let e = engine();
        let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            e.search(
                &SeqSource::Memory(vec![Sequence::new("q", core.clone())]),
                &SeqSource::Memory(vec![Sequence::new("s", subject.clone())]),
            )
        }));
        println!(
            "  left flank {left_flank:>2} -> {}",
            match outcome {
                Ok(Ok(h)) => format!("ok, {} hit(s)", h.len()),
                Ok(Err(e)) => format!("error: {e}"),
                Err(_) => "PANIC (gapped.rs:184)".to_string(),
            }
        );
    }
}
