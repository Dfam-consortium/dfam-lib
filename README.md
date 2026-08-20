# dfam-lib

A Rust alignment stack for the `acons` / `autocons` / `amaln` port, designed to
be adopted by `dfam-curator` without a second, competing set of core types.

Status: **foundation complete and validated; parasail SIMD backend working.**

```
dfam-lib/
├── aln-core/        sequences, matrices, alignments, divergence, MSA,
│                    consensus calling, output formats, FASTA I/O
├── aln-engine/      the two traits + parallel drivers
├── aln-parasail/    parasail 2.6.2 SIMD backend (vendored subset)
├── aln-rmblast/     SearchEngine over the rmblast Rust port (excluded — see below)
├── aln-reference/   plain O(mn) Gotoh — the correctness arbiter
└── autocons/        the first tool: N best consensus sequences
```

`aln-rmblast` is **excluded from the workspace** and stands alone
(`aln-rmblast/Cargo.toml` carries its own `[workspace]` table). It depends on
[RMBlast](https://github.com/Dfam-consortium/RMBlast) — currently a private
repository — pinned to tag `3.0.2`. An unresolvable dependency fails *workspace
resolution*, so were it a member, nobody without access to that repository could
build any crate here. As it stands the workspace loads from a bare clone, and
`cargo build` inside `aln-rmblast/` is a separate, opt-in step that needs SSH
access to RMBlast (`aln-rmblast/.cargo/config.toml` sets
`net.git-fetch-with-cli` so your normal git credentials apply). To develop
against a local `rmblast-port` checkout instead, swap the `git`/`tag` dependency
for the commented-out `path` line beneath it — a `[patch]` section will not do,
because cargo loads the original source before applying a patch.

## The three conventions

Everything else follows from these. Each is documented at length in its module.

| | rule | why it bites |
|---|---|---|
| **Coordinates** | 0-based, half-open, forward-strand, `start < end` on both sides | `rmblast-lib` is 0-based half-open; `dfam-curator`, RepeatMasker and BLAST tabular are 1-based closed. Conversion happens *only* at I/O boundaries (`Alignment::query_one_based`). |
| **Gap vs. padding** | `-` means present-but-deleted; `' '` means not present | GIRI uses `<`/`>` for padding instead; `aln_core::seq` converts. Counting padding as a gap silently corrupts coverage and divergence denominators. |
| **Matrix orientation** | `matrix[subject][query]`, subject = consensus, query = genomic | Arian Smit's matrices are **not symmetric**. `14p35g` scores `G/A` as −7 and `A/G` as −10. A swapped lookup changes scores without failing. |

## Two traits, not one

`parasail` and `rmblast` do not belong behind the same interface:

- **`PairwiseAligner`** — two sequences in, one alignment out. Full DP, no
  seeding, no database. parasail, Farrar's striped SSE2, Monardo, SWAT. The
  analogue of GIRI's `SWAligner`, and the **only** trait `acons`/`autocons` need.
- **`SearchEngine`** — query + subject database in, many HSPs out, with seeding,
  masking and score cutoffs. rmblast, crossmatch, HMMER. The analogue of
  RepeatMasker's `SearchEngineI`.

`PairwiseAligner` has a `prepare_subject` step because striped SIMD aligners
build a reusable profile (GIRI `getProfile`, parasail
`parasail_profile_create_sat`); rebuilding it per pair throws away most of the
speedup in the one-to-many pass `acons` actually runs. Associated types are not
object-safe, so `DynAligner` (blanket-implemented) exists for CLI-time backend
selection.

`aln_engine::driver` replaces `ThreadedAligner`: `one_to_many`, `all_vs_all`,
`align_pairs`, on rayon, with **no global state** and deterministic output
ordering. GIRI's `MultipleAlignment::setAligner`,
`PairwiseAlignment::setScoreMatrix` and `ThreadedAligner::isMultithreaded` are
all static; none of that survives here.

## The insertion policy

Merging N pairwise alignments into one MSA forces a choice the two lineages make
differently, and it is not cosmetic:

- `InsertionPolicy::Drop` — `dfam-curator`'s `hits_to_multialign`. Width stays
  equal to the reference; insertion evidence is destroyed.
- `InsertionPolicy::GrowIncremental` — a port of GIRI's `adjustReference`. Merges
  each member into an accumulating reference left to right, **reusing** gap
  columns earlier members opened. Narrower, and it coalesces one indel event the
  aligner placed inconsistently — but it can imply homology that is not there.
  Order-sensitive by construction. **The default in `autocons`**, so output is
  comparable with `bin/autocons`.
- `InsertionPolicy::GrowPerSlot` — opens, at each inter-base slot, as many
  columns as the largest single insertion seen there. Members stay independent,
  so a shared column always means positional homology. Slightly wider.

`assemble_msa` takes it as an argument so the choice is visible at the call site;
`autocons --insertions {incremental,per-slot,drop}` exposes it.

One subtlety when porting `adjustReference`: where **both** the accumulated
reference and the member's row carry a gap, the C++'s `*sptr1 != *sptr2` test is
false and the existing column is *reused*. Opening a fresh column there instead
made a test alignment 543 columns wide rather than 539.

## Validation

Every port is pinned against the original, not just unit-tested:

| test | pins against |
|---|---|
| `aln-core/tests/matrix_lambda.rs` | `Matrix.pm::_calculateLambda`, every shipped crossmatch matrix |
| `aln-core/tests/rescore_perl_parity.rs` | `SearchResult.pm::rescoreAlignment`, 11 cases incl. CpG scoring + complexity adjustment |
| `aln-core/tests/divergence_perl_parity.rs` | `calcKimuraDivergence` + `calcK2PGapDivergence`, 13 cases |
| `aln-core/tests/format_perl_parity.rs` | `_toCrossMatchFormat` + `_toOUTFileFormat`, byte for byte |
| `aln-core/tests/caf_perl_parity.rs` | `_toCAF` byte for byte; pins `_toCIGAR`'s broken output |

Reference values were generated by running the Perl directly; each test file
carries the script that regenerates it. The lambda test skips itself when
`/usr/local/RepeatMasker` is absent; the others are self-contained.

### Perl behaviours reproduced deliberately

These are faithful, not bugs introduced in translation:

- **K2P saturation.** `calcKimuraDivergence` returns its `100.00` initialiser;
  `dfam-curator` returns `NaN`. Here the value is `Option<f64>`, with
  `Divergence::or_repeatmasker_default` for byte-identical Perl behaviour.
- **K2P-Gap saturation returns `-50.0`, not a sentinel** — it substitutes a
  literal `1` for the inner log term rather than bailing out. Its `100.00`
  initialiser survives only when there are zero well-characterised bases, and is
  then multiplied by 100, so callers see `10000`.
- **CpG site counts are inconsistent between routines.** The two divergence
  functions count them only when `div_cpg_mod` is on; `rescoreAlignment` counts
  them unconditionally. Same alignment, different `cpg_sites`.
- **The CpG discount ladder uses exact float equality.** A `prev_trans` of `1.1`
  — reachable in overlapping CpG runs — matches neither `== 2` nor `== 1` and
  passes through untouched.
- **xDrop re-bases before it splits.** The `adjScore < 0` reset is evaluated
  before the drop test, so an alignment decaying through zero re-bases rather
  than emitting a fragment at its peak.

### One deliberate departure

`aln-reference` uses the standard `H`/`E`/`F` recurrence, where **both** gap
states are opened from `H`. An earlier three-state version opened gaps only from
the match state, which forbids an insertion abutting a deletion and scored a
minority of alignments 1–2 points lower. That is a real variant but not the one
parasail, Farrar, crossmatch or BLAST implement — and since this module's job is
to arbitrate the others, it follows the mainstream recurrence. The differential
suite caught this on the first run.

## GIRI's striped traceback

`docs/giri-farrar-traceback.md` documents how GIRI gets a traceback out of a
striped SIMD aligner — Farrar's published algorithm is score-only — by storing
one `short` per cell encoding a whole gap **run length** rather than a direction,
so traceback decompresses runs instead of stepping cells.

It is a real trade-off, not an accident: two bytes per cell against parasail's
separate `_trace_` kernels and their `O(mn)` tables. The cost is that the
reconstructed path is not guaranteed to be the maximum-scoring one, because the
recorded run lengths track the current best `E`/`F` state rather than the run the
optimal path takes, and the lazy-F pass rewrites `F` afterwards.

Measured on a real pair: `autocons` reported **918** while the alignment it
emitted is worth **364** under the same matrix. So `autocons` ranks candidate
references by DP maxima while building its MSA from the reconstructed paths —
ranking and content come from different alignments.

## Comparing backends

Compare **scores** exactly. Do *not* require identical tracebacks — on ties,
which path a backend reports depends on its tie-breaking order, and striped
implementations legitimately differ from a row-major scalar one. A traceback is
correct if re-scoring it reproduces the reported score, which is what
`aln_core::stats::rescore` is for. `aln-reference` exists to be that arbiter.

`aln-parasail/tests/differential.rs` runs 200 randomised pairs per mode against
`aln-reference`, asserting exact score equality and self-consistent tracebacks on
both sides.

## aln-parasail

Vendors ~30 of parasail's 595 generated `.c` files — the striped-traceback
kernels for `sw`/`nw`/`sg` across SSE2/SSE4.1/AVX2 at 8/16/32-bit lanes, plus
allocation and CIGAR support. Deliberately **not** vendored:

- `satcheck.c` (2.1 MB) — the `_sat` wrappers. Saturation fallback (8 → 16 → 32
  bit) runs in Rust instead, which avoids pulling in nearly the whole library.
- `sw_dispatch.c` and the cpuid dispatcher — ISA selection uses
  `is_x86_feature_detected!`, which is already correct about OS-level AVX state.

Each ISA compiles to its own static library with its own `-m` flag, so the
compiler can never hoist a wider instruction into a narrower kernel.
`build.rs` drives `cc` directly; CMake is not required.

### Two inversions the code has to get right

parasail names the profiled sequence `s1` and treats it as the query. This crate
profiles the **subject** (the reusable side in `acons`'s one-to-many loop), so
`s1 = subject` and `s2 = query`. Therefore:

1. **The matrix is transposed** — parasail looks up `matrix[s2][s1]`, `aln-core`
   defines `matrix[subject][query]`.
2. **The CIGAR opcodes are swapped** relative to SAM — parasail's `'I'` consumes
   `s1` alone, so here it is a gap in the *query*. `EditScript::from_cigar` must
   not be used on parasail output.
3. **The semi-global variant names invert** — parasail's `q` is our subject, so
   free *query* ends select `sg_dx` and free *subject* ends select `sg_qx`.

All three are pinned by tests that fail if the mapping is reversed.

### Known limitation: all-ends-free semi-global traceback

Measured over 400 randomised pairs per mode (`tests/traceback_consistency.rs`):

| mode | score matches reference | traceback self-consistent |
|---|---|---|
| `Local` | 400/400 | 400/400 |
| `Global` | 400/400 | 400/400 |
| `SemiGlobal` free subject ends | 400/400 | 400/400 |
| `SemiGlobal` free query ends | 400/400 | 400/400 |
| `SemiGlobal` **all ends free** | 400/400 | ~390/400 |

parasail's plain `sg` returns an exact score, but `parasail_result_get_cigar`
sometimes reconstructs a path that score did not come from — typically one
ending short of *both* sequences, which is not a legal semi-global endpoint.
`ParasailAligner` re-scores the traceback in that mode only and returns an error
rather than handing back a wrong alignment. Use a one-sided semi-global mode, or
`aln-reference`, if you need the path itself. `Local` — what `acons` uses — is
unaffected.

## aln-rmblast

The `SearchEngine` side: seeded database search over the Rust rmblast port.
`Hsp` already uses 0-based half-open offsets with plus-strand subject
coordinates, so converting to `Alignment` is a field copy. Two things are not:

**The matrix must be transposed.** RepeatMasker ships each matrix twice —
`Matrices/crossmatch/14p35g.matrix` and `Matrices/ncbi/nt/14p35g.matrix` are
transposes of one another, because crossmatch indexes `matrix[subject][query]`
and NCBI indexes `matrix[query][subject]`. Each is correct when fed to its own
tool; converting between the in-memory forms means transposing.
`aln_rmblast::matrix` renders `SubstMatrix` into NCBI layout and lets rmblast's
own reader parse it, so symbol mapping, gap-slot handling and lambda estimation
all happen exactly as they would for a file on the command line. rmblast also
only recognises frequencies on a *comment* line (`# FREQS`), so the bare `FREQS`
line the crossmatch files carry must be commented on the way through.

**Gap costs use a different convention.**

```text
crossmatch / aln_core::stats::rescore :  open + (k-1) * extend
NCBI / rmblast                        :  open +  k    * extend
```

So a crossmatch `gap_init` of −25 with `gap_ext` −5 becomes rmblast's `20/5` —
the same pair RepeatMasker and `dfam-curator` pass on the command line.
`reported_scores_survive_rescoring` is the end-to-end check: it re-scores every
returned HSP under `aln-core`'s model and requires the reported score back, so a
slip in the transpose, the gap conversion or the edit-script mapping all surface
in one place.

### Upstream bug: left extension underflows

`rmblast-lib/src/search/gapped.rs:184` initialises the `b` pointer
unconditionally in the left-extension (`REVERSE`) pass, before the loop bound
that would make it safe:

```rust
unsafe { b.as_ptr().add(n - 1 - first_b_index) }
```

When `first_b_index >= n` the `usize` subtraction underflows. It panics under
`debug_assertions`, and in release wraps to an out-of-range pointer that is never
dereferenced — so release results are correct, but forming the pointer is
undefined behaviour and is not guaranteed to stay benign.

Trigger: any hit that does not start at the very beginning of the subject.
**Five bases of subject left-flank is enough** — the ordinary RepeatMasker shape.
`cargo run -p aln-rmblast --example left_flank_panic` reproduces it and suggests
a one-line fix. Two tests are `#[cfg_attr(debug_assertions, ignore)]` for this
reason; they pass under `cargo test -p aln-rmblast --release`.

## Output formats

`aln_core::fmt` ports `SearchResult.pm`'s two writers. Both are compared
byte for byte — trailing whitespace included — against fixtures generated by the
Perl itself (`aln-core/tests/data/`), because column alignment is the entire
point of these formats.

- `to_out_line` — the `.out` annotation line (`_toOUTFileFormat`).
- `to_crossmatch` — annotation line plus the wrapped alignment block
  (`_toCrossMatchFormat`), in `NoAlign` / `WithQuerySeq` / `WithSubjSeq` modes.
- `to_caf` — Compressed Alignment Format (`_toCAF`).
- `to_cigar_record` — the CAF record header with a CIGAR payload (`_toCIGAR`),
  with one deliberate divergence; see below.

`aln_core::SearchResult` carries the reporting annotation (percent divergence,
class/family, cluster id, Kimura) that `Alignment` deliberately omits, since
`autocons` holds `O(n²)` alignments at once and must stay lean.

Two faithfully-reproduced quirks:

- **Minus-strand coordinates print descending**, with the remaining-base count
  moved to the front: `C name (left) end begin`. `Alignment` always stores the
  span ascending; the reordering happens only in the writer.
- **The gap-opening counter reads `$qChars[($j - 1) | 0]`.** At `j == 0` that is
  index −1, which in Perl is the *last* character of the current 50-column
  block, not the previous one. Reproduced so gap counts match exactly.

- **CAF leaves a trailing gap run unclosed.** The Perl's loop simply ends, so
  `ACGTACGT--` against `ACGTACGTAC` yields `ACGTACGT+AC` with no closing `+`.

There is deliberately **no `out_header()`**: RepeatMasker's `.out` header comes
from `ProcessRepeats`, not `SearchResult.pm`, and its column widths are computed
from the widest value present in the file. A fixed header string would line up
only by coincidence.

### The one deliberate divergence: `_toCIGAR` is broken

`to_cigar_record` does **not** reproduce `_toCIGAR`'s output, because that output
is unusable. The Perl pairs each run's length with the *following* run's
operator and opens with a zero-length run:

```text
query   AAGACTT---A
subject AAT--CTAATA

  _toCIGAR emits:            0M3D2M2I3M1M
  its own doc comment says:  3M2D2M3I1M      <- what we emit
```

Ten aligned columns with no gaps come out of the Perl as `0M10M`. Nothing can
parse that, so this emits what RepeatMasker documents rather than what it does.
`caf_perl_parity.rs` asserts *both* strings, so the divergence stays visible and
will surface if the Perl is ever fixed. The shared record header is still byte
for byte.

Separately, RepeatMasker's CIGAR uses `I` for a gap in the **query** and `D` for
a gap in the **subject** — the inverse of SAM, and therefore of
`EditScript::to_cigar`. That is the format's definition rather than a defect, so
it is preserved; a test asserts the two encodings stay mirror images.

**PSL is not implemented** — RepeatMasker declares the constant but
`toStringFormatted` croaks on it, so there is nothing to port.

## autocons

The first tool on the stack, a port of `acons/src/autocons.cpp`.

```sh
autocons family.fa --format fasta -n 3 --aln out
```

Two phases:

1. **Reference selection** — every input sequence is tried as the reference:
   align the rest to it, assemble an MSA, call a consensus. Candidates rank by
   **total alignment score**, not by any property of the consensus
   (`MultipleAlignment::autoConsensus`'s `outscore` is `construct`'s return).
2. **Refinement** — each of the top N consensi is re-aligned against the whole
   input and re-called, up to `--iterations` further passes, stopping early once
   the consensus stops changing.

One subtlety, faithful to the C++: **the two phases treat the reference row
differently.** Phase 1 calls with `withRefSeq = true` so the reference *is*
counted in the profile; phase 2 does `maln.erase(maln.begin())` first so it is
*not*. Phase 1's consensus is the starting point for phase 2, so the difference
is load-bearing.

Reference selection skips self-comparison **by index**. The C++ skips by pointer
identity, which silently fails to skip a duplicated sequence held as a distinct
object.

### Parallelism, and how to size a run

Rayon throughout; no MPI. The C++'s three modes (`process_locally`,
`process_with_pthreads`, `process_with_mpi`) collapse to one, because rayon
work-steals over a single pool: nesting parallel loops cannot explode the thread
count, so the C++'s `ThreadedAligner::setMultithreaded(false)` bookkeeping is
unnecessary.

The axis matters, though, and follows the C++ default (`process_with_pthreads`):

- **Phase 1** — one task per candidate reference, inner loop sequential. Each
  task also does its own MSA assembly and consensus call, so those parallelise
  too; with the loop the other way round they are a sequential barrier per
  candidate. Measured on 120 × 2.5 kb sequences, 64 cores: **15.0 s → 13.1 s**,
  byte-identical output.
- **Phase 2** — only one reference in flight, so the inner loop is parallel.

**Peak memory scales with thread count, not input size.** parasail's traceback
kernels allocate an `m × n` matrix per in-flight alignment; each worker holds one.
Measured:

| input | threads | wall | peak RSS |
|---|---|---|---|
| 120 × 2.5 kb | 1 | 73.3 s | 41 MB |
| | 4 | 20.0 s | 93 MB |
| | 8 | 10.7 s | 158 MB |
| | 16 | 10.9 s | 248 MB |
| 20 × 10 kb | 4 | 55.8 s | 1.56 GB |
| | 16 | 36.3 s | 5.40 GB |

Two things to take from that. Throughput **stops improving past ~8 threads** on
2.5 kb input — the work is memory-bandwidth-bound well before it is core-bound,
and more workers only add allocation churn. And at 10 kb the cost is roughly
**350 MB per worker**, quadratic in sequence length, so long families need the
thread count chosen deliberately. `autocons` estimates this at startup and warns
with a concrete `--threads` suggestion when it would exceed half of
`MemAvailable`.

If phase 1 ever becomes the bottleneck on large inputs, the fix is not more
threads: it is to replace its `O(n²)` full dynamic programming with a seeded
search — which is what `aln-rmblast` exists to provide.

### End-to-end validation

`autocons/tests/alu_recovery.rs` takes the real 311 bp AluY consensus, diverges
it into 20 copies at 15% substitution and 3% indel, buries each in random flanks
so hits are embedded rather than edge-aligned, and requires `autocons` to get it
back. It recovers AluY at **100% identity over 310 of 311 bases** — and still
clears 95% at 22% divergence.

That single test exercises the whole stack: FASTA I/O, parasail striped SIMD
alignment, MSA assembly under `GrowReference`, the Dfam consensus caller with CpG
restoration, and refinement to convergence.

Note the emitted consensus runs somewhat longer than the true one (361 vs 311 bp
above) because local alignment picks up flanking sequence at the edges. The C++
behaves the same way; trim as needed.

## Consensus calling

`aln_core::consensus` is the Dfam caller — Perl `MultAln.pm::
buildConsensusFromArray` by way of the verified Rust in `dfam-curator`. Public
names match `dfam_curator::consensus` so that migration is a delete-and-re-export.
It is the default in `acons`/`autocons`; the original GIRI caller is `--orig`
there and is not yet ported.

Two passes: per-column argmax over all 18 symbols (gap included, `N` preferred on
ties), then CpG restoration, which compares the called dinucleotide against a
hypothetical `CG` under a deamination model and overwrites both positions if `CG`
wins.

A caution for anyone writing tests here: because `C` and `G` carry the highest
self-scores in the matrix, **a balanced mix of `TG` and `CA` already calls `CG`
from the per-column pass alone** — such a dataset proves nothing about
restoration. `deaminated_cpg_is_restored` uses a skewed set where the plain call
gives `TG` and only the model recovers `CG`; `a_balanced_tg_ca_mix_needs_no_
restoration` pins the trap.

## Licensing

Code written for `dfam-lib` is **CC0-1.0** (`LICENSE`).

`aln-parasail` additionally redistributes a subset of **parasail** under
Battelle's BSD-3-Clause-plus-citation terms, kept in
`aln-parasail/third_party/parasail/` alongside its unmodified `LICENSE.md`;
the crate declares `license = "BSD-3-Clause AND CC0-1.0"`. `THIRD-PARTY.md` is
the index of everything vendored and what each licence actually requires — read
it before adding a second upstream, because the layout rule (one directory, one
upstream, one licence) is what keeps attribution mechanical. Anything published
from parasail-backed output should cite Daily (2016).

Vendoring parasail is also what makes a C compiler a build requirement, which is
why the intended arrangement is for consumers to gate it behind a Cargo feature
rather than depend on it unconditionally:

```toml
[features]
default = []                     # a library should not imply a C toolchain
parasail = ["dep:aln-parasail"]

[dependencies]
aln-parasail = { workspace = true, optional = true }
```

with a single `#[cfg(feature = "parasail")]` in one backend-selection function
returning a `Box<dyn DynAligner>`, and `aln-reference` as the fallback. Binaries
can flip `default` on; there is no consumer crate in the workspace yet, so the
gate is not wired anywhere today. `aln-parasail` vendors only the x86 kernels —
on Apple Silicon `build.rs` fails with a message naming the feature to turn off.

## Not yet built

- `aln-farrar` — FFI over `libalign/src/swsse2`. It is the default `acons`
  aligner (`autocons.cpp:475-481`), so any bit-exactness claim against current
  output is a claim about Farrar. **Blocked on licensing**: those sources carry
  "Copyright 2006 by Michael Farrar. All rights reserved. This program may not be
  sold or incorporated into a commercial product, in whole or in part, without
  written consent" — incompatible with redistributing this workspace under a
  permissive licence. parasail's striped kernels implement the same algorithm
  under BSD-3-Clause and are already wired up; they are algorithmically
  equivalent but will not be bit-identical to GIRI's modified variant.
- `ProcessRepeats`' adaptive-width `.out` table (a two-pass column-width
  calculation over a whole file, not a per-record writer).
- `dfam-curator` migration onto `aln-core` (planned for a later session).
- **Open, deferred:** the `rmblast-lib` left-extension underflow above — decide
  whether to patch `rmblast-port` and drop the two `#[ignore]` attributes.
