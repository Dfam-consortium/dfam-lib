# dfam-lib

A Rust library to support tools used in the curation of transposable element
sequences.  This library includes interfaces to pairwise alignment programs 
(parasail/rmblast etc), data structures for storing alignment data, file parsers
/exporters, and general analysis functions.  Aside from a handful of new
features, this library includes many existing concepts from the RepeatMasker, 
RepeatModeler, Dfam and GIRI/Repbase projects.

```
dfam-lib/
├── aln-core/        sequences, matrices, alignments, divergence, MSA,
│                    consensus calling, output formats, FASTA / 2bit /
│                    crossmatch I/O
├── aln-engine/      the two traits + parallel drivers
├── aln-parasail/    parasail 2.6.2 SIMD backend (vendored subset)
├── aln-reference/   plain O(mn) Gotoh implementation for comparison
├── aln-rmblastn/    SearchEngine over the external rmblastn binary
├── aln-rmblast/     SearchEngine over the rmblast Rust port (git dependency)
├── dfam-stk-io/     Stockholm 1.0, and the bridge to MultiAlign
├── docs/            notes too long to sit in a doc comment
└── tools/           shell probes for the C++ binaries this stack replaces
```

## Building

```sh
git clone https://github.com/Dfam-consortium/dfam-lib
cd dfam-lib
cargo test --workspace
```

What that needs on the machine:

- **Rust 1.86 or newer.** `rust-version` in the root manifest is the oldest
  toolchain the workspace is tested on; cargo refuses older ones with a clear
  message rather than a page of errors.
- **A C compiler.** `aln-parasail` compiles the vendored parasail kernels
  with `-msse2`, `-msse4.1` and `-mavx2` and picks one at run time by `cpuid`,
  so the compiler has to accept those flags (gcc and clang on x86-64 do).
- **A C++ compiler.** `aln-rmblast` depends on `rmblast-lib`, whose default
  `alp-fit` feature compiles the ALP sources vendored in the RMBlast
  repository.
- **Network access on the first build.** Two dependencies come from GitHub by
  tag rather than from crates.io: `rmblast-lib` (RMBlast `3.0.6`) and
  `smitten` (`1.0.2`). Cargo caches them; later builds are offline. This is
  also why the crates here are not on crates.io: a published crate cannot
  depend on a git repository.
- **Optionally, an `rmblastn` binary on `PATH`.** `aln-rmblastn`'s live tests
  run against it and skip with a message when it is absent. Nothing else
  needs it.

Nothing here has to be installed; consumers depend on the crates by git tag
(see [Consumers](#consumers)).

## The workspace

The root `Cargo.toml` has a `[workspace]` table and no `[package]`: a *virtual*
workspace, a directory that owns member crates but builds nothing itself. The
members share one `target/`, one `Cargo.lock` and one resolved dependency graph,
so `cargo test` at the root compiles `aln-core` once however many members depend
on it, and a dependency declared once resolves once for all of them.

Settings that must not drift live in `[workspace.package]` and
`[workspace.dependencies]`; members write `version.workspace = true` and
`thiserror.workspace = true` rather than repeating a value someone then has to
keep in sync by hand. Path entries in `[workspace.dependencies]` are what wire
the members to each other.

Useful flags: `-p aln-core` scopes a command to one member and `--workspace`
widens it to all of them, while `--all-targets` reaches the examples, benches and
test targets a bare `cargo check` leaves alone.

`aln-rmblast` is the only member with a dependency outside crates.io: it takes
[RMBlast](https://github.com/Dfam-consortium/RMBlast), the Rust port of the
search engine, as a git dependency pinned to tag `3.0.6`. A build therefore
fetches from GitHub the first time, and `rmblast-lib`'s default `alp-fit`
feature compiles the ALP sources vendored in that repository, which needs a C++
compiler. `aln-parasail` already compiles vendored C, so the C++ compiler is the
only new build requirement.

Bump the tag in `aln-rmblast/Cargo.toml` when RMBlast releases, then run
`cargo update -p rmblast-lib` to move the lockfile with it.

To develop against a local `rmblast` checkout, swap the `git`/`tag`
dependency for the commented-out `path` line beneath it. A `[patch]` section
will not serve here: cargo resolves the git source before applying the patch, so
it still fetches the tag you were trying to bypass.

## The three conventions

Everything else follows from these. Each is documented at length in its module.

| | rule | why it bites |
|---|---|---|
| **Coordinates** | 0-based, half-open, forward-strand, `start <= end`; a range is an `aln_coord::Span`, an absent range is `Option<Span>` | The file formats around us (Smitten identifiers, Stockholm, RepeatMasker `.out`, BLAST tabular) are 1-based closed, and parasail's accessors are 0-based closed. A bare `u64` with the convention in a comment survives every test whose fixtures are self-consistent; `Span` makes the conversion a named call. |
| **Gap vs. padding** | `-` means present-but-deleted; `' '` means not present | GIRI uses `<`/`>` for padding instead; `aln_core::seq` converts. Counting padding as a gap silently corrupts coverage and divergence denominators. |
| **Matrix orientation** | crossmatch's: `matrix[subject][query]`, rows = subject = consensus, columns = query = genomic. NCBI's files are the transpose | Symmetry is a property of each matrix, not of the format. Arian Smit's `##p##g` matrices are asymmetric (`14p35g`: `G/A` −7, `A/G` −10), and so is the consensus caller's. A swapped lookup on those changes scores without failing. |

### Coordinates in more detail

The default is the unmarked one. A field called `start`, `end` or `span`, with
no suffix, is 0-based half-open on the forward strand. Strand lives in its own
field (`Alignment::strand`, `SequenceRow::orient`); a reverse-strand interval
is still stored ascending.

Anything else says so in its name. `Alignment::query_one_based()` and
`Span::as_1b_closed()` return 1-based closed pairs for writers;
`Span::from_1b_closed()` is the constructor parsers use. dfam-curator's
`dfam-coord` audits identifiers in their published form and keeps
`start_1b`/`end_1b` for that reason. The only fields that are neither are
positions rather than ranges: RepeatAfterMe's DP engine walks single indices
(`left_seq_pos`, `upper_seq_bound`) and reads them back as ranges through
`CoreAlignment::core_span()`.

`SequenceRow::col_start`/`col_end` look like a third convention. They are
column indices into the gapped row, half-open like everything else. Rust
slicing is half-open too, which is why the house convention is:
`&seq[span.range_usize()]` needs no arithmetic.

`docs/coordinate-migration.md` records how the crates got here and what each
consumer had to change.

### Matrix orientation in more detail

Two orientations are in use, and RepeatMasker ships every matrix in both:

| | rows | columns | files |
|---|---|---|---|
| crossmatch | subject (consensus) | query (genomic) | `Matrices/crossmatch/*.matrix` |
| NCBI | query | subject | `Matrices/ncbi/nt/*.matrix` |

`SubstMatrix::score(subject, query)` uses crossmatch's, which is also
`SearchResult.pm::rescoreAlignment`'s. `aln_rmblast::matrix` transposes on the
way out to rmblast, and `aln-parasail` transposes on the way in, because
parasail looks up `matrix[s2][s1]`. `aln-rmblastn` does not: it takes a path
to an NCBI-format file and cannot check which layout it was handed.

Whether the transpose matters depends on the matrix:

- The `##p##g` series is asymmetric throughout. `14p35g` scores `G/A` as −7
  and `A/G` as −10.
- RepeatModeler's `comparison.matrix` is symmetric over `A C G T` and
  asymmetric over its IUPAC rows: `C/Y` is 2, `Y/C` is 1.
- The consensus-calling matrix (`aln_core::consensus::MATRIX`, from GIRI's
  `DNACONMATRIX`) is asymmetric by design. Rows are the candidate base,
  columns the observed one: candidate `A` against observed `G` is −8, the
  reverse −4.

So a swapped lookup is a silent error for all three, and only looks harmless
on `comparison.matrix` while no ambiguity codes are present.

## Aligner Traits

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

### Memory scales with thread count, not input size

The traceback kernels allocate an `m × n` matrix per in-flight alignment, and
each worker holds one. Measured through `autocons`, which is the heaviest
consumer:

| input | threads | wall | peak RSS |
|---|---|---|---|
| 120 × 2.5 kb | 1 | 73.3 s | 41 MB |
| | 4 | 20.0 s | 93 MB |
| | 8 | 10.7 s | 158 MB |
| | 16 | 10.9 s | 248 MB |
| 20 × 10 kb | 4 | 55.8 s | 1.56 GB |
| | 16 | 36.3 s | 5.40 GB |

Two things follow. Throughput **stops improving past ~8 threads** on 2.5 kb
input — the work is memory-bandwidth-bound well before it is core-bound, and
more workers only add allocation churn. And at 10 kb the cost is roughly
**350 MB per worker**, quadratic in sequence length, so long families need the
thread count chosen deliberately rather than left at the core count. A caller
holding `O(n²)` alignments should estimate this at startup; `autocons` does, and
warns with a concrete `--threads` suggestion when the estimate would exceed half
of `MemAvailable`.

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

### Fixed upstream: the left extension underflowed

Through RMBlast 3.0.4, `rmblast-lib`'s left-extension (`REVERSE`) pass in
`search/gapped.rs` formed its `b` pointer as `b.as_ptr().add(n - 1 -
first_b_index)` ahead of the loop bound that made it safe. When `first_b_index`
reached `n` the `usize` subtraction underflowed: a panic under
`debug_assertions`, and in release an out-of-range pointer that nothing ever
dereferenced, so results stayed correct while forming the pointer was undefined
behaviour. Any hit that did not start at the beginning of the subject triggered
it, and five bases of left flank sufficed, which is the ordinary RepeatMasker
shape.

RMBlast 3.0.5 clamps the pointer when the loop will not run.
`embedded_hit_lands_at_the_right_offset` and
`minus_strand_keeps_plus_strand_subject_coordinates` carried
`#[cfg_attr(debug_assertions, ignore)]` for this and now run in every profile.

### Searching in batches

Driving a search engine one pair at a time rebuilds a query lookup table per
call and scans a single subject, which throws away the point of seeding.
`RmblastEngine::one_to_many` and `::all_vs_all` issue **one** search and
demultiplex the results — the two shapes `autocons` needs for refinement and for
reference selection.

Unlike a pairwise aligner, rmblast can return several HSPs for one
`(query, subject)` pair: an instance interrupted by an insertion, or one
matching on both strands. `mask_level` is the filter — an HSP is dropped when a
higher-scoring HSP already covers more than that percentage of its **query**
span. RepeatMasker's default is 80, and 101 disables it. Survivors are all
returned, so one instance can contribute more than one row to an MSA, which is
GIRI's own `FRAGMENT` model and needs nothing downstream to change.

**Masking is applied per pair here, not across subjects.** rmblast's own
`apply_mask_level` follows NCBI's `Blast_HSPResultsApplyMasklevel` and ranks a
query's HSPs across *all* subjects together, so a hit against one subject can
suppress a hit against another. That is right for genome annotation, where you
want one family per region, and wrong here: in an all-against-all search every
sequence hits itself perfectly over 100% of its own query span, and cross-subject
masking then discards essentially every real cross-hit. Measured: with
cross-subject masking at 80, 29 of 30 families produced no consensus at all. So
these two entry points search with masking disabled and apply the same rule per
pair afterwards.

## aln-rmblastn

The same trait over the **external `rmblastn` binary** — what RepeatMasker
actually runs. It spawns a process and writes temp files, so it is slower on
small inputs and cannot be unit-tested, but its output is the reference for
anything that has to agree with published annotation, and it is the only path
that can search a prepared BLAST database. Choosing between it and the
in-process port is a deployment decision, not a correctness one, which is why
both sit behind `SearchEngine`.

Flags follow `NCBIBlastSearchEngine.pm`, so the gap-cost conversion above
applies here too: crossmatch's `-25`/`-5` goes out as `-gapopen 20 -gapextend 5`.
X-drop is not an independent knob — NCBI derives all three cutoffs from the
score floor (`min_score × 2`, `÷ 2`, `× 1`), so lowering `min_score` also
shortens extensions.

Two things it refuses rather than guesses:

- **It will not synthesise a matrix.** `SearchParams::matrix` holds a *parsed*
  matrix in crossmatch layout and NCBI's files are the transpose, so
  `RmblastnOptions::matrix_path` must point at a real NCBI-format file
  (RepeatMasker's `Matrices/ncbi/nt/`). Constructing the engine with one set and
  not the other is an error rather than a fallback. rmblastn also requires a
  bare filename, resolved through `BLASTMAT`, which the engine sets from that
  file's parent directory.
- **`makeblastdb` comes from `rmblastn`'s own directory**, not from `PATH`. One
  machine can carry several BLAST+ installs (here `/usr/local/rmblast/bin` is
  2.17.1 while `/usr/bin` is 2.12.0), and indexing with one version then
  searching with another corrupts silently instead of failing.

`tests/live.rs` covers what the unit tests cannot: flag construction, the
`makeblastdb` step, and whether the tabular output asked for is the output that
comes back. It skips itself when no binary is installed.

## dfam-stk-io

Stockholm 1.0 — Dfam's seed-alignment format — together with the conversion to
and from `aln_core::msa::MultiAlign`. Both halves are in one crate so a second
tool cannot reimplement the bridge and drift away from this one.

`StkRecord` is the format itself: `#=GF` file annotation, the `#=GC RF`
reference line, sequences interleaved across multiple blocks, `//` terminators,
and a streaming `iter_records` so a multi-record file need not be held in memory
at once. `msa::read_select` pulls a single record by 1-based number or by its
`#=GF ID`.

Identifiers go through [Smitten](https://github.com/Dfam-consortium/Smitten)
rather than a local regex, so recursive ranges normalise the way they do
everywhere else in the toolchain; rows whose names do not parse, bare consensus
labels among them, are still stored without coordinates. Smitten is
pinned to tag `1.0.2`: an untagged git dependency resolves to whatever the
default branch holds at fetch time, so two clones of the same dfam-lib tag would
otherwise build different code.

Four gap characters are read — `-`, `.`, `_`, `~` — and `.` is written.

## Reading what other tools wrote

Two readers sit in `aln-core` because the alternative was a third and a fourth
copy of each.

**`twobit`** — random-access UCSC `.2bit`. Consolidated from three
implementations (dfam-curator, RepeatAfterMe, this crate); the one kept is
RepeatAfterMe's, which rejects inverted ranges and short-circuits empty ones
where the others underflow `end - start` on `u64` and panic in
`Vec::with_capacity`. The index is built on open, and `fetch` reads only the
packed bytes covering the request via `pread`, so the reader is `Send + Sync`
and `fetch` takes `&self`. It returns uppercase `ACGTN`; mask blocks are
ignored, matching the C loader's unconditional `toUpperN`. Both endiannesses,
versions 0 and 1, Unix only.

**`crossmatch`** — a reader for the `.align` files `fmt::to_crossmatch` has been
able to write all along. The parser lived in `dfam-curator`, so anything else
that wanted to read RepeatMasker alignment output had to depend on a curation
tool to get at it.

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

## Consensus calling

`aln_core::consensus` is the Dfam caller — Perl `MultAln.pm::
buildConsensusFromArray` by way of the verified Rust in `dfam-curator`. Public
names match `dfam_curator::consensus` so that migration is a delete-and-re-export.
It is the default in `acons`/`autocons`; the original GIRI caller is `--orig`
there, and is described below.

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

### The GIRI caller

`aln_core::giri` is the other one: `MultipleAlignment::getConsensus()` from
`giri_cpp_lib`, by way of the verified Rust in `dfam-curator`. It is a different
caller rather than a variant of the same one:

- **12 symbols** (`A R G C Y T K M S W N X`) against the Dfam caller's 18 — it
  cannot emit `B/D/H/V/Z`.
- **Gap weights** `-5` base-vs-gap and `+2` gap-vs-gap, against `-6`/`+3`, and
  an `N` penalty of `-5` against a softened `-2`.
- **A minimum-coverage gate** (`acons --min`) with no Dfam equivalent.
- **Ties resolve the other way.** Candidates are scanned in reverse with a
  strict `>`, so gap beats `X` beats `N` beats … beats `A`. The Dfam caller
  prefers `N`.
- **CpG restoration is a separate pass**, applied only when a species is given
  (`--mam`), under a probabilistic model rather than the Dfam caller's
  deterministic bonus.

`giri::fixed` holds a corrected restoration: GIRI's `firstBase` lookup table is
inert, and its `CG` mask neither skips gaps nor stops at the motif end. It is
kept beside the faithful version rather than folded into it, because reproducing
`acons` is what the faithful version is for.

## Consumers

This repository is libraries only. The tools live in their own repositories and
depend on it:

- **`dfam-curator`** — `autocons`, `cons-core` (the consensus pipeline `autocons`
  and `te-composer` share), and the curation tooling around them.
- **`RepeatAfterMe`** — `ram-core`, which contributed the 2bit reader now in
  `aln-core`, and `ram-cli`.

Both depend on dfam-lib **by git tag**, never by path:

```toml
[workspace.dependencies]
aln-core = { git = "https://github.com/Dfam-consortium/dfam-lib", tag = "0.0.3" }
```

A path dependency ties a consumer to one directory layout on one machine: when
dfam-lib was relocated, every crate in dfam-curator stopped building until
someone edited the paths. A tag resolves the same way from any clone.

For work against a live checkout, a commented-out `[patch]` at the foot of the
consumer's root manifest substitutes the working copy for the tag without any
committed manifest changing:

```toml
#[patch."https://github.com/Dfam-consortium/dfam-lib"]
#aln-core    = { path = "../dfam-lib/aln-core" }
#dfam-stk-io = { path = "../dfam-lib/dfam-stk-io" }
```

`[patch]` rewrites the *source*, not one dependency edge, so it applies to the
whole graph at once: transitive users of dfam-lib follow the working copy too,
and the build stays on a single copy of `aln-core`. That property is why they
must not path-depend on dfam-lib themselves: a path dependency escapes the
patch, and two `aln-core`s in one graph give you types that will not
interconvert. Uncomment to work, comment out again and let `Cargo.lock`
re-resolve before committing.

### autocons

A faithful port of `acons/src/autocons.cpp`, and the heaviest user of this
stack. It lives in `dfam-curator` now, so what is left here is what it asks of
the library.

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

#### Parallelism, and how to size a run

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

If phase 1 ever becomes the bottleneck on large inputs, the fix is not more
threads: it is to replace its `O(n²)` full dynamic programming with a seeded
search — which is what `aln-rmblast` exists to provide.

#### End-to-end validation

`cons-core/tests/alu_recovery.rs`, in dfam-curator, takes the real 311 bp AluY
consensus, diverges it into 20 copies at 15% substitution and 3% indel, buries
each in random flanks so hits are embedded rather than edge-aligned, and requires
`autocons` to get it back. It recovers AluY at **100% identity over 310 of 311
bases** — and still clears 95% at 22% divergence.

That single test exercises the whole stack: FASTA I/O, parasail striped SIMD
alignment, MSA assembly under `GrowPerSlot`, the Dfam consensus caller with CpG
restoration, and refinement to convergence.

Note the emitted consensus runs somewhat longer than the true one (361 vs 311 bp
above) because local alignment picks up flanking sequence at the edges. The C++
behaves the same way; trim as needed.

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
- The rest of the `dfam-curator` migration. Its consensus caller, MSA types,
  Stockholm reader and 2bit reader now come from here — `io/twobit.rs` is a
  re-export and nothing else. What is left is `blast.rs`, a second rmblastn
  wrapper that predates `aln-rmblastn` and does the same work behind a different
  interface.
