# Standardising coordinates on 0-based half-open

Working notes for a migration in progress. Phase 0 shipped in `0.1.1`; the
current tag is `0.1.2`, which only fixed a version field. Phases 1–5 shipped in
`0.2.0`. This records the reasoning and the dead ends, not just the steps,
because the next person will not have been in the room.

## Why

Four conventions are in use across the Dfam and GIRI lineages:

| convention | who uses it | `chr1` bases 101..200 |
|---|---|---|
| 0-based half-open | `rmblast-lib`, `aln_core::Alignment`, `aln-coord` | `100, 200` |
| 1-based fully closed | Smitten IDs, RepeatMasker `.out`, BLAST tabular, Stockholm | `101, 200` |
| 0-based fully closed | parasail's `end_query`/`end_ref`; RepeatAfterMe `glocal`/`library` until Phase 5 | `100, 199` |
| 1-based half-open | nobody | n/a |

RepeatAfterMe carries all four in 3617 lines: `ram-core/src/library.rs:17`
(half-open), `library.rs:56,527` (0-based closed), `glocal.rs:26-27` (0-based
inclusive), `ram-cli/src/main.rs:553` (1-based closed on output).

Dfam already pays for this. `dfam-curator/dfam-coord` is a 1082-line crate whose
`validate_sequences` (`src/lib.rs:424`) **guesses which convention an identifier
was written in** by loading the genome and testing whether the sequence
matches:

```rust
if range_length == fasta_sequence_length {
    // "Detected half-open coordinates ... Converting to one-based fully closed."
    validation_str.push_str("_halfopen");
    record.start = record.start.map(|start| start + 1);
}
```

and brute-forces shifts of ±1, ±2, ±3 when that fails. It has to guess because
the convention lives in a doc comment attached to a `pub u64`.

## Why half-open and not 0-based closed

0-based closed was considered and rejected. Four reasons:

1. **Empty ranges become unrepresentable.** `dfam-stk-io/src/msa.rs:104` returns
   `(0, 0)` for an unparseable name; under 0-based closed that reads as base 0.
2. **`Alignment::new` would underflow.** `align.rs` computes
   `query_end = query_start + edits.query_consumed()`. Closed makes it
   `+ consumed - 1`, which panics on an empty script.
3. **Rust slicing is half-open.** Closed coordinates move the arithmetic out of
   a handful of tested parser boundaries and into every untested slice site.
4. **Adjacency.** `a.end == b.start` means "abuts". RepeatAfterMe's
   `boundaries[i]` (`library.rs:84`) is already documented exclusive; closed
   bounds would need a `+1` at every concatenation seam.

parasail is not an argument for closed: `aln-parasail` never calls
`parasail_result_get_end_query`/`_end_ref` (declared at `aln-parasail/src/ffi.rs:126-128`,
zero call sites). It reads `beg_query`/`beg_ref` from the CIGAR and derives ends
from the consumed counts. `rmblast-lib`'s `Hsp`, the other engine, is already
0-based half-open.

Counted by grep: migrating **to** 0-based closed touches ~185 sites in dfam-lib;
migrating to half-open touches ~75, all in parsers and one struct.

## The two rules

**Never change what a field means under its existing name.** Rename in the same
commit so stale callers fail to compile. Renaming is the substitute for a
compile-time check, and it is what makes each phase safe. A type change
(`u64` → `Option<u64>`, `u64` → `Span`) is self-enforcing and needs no rename.

**Annotate the deviation, not the default.** Bare `start`/`end` mean 0-based
half-open and stay unmarked. Anything else carries a suffix (`_1b`, `_incl`).
Most code then needs no change at all.

## What shipped (Phase 0, code in `0.1.1`)

**`aln-coord`**, a new leaf crate, zero dependencies, holding `Span`: a 0-based
half-open forward-strand interval with private fields.

```
Span::new(start, end)            0-based half-open, the house convention
Span::from_1b_closed(start, end) Smitten, .out, .align, Stockholm, BLAST tabular
Span::start() / end()            0-based
Span::as_0b_half_open()          -> (u64, u64)
Span::as_1b_closed()             -> Option<(u64, u64)>, None when empty
Span::range_usize()              -> Range<usize>, for &seq[span.range_usize()]
Span::len() / is_empty()
aln_coord::calibration::{Case, check_all}
```

Nine methods, deliberately. The `# Scope` section in the module docs records why
each absence is deliberate; read it before adding anything. Adding a method later
is a compatible change; removing one breaks callers, and the repair has to move
three repos together.

It is a separate crate rather than a module in `aln-core` so that anything
needing only coordinates does not compile 7,000 lines of alignment machinery, and
so Smitten could take it later. Smitten ships Perl and Python implementations of
the same grammar, so it can only depend on a leaf.

**`Alignment::validate` rejects zero-length alignments** (`aln-core/src/align.rs`).
An empty edit script, or one consuming nothing on one side. Nothing in the
workspace produced either; the full suite passed unchanged. This is what makes
`Span::as_1b_closed()`'s `None` case unreachable for a validated `Alignment`, so
writers can `.expect()` it away.

**`crossmatch::PairwiseHit` holds `query: Span` and `subj: Span`**, the first
real boundary converted, and the pilot for the rest.

### What that pilot caught

The rename broke two tests, which failed with:

```
subject coordinates in "254 28.0 4.0 9.0 seq-13 3873751 3873963 (114) C L2d 3331 3075 (6)":
1-based coordinates 0-3075 include 0
```

The fixture was hand-written with `(sLeft)` in the forward-strand position. Real
RepeatMasker C lines put it before the coordinates:

```
243 30.00 0.00 0.00 seq 970 1069 (22196) C L2-1_AMi#LINE/L2 (679) 556 457 m_b1s551i0 3
```

The parser was right; the fixture was wrong, so `subj_start` had been parsing as
0 for as long as the test existed. No test asserted on that field.

`crossmatch::span_1b` now also rejects a header covering no bases, for the same
reason `Alignment::validate` does.

**dfam-curator** consumed it in `src/build.rs`: four accessors became
`ref_span`/`inst_span` plus a `one_based` converter. No call sites changed.

## Phases 1–4 (shipped in `0.2.0`)

The plan below was written as five phases, to be shipped one tag at a time.
They landed together instead, in one change across the three repos. The reason: every phase was a change of both name and type
(`seq_start: u64` → `span: Option<Span>`), so the compiler found every site
in one pass, and the intermediate state the phases were protecting against
(`Option<u64>` still holding 1-based values under a type that does not say
so) never had to exist. The rules held: nothing changed meaning under its old
name, and only the deviations carry a suffix.

### What changed

**`aln_core::msa::SequenceRow`**

| before | after |
|---|---|
| `start: usize`, 0-based column, inclusive end | `col_start: usize` |
| `end: usize` | `col_end: usize`, exclusive; `(0, 0)` for an all-padding row |
| `seq_start: u64`, `seq_end: u64`, 1-based closed, `(0, 0)` = none | `span: Option<Span>` |

`MsaMember` lost `seq_start`/`seq_end` for `span: Option<Span>` the same way.

**`dfam_stk_io::SeqRow`**: `seq_start`/`seq_end: Option<u64>` → `span:
Option<Span>`, converted from Smitten's 1-based `Range` at parse time. The
comment that called the old fields "0-based" was wrong, as the plan noted.

**`dfam_stk_io::msa::parse_seq_name_coords`** is now `pub` and returns
`(String, Option<Span>, Strand)`. dfam-curator's FASTA and clustal readers each
carried a copy of it; both now call this one. A name whose coordinates include
a 0 (`chr1:0-100`) is no longer split: 0 cannot occur in a 1-based name, so
the whole string is kept as the name with no span. Previously it parsed as
`seq_start = 0`.

**dfam-curator**

- `blast::BlastHit`: `query_start`/`query_end`, `subj_start`/`subj_end` →
  `query: Span`, `subj: Span`, built with `from_1b_closed` in `parse_hits`, so
  a malformed BLAST line is an error at parse time instead of a bad row later.
- `build.rs`: `one_based` and its four wrappers are gone. `ref_min`/`ref_max`
  are 0-based half-open, and the reference row now gets a full span (it used
  to get `seq_start = ref_min` and `seq_end = 0`).
- `dfam-coord::SequenceRecord`: `start`/`end` → `start_1b`/`end_1b`. This
  crate audits identifiers as written, so it keeps the 1-based closed values
  and says so in the name. The `_halfopen` repair logic is untouched.
- `cons-core/src/lowqual.rs`: seven column comparisons moved from inclusive to
  exclusive end. An all-padding row used to look like it covered column 0.
- `te-composer/src/extend.rs`: `RowSpan` holds `span: Option<Span>`; a copy
  with a coordinate-less row is reported as such rather than as
  "extent 0-0 is not a valid range".
- `linup` and `linup_fmt` print `0-0` for a row with no coordinates, which is
  what Perl printed. Nothing stores that value; it is only what the report
  prints.
- `tsd`: skips rows with no span or an empty one; `lb`/`rb` come straight from
  `span.start()`/`span.end()`, and the `- 1` is gone.

**RepeatAfterMe** `ram-cli/src/stk.rs`: the `s0 - 1` conversion is gone, since
`SeqRow.span` is already half-open. Its calibration test passes unchanged.

### Verification done

- dfam-lib, dfam-curator (244 tests), RepeatAfterMe: all tests pass.
- New calibration tests: `dfam_stk_io::msa::tests::name_coordinates_are_converted_to_half_open`
  (`chr1:101-200_+` → `100..200`; `gi|57:120437225-120436960` → minus strand),
  `written_labels_are_one_based_closed` (write-side),
  `dfam_coord::coordinate_tests::record_coordinates_are_as_written_in_the_identifier`
  (round trip back to 1-based).
- `linup --format {linup,stockholm,stats,consensus}` over six real Stockholm
  files (`ex1`, `foo`, `timea`, `r1f2`, `iRic2.1.1469`, `AnERV11b_Rmuscosa`):
  output byte-identical to the pre-migration binary.
- Not done: the `validate_sequences` corpus run, which needs genomes on disk.

### Build state

Released. dfam-lib `0.2.0` (`34f2134`), RepeatAfterMe `RepeatAfterMe_V0.2.0`
(`8fbc350`), dfam-curator `8c8fa56` on `main` pinning both. The `[patch]`
sections in the two consumers are commented out again.

## Phase 5: RepeatAfterMe (shipped in `RepeatAfterMe_V0.2.0`)

Done in the same release. The engine (`engine.rs`) is a banded DP that
walks single positions in the concatenated library buffer with the C's
wrapping `uint64_t` arithmetic; `left_seq_pos`, `right_seq_pos`,
`lower_seq_bound` and `upper_seq_bound` on `CoreAlignment` are positions, not
ranges, and a position has no open/closed convention to annotate. They stay as
they are, with their docs rewritten to say so. Every *range* now crosses the
crate boundary as a `Span`:

| type | before | after |
|---|---|---|
| `library::RangeRecord` | `start: i64`, `end: i64` (file is half-open) | `span: Span`; `read_ranges` rejects a negative or inverted pair |
| `library::CoreEdge` | `start`/`end: u64`, 0-based closed | `span: Span` on the source sequence |
| `glocal::GlocalResult` | `query_start`/`query_end: i32`, `subj_start`/`subj_end: u64`, 0-based closed | `query: Span`, `subj: Span` |

Two accessors replace the hand arithmetic: `CoreAlignment::core_span()` and
`extended_span()` (library coordinates), and `SequenceLibrary::to_source(idx,
span)`. `ram-cli/src/main.rs`'s report, TSV and FASTA writers now go
`to_source(...).as_1b_closed()` instead of repeating `pos - seq_lower + 1 +
subseq_offset` in each writer. The C-compat TSV anchor and the FASTA
flank clamps still work on positions, because the C quirks they reproduce
(clamping to 1, reading the boundary base) are stated on positions.

`ram_core::Span` re-exports `aln_coord::Span` so ram-cli and te-composer can
name it without a second dependency.

te-composer's `Anchor` (`extend.rs`) holds `span: Span` instead of
`start`/`end: i64`; its core-edges report keeps printing 0-based closed
because it exists to be diffed against the C tool's `printCoreEdges`.

### Verification done

- `harness/diff-c-rust.sh`: all 12 C-vs-Rust comparisons pass (three ce10
  families, four parameter sets), before and after.
- The harness does not compare the stdout report, so that was diffed
  separately against the pre-change `ram-extend` from
  `~/projects/Claude/RepeatAfterMe/target/release`: stdout, TSV and FASTA
  byte-identical on all three families under two parameter sets.
- `glocal`'s C golden values are restated as spans in the test
  (`67..=267` → `Span::new(67, 268)`), which is the calibration for that
  boundary.

## Smitten does not change

Its `Range` is a faithful decode of a string that is 1-based fully closed by
published spec. Converting at the `dfam-stk-io` boundary drops 29 sites and the
one cross-repo edit with two consumers and no compile-time net. The identifier
string format never changes.

## Verification

**Calibration tests.** One per boundary, with literal numbers worked out by hand,
using `aln_coord::calibration::{Case, check_all}`. `aln-coord/tests/calibration.rs`
has four worked examples. The pattern comes from RepeatAfterMe's
`coordinates_match_stk2ranges_convention` (`ram-cli/src/stk.rs:279`), which
asserts `chr1:101-200_+` → `(100, 200)`.

**Acceptance test you already own.** Run `dfam-coord`'s `validate_sequences` over
a corpus before and after each phase. The `fixed_halfopen` count must not
increase. It loads the genome and checks the sequence actually matches, so it
catches what unit tests cannot.

**Layout.** `aln-coord/tests/layout.rs` pins `size_of::<Span>() == 2 *
size_of::<u64>()`. `Alignment` is held `O(n^2)` live during all-against-all, so a
widened `Span` would be a real cost. Measured: the coordinate block is 64 bytes
with or without `Span`, `Option<Span>` (24) is narrower than two `Option<u64>`
(32), and LLVM folds a `Span` loop and a raw tuple loop to the same address.

## Release procedure

Every repo pins dfam-lib by git tag, so two tags are two sources and cargo
compiles two copies of `aln-core`. Ship in this order:

1. **dfam-lib**: bump `[workspace.package] version` **and** move the tag in the
   same commit. All eight crates use `version.workspace = true`; check none has
   regressed to a hardcoded version.
2. **RepeatAfterMe**: bump the three pins (`ram-core/Cargo.toml:11`,
   `ram-cli/Cargo.toml:19-20`), rebuild so the **tracked** `Cargo.lock` records
   the new rev, then tag.
3. **dfam-curator**: bump the seven dfam-lib pins and the `ram-core` pin in one
   commit. Never partially: dfam-lib at the new tag with `ram-core` at the old
   one gives

   ```
   error[E0308]: expected `ram_core::twobit::TwoBitReader`,
                    found `aln_core::twobit::TwoBitReader`
   note: two different versions of crate `aln_core` are being used
   ```

   at `te-composer/src/extend.rs:558`, because `ram-core/src/twobit.rs:8` is
   `pub use aln_core::twobit::*`.

dfam-curator is the leaf and needs **no tag**. Its `.github/workflows/release.yml`
fires on `push: tags: 'v[0-9]+.[0-9]+.[0-9]+*'`, so tagging it triggers a rebuild
nobody asked for. Commit and push to `main`.

Check the end state with `cargo tree -i aln-core`: exactly one entry.

**After moving a tag**, run `cargo update -p aln-core` anywhere already resolved.
`Cargo.lock` pins the rev, not the tag, so a stale lock keeps building the old
commit without saying so.

To develop across repos before tagging, uncomment the `[patch]` block at the foot
of `dfam-curator/Cargo.toml`. It rewrites the git source for the whole graph, so
RepeatAfterMe follows the working copy too and everything stays on one
`aln-core`. Comment it out before committing.

## Local checkouts

The working clones for this migration are the three side by side under
`~/projects/Claude/combined/`: `dfam-lib`, `dfam-curator`, `RepeatAfterMe`
(with `Smitten` and `RMBlast` alongside for reference). The `[patch]` paths
assume that layout. Older clones of the same repos sit one level up
(`~/projects/Claude/RepeatAfterMe`, `~/projects/Claude/dfam-curator`,
`~/projects/Claude/dfam-lib-project/dfam-lib`); they are behind. Duplicates of the same repos at the
same commit exist elsewhere (`~/projects/RepeatAfterMe`,
`/u3/local/src/dfam-curator`, and various `attic`/backup copies). Edit the wrong one and it builds and tests clean, and then `git status` in the
real one shows nothing. `/u3/home/rhubley` and `/home/rhubley` reach the same files.

`/u3/local/src/dfam-curator` is stale (May 2026) and still contains the removed
`src/alignment.rs`; do not read it for current structure.

## State at handoff

| repo | state |
|---|---|
| dfam-lib | tag `0.2.0`, pushed |
| RepeatAfterMe | tag `RepeatAfterMe_V0.2.0`, pushed |
| dfam-curator | `main` at `8c8fa56`, pushed, pins `0.2.0` / `RepeatAfterMe_V0.2.0` |

All five phases are shipped. Open: `#[non_exhaustive]` on the public structs
before anything outside these three repos depends on them.
