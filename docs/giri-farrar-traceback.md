# How GIRI gets a traceback out of a striped SIMD aligner

Notes from reading `giri_cpp_lib/src/libswaligner/{FarrarSWAligner,SWAligner}.cpp`
while chasing why `autocons` and the Rust port disagree. Recorded because the
design is a deliberate, non-obvious trade-off, not an accident — and because its
consequences show up in real output.

## The problem it solves

Farrar's published striped Smith-Waterman is **score-only**. The striped layout
processes query positions out of order (segment-interleaved), and the lazy-F
correction pass rewrites `F` values *after* a column is nominally finished, so
there is no natural point at which a per-cell traceback direction is final.
parasail's answer is a separate `_trace_` kernel family that materialises a full
`m × n` traceback matrix and costs both time and memory.

GIRI took a different route: emit traceback from the striped kernel itself, using
one `short` per cell.

## The encoding

`pathMatrix` holds one `short` per DP cell (`SWAligner::traceBack`):

| value | meaning |
|---|---|
| `-0x8000` | local-alignment termination — stop here |
| `0` | diagonal step (aligned pair) |
| `v > 0` | a gap of length `v` in the **top** sequence |
| `v < 0` | a gap of length `-v` in the **bottom** sequence |

The critical part is that a cell stores a whole **run length**, not a direction.
Traceback therefore *jumps* an entire gap run in one step:

```c
if (pathMat[i] > 0) {
   short gaplen = pathMat[i];
   while (gaplen-- > 0) {           // consume the whole run at once
      topAln.push_back(GAPCHAR);
      botAln.push_back(index2char[botSeq[--x]]);
   }
}
```

Run lengths are maintained in SIMD alongside the scores — `vEGL` for gaps in the
query, `vFGL` for gaps in the database:

```c
vMSK = _mm_cmpgt_epi16(vE, vH);   // still extending?
vEGL = _mm_and_si128(vEGL, vMSK); // zero the counter where the gap re-opens
vEGL = _mm_add_epi16(vEGL, vONE); // otherwise lengthen it
```

and the lazy-F correction folds the F run length into the direction word with a
neat trick — subtracting so that "F wins" yields `0 - vFGL`, i.e. the negative
encoding, while "F loses" keeps the previous value:

```c
vMSK = _mm_cmpgt_epi16(vF, vH);
vDIR = _mm_and_si128(vFGL, vMSK);              // vFGL where F wins, else 0
vMSK = _mm_andnot_si128(vMSK, *(pvDIR + j));   // old DIR where F loses, else 0
vDIR = _mm_subs_epi16(vMSK, vDIR);             // oldDIR - vFGL
```

## Why this is a reasonable choice

- **Two bytes per cell**, against parasail's separate trace kernels and their
  `O(mn)` trace tables — at 10 kb × 10 kb that difference is hundreds of MB.
- **Traceback becomes run-length decompression** — a handful of branches, no
  re-derivation of the DP.
- Most importantly, it gets a traceback out of the striped layout *at all*,
  without the second pass Farrar's original design forces.

## What it costs

The reconstructed path is not guaranteed to be the maximum-scoring one.

The run-length-per-cell encoding is only faithful if the recorded run length at a
cell matches the run the optimal path actually takes through it. Two things
undermine that:

1. **The counters record the current best `E`/`F` run, not the run on the final
   optimal path.** A cell whose best `E` state came from a length-5 run can be
   entered by an optimal path that used a length-2 run then a diagonal; the
   traceback jumps 5 regardless.
2. **The lazy-F pass rewrites `F` after the fact**, and `vFGL` is shifted between
   segments (`vFGL = _mm_slli_si128(vFGL, 2)`) with a zero shifted in — unlike
   `vF`, which gets `_mm_or_si128(vF, vMin)`. A gap run crossing a segment
   boundary loses its recorded length.

Because traceback jumps rather than stepping cell by cell, it cannot detect or
recover from either.

## Which score feeds reference selection

`swStripedWord` ends:

```c
double score = (short) iMaxScore + SHORT_BIAS;   // the DP maximum
if (score) {
   traceBack(topSeq, botSeq, pathMatrix, iMaxX, iMaxY, top, bot);
}
return score;                                     // traceBack's return is discarded
```

`SWAligner::align(PairwiseAlignment&)` stores that DP maximum via
`aln.setScore(...)`; `alignOneToMany` accumulates `args->score += aln.getScore()`;
`MultipleAlignment::align` returns the sum; `construct` returns it;
`autoConsensus` assigns it to `outscore`; `autocons` ranks candidate references
on it.

**So reference selection uses the Smith-Waterman DP maximum, never the score of
the traceback that is actually placed in the MSA.** Ranking and MSA content come
from two different alignments whenever the two diverge.

## Observed consequence — and which aligner it applies to

Instrumenting `ThreadedAligner::alignOneToMany`
(round-4/family-1339, reference `gi|2` 1346 bp against `gi|1` 1367 bp,
`lightWeightMatrix`, gap 10/2):

| aligner | columns | reported score | actual score of the emitted alignment |
|---|---|---|---|
| **Monardo** (`--mon`, and the *default* — see below) | 1478 | **918** | **364** |
| **Farrar** (`--sse2`) | 1455 | 531 | ~566 (within hand-scoring noise) |

The large divergence belongs to **Monardo**, not to the striped Farrar code
described above. Farrar's reported score and emitted alignment agree closely.
The run-length traceback design documented here is a real trade-off, but on this
evidence it is not producing badly suboptimal paths.

## The default aligner is Monardo, not Farrar

`autocons.cpp` selects:

```c
if(pAligner.get() == NULL) {
#ifdef HAVE_SSE2
   pAligner.reset(new FarrarSWAligner());
#else
   pAligner.reset(new MonardoSWAligner());
#endif
}
```

but `config.h` is only included under a guard:

```c
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
```

and `build-acons.sh` compiles `autocons` **without `-DHAVE_CONFIG_H`**. So
`HAVE_SSE2` is never defined, the `#else` branch is taken, and the shipped binary
defaults to the scalar **Monardo** aligner — even though `--help` says of
`--sse2`, *"This is the default alignment."*

Verified by instrumentation: default and `--mon` produce byte-identical results
(918 / total 7038); `--sse2` matches a standalone Farrar harness exactly
(531 / total 4401).

Two consequences:

- Comparisons against `autocons` **must pass `--sse2`** to compare like with
  like. Against the default one is comparing different algorithms. Over all 791
  families, byte-identical consensus went from **95/791 (12.0%)** against the
  Monardo default to **251/791 (31.7%)** with `--sse2` — and the systematic
  "port produces longer consensi" bias vanished (p90 length ratio 1.219 → 1.000).
- Part of `autocons`'s runtime is the scalar Monardo aligner, not SIMD. Measured
  speedup of the Rust port drops from 4.7× against the default to **1.6×**
  against `--sse2`.

**This looks like a build regression, not a design decision** — the aligner
selection clearly intends Farrar where SSE2 is available, and `--help` documents
`--sse2` as the default. The fix is to add `-DHAVE_CONFIG_H` to the `autocons`
compile in `build-acons.sh` (it already passes `-I` to where `config.h` lives),
or to drop the guard.

## Checking any autocons binary

`tools/which-aligner.sh <path-to-autocons>` reports which aligner a binary
defaults to. Exit status 0 = Farrar (intended), 1 = Monardo (regression),
2 = inconclusive.

```sh
$ tools/which-aligner.sh /usr/local/bin/autocons
Monardo        <- BUILD REGRESSION: config.h not included, so HAVE_SSE2 undefined
                  rebuild with -DHAVE_CONFIG_H (see build-acons.sh)
```

It runs the binary three ways on `tools/aligner-probe.fa` — default, `--mon`,
`--sse2` — and reports which the default matches byte for byte.

**Why not `nm` or `strings`.** Both aligner classes are linked into every build,
since each is reachable through its own flag, so symbol presence discriminates
nothing. The default, however, must match exactly one of them.

The probe family is chosen so `--mon` and `--sse2` genuinely disagree; the script
checks that and reports `INCONCLUSIVE` rather than guessing if they happen to
agree on some other input you pass it. Validated against both a binary built by
`build-acons.sh` (reports Monardo) and one built with `-DHAVE_CONFIG_H`
(reports Farrar).

A weaker corroborating signal: Monardo is scalar and roughly 3x slower, so a
default run that takes far longer than the same run with `--sse2` is suspicious.

## Instrumenting this yourself

`alignOneToMany` is a template in a header, so patch **copies** and shadow them
with `-I`:

```sh
cp giri_cpp_lib/src/libbio/{ThreadedAligner,MultipleAlignment}.h inc/
# patch inc/ThreadedAligner.h, then build with -Iinc first
```

Both headers must be copied. A quoted `#include "ThreadedAligner.h"` inside
`MultipleAlignment.h` resolves relative to *that file's* directory first, so
shadowing only `ThreadedAligner.h` silently has no effect — the build succeeds
and prints nothing.

## Resolved

The earlier "same aligner, different alignment inside autocons" puzzle was not a
puzzle: it was never the same aligner. `autocons` was running Monardo while the
standalone harness ran Farrar. Ruled out along the way, and worth not
re-checking: the sequence reader (`IstreamWrapper` gives identical bytes),
threading mode, optimisation level `-O0`…`-O3`, cross-alignment buffer
contamination, and post-read sequence mutation.
