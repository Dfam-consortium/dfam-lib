# Vendored parasail subset

From **parasail 2.6.2** (`parasail-v2.6.2.tar.gz`), upstream
<https://github.com/jeffdaily/parasail>.

Copyright (c) 2015, Battelle Memorial Institute. Licensed under BSD-3-Clause
terms *plus* a name-use restriction and a citation request — see `LICENSE.md`,
retained here as that licence requires, and `../../../THIRD-PARTY.md` for what
the extra clauses mean in practice. `README.upstream.md` is parasail's own
README, kept for the function-naming table.

This directory holds one upstream project and one licence. A second vendored
upstream (a NEON shim, say) gets its own sibling directory under
`aln-parasail/third_party/`; do not add files here that `LICENSE.md` does not
cover.

## What is here, and why so little

parasail generates 595 `.c` files covering every algorithm × ISA × lane width.
Only the striped **traceback** kernels are needed:

| files | purpose |
|---|---|
| `src/{nw,sg,sw}_trace_striped_{sse2_128,sse41_128,avx2_256}_{8,16,32}.c` | the 27 kernels |
| `src/sg_helper.h` | expands the one `sg` source into all 11 semi-global variants |
| `src/memory.c`, `src/memory_sse.c`, `src/memory_avx2.c` | matrices, profiles, aligned allocation |
| `src/cigar.c`, `src/cigar_template.c` | traceback → CIGAR |
| `parasail.h`, `parasail/*.h` | public and internal headers |
| `config.h` | hand-written; see below |

Deliberately excluded:

- **`satcheck.c`** (2.1 MB) — parasail's `_sat` wrappers. Including it would pull
  in the 8- and 16-bit variants of every algorithm. The 8 → 16 → 32 bit fallback
  runs in Rust instead (`src/lib.rs`).
- **`sw_dispatch.c` / `cpuid.c` / `isastubs.c`** — runtime ISA dispatch. Rust
  uses `is_x86_feature_detected!`, which also handles OS-level AVX state
  correctly.
- **`matrix_lookup.c` and `parasail/matrices/`** — the built-in BLOSUM/PAM
  tables. Matrices come from RepeatMasker `.matrix` files via
  `aln_core::SubstMatrix`.
- Everything for PowerPC (`altivec`) and ARM (`neon`), the CLI apps, and the
  test suite.

## config.h

parasail normally generates this with CMake or autotools. Vendoring a fixed
subset makes the probe results knowable up front, so `config.h` is hand-written
and `build.rs` never has to run CMake. It declares all three x86 SIMD levels
because each kernel is compiled with its own `-m` flag; which one *runs* is a
runtime decision.

## Re-vendoring

To add kernels (banded, `scan`, `stats`, or the 64-bit widths), copy the matching
`src/*.c` files into `src/` here and extend the `ALGORITHMS` / `WIDTHS` / `ISAS` tables in
`build.rs`. To move to a new parasail release, re-copy the same file list and
re-check `config.h` against the upstream `cmake/config.h.in`.
