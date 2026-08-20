/* Hand-written replacement for parasail's cmake/autotools-generated config.h.
 *
 * parasail normally probes the host with CMake or autotools.  Vendoring only
 * the striped-traceback kernels means the probe results are knowable up front,
 * so this file states them directly and `build.rs` never has to run CMake.
 *
 * Scope: the three x86 SIMD levels are compiled unconditionally; which one runs
 * is chosen at *run time* by `is_x86_feature_detected!` on the Rust side, not
 * by parasail's own cpuid dispatcher (which is not vendored).  Each kernel is
 * built with its own `-m` flag, so declaring all three here is correct even on
 * a host that supports only SSE2.
 *
 * If you vendor additional kernels, revisit HAVE_ZLIB / HAVE_GETOPT / the
 * timing probes — nothing currently compiled needs them.
 */

#ifndef PARASAIL_VENDORED_CONFIG_H
#define PARASAIL_VENDORED_CONFIG_H

/* Aligned allocation.  build.rs compiles with -std=gnu99 and _POSIX_C_SOURCE
 * high enough for posix_memalign, which is available on every platform this
 * crate targets (Linux, macOS).  Windows would want HAVE__ALIGNED_MALLOC. */
#ifdef _WIN32
#define HAVE__ALIGNED_MALLOC 1
#else
#define HAVE_POSIX_MEMALIGN 1
#endif

#define SIZEOF_INT 4

/* SIMD levels.  Each vendored kernel is compiled with a matching -m flag, so
 * the intrinsics are always available to the translation unit that needs them.
 * Runtime selection happens in Rust. */
#define HAVE_SSE2 1
#define HAVE_SSE2_MM_SET1_EPI64X 1
#define HAVE_SSE2_MM_SET_EPI64X 1

#define HAVE_SSE41 1
#define HAVE_SSE41_MM_INSERT_EPI64 1
#define HAVE_SSE41_MM_EXTRACT_EPI64 1

#define HAVE_AVX2 1
#define HAVE_AVX2_MM256_SET1_EPI64X 1
#define HAVE_AVX2_MM256_SET_EPI64X 1
#define HAVE_AVX2_MM256_INSERT_EPI64 1
#define HAVE_AVX2_MM256_INSERT_EPI32 1
#define HAVE_AVX2_MM256_INSERT_EPI16 1
#define HAVE_AVX2_MM256_INSERT_EPI8 1
#define HAVE_AVX2_MM256_EXTRACT_EPI64 1
#define HAVE_AVX2_MM256_EXTRACT_EPI32 1
#define HAVE_AVX2_MM256_EXTRACT_EPI16 1
#define HAVE_AVX2_MM256_EXTRACT_EPI8 1

/* Not vendored: PowerPC and ARM kernels. */
#define HAVE_ALTIVEC 0
#define HAVE_NEON 0

/* Not vendored: the cpuid dispatcher.  Rust does feature detection. */
#define HAVE_XGETBV 0

/* Not vendored: file I/O, the CLI apps, or the benchmark timers. */
#define HAVE_ZLIB 0
#define WORDS_BIGENDIAN 0

#define HAVE_UNISTD_H 1

#define INT64_LITERAL_SUFFIX_I64 0
#define INT64_LITERAL_SUFFIX_LL 1

#endif /* PARASAIL_VENDORED_CONFIG_H */
