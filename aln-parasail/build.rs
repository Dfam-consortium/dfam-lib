//! Build the vendored parasail subset.
//!
//! parasail ships CMake, autotools and meson builds.  None of them are used
//! here: only the striped kernels are vendored (under `third_party/parasail/`)
//! — traceback and score-only —
//! and their build
//! requirements are simple enough to state directly, which keeps CMake off the
//! list of things a contributor has to install.
//!
//! Two things this has to get right:
//!
//! 1. **Per-file ISA flags.**  Each kernel uses intrinsics for exactly one
//!    instruction set, so `sse41` files get `-msse4.1` and `avx2` files get
//!    `-mavx2`.  Compiling the whole library with `-mavx2` would let the
//!    compiler emit AVX2 into the SSE2 kernels and segfault on older hardware.
//!    They are therefore built as three separate static libraries.
//! 2. **No `-march=native`.**  The binary must run on any x86-64 host; which
//!    kernel executes is decided at run time in `src/isa.rs`.

use std::path::{Path, PathBuf};

/// Kernels vendored for each `(algorithm, ISA, lane width)` combination.
const ALGORITHMS: &[&str] = &["nw", "sg", "sw"];
const WIDTHS: &[&str] = &["8", "16", "32"];

/// `(vendor infix, cc flag, output library name)`
const ISAS: &[(&str, &str, &str)] = &[
    ("sse2_128", "-msse2", "parasail_sse2"),
    ("sse41_128", "-msse4.1", "parasail_sse41"),
    ("avx2_256", "-mavx2", "parasail_avx2"),
];

/// Shared support code.  Compiled at baseline SSE2 — `memory_sse.c` uses SSE2
/// intrinsics, `memory_avx2.c` needs AVX2 and is built with the AVX2 group.
const SUPPORT_BASE: &[&str] = &["memory.c", "cigar.c"];

fn main() {
    let vendor = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("third_party")
        .join("parasail");
    let src = vendor.join("src");

    if !src.join("memory.c").exists() {
        panic!(
            "vendored parasail sources are missing from {} — see \
             third_party/parasail/README.md",
            src.display()
        );
    }

    let target_arch = std::env::var("CARGO_CFG_TARGET_ARCH").unwrap_or_default();
    if target_arch != "x86_64" && target_arch != "x86" {
        panic!(
            "aln-parasail vendors only the x86 SIMD kernels, and this is \
             {target_arch}.  Build without this backend — disable the `parasail` \
             feature of the consuming crate, or drop aln-parasail from its \
             dependencies — or add the neon/altivec kernels under \
             third_party/parasail/ and extend the ISAS table below."
        );
    }

    // Baseline group: allocation and CIGAR helpers, plus the SSE2 profile
    // builders that every ISA path can fall back to.
    let mut base = cc::Build::new();
    configure(&mut base, &vendor);
    base.flag("-msse2");
    for f in SUPPORT_BASE {
        base.file(src.join(f));
    }
    base.file(src.join("memory_sse.c"));
    base.compile("parasail_support");

    // One library per ISA so the compiler can never hoist a wider instruction
    // into a narrower kernel.
    for (infix, flag, libname) in ISAS {
        let mut build = cc::Build::new();
        configure(&mut build, &vendor);
        build.flag(flag);

        if *infix == "avx2_256" {
            build.file(src.join("memory_avx2.c"));
        }
        for alg in ALGORITHMS {
            for w in WIDTHS {
                // Trace kernels write an O(mn) traceback matrix; the score-only
                // kernels keep O(m) column state and are what a score prepass
                // wants.  Both are vendored so `AlignParams::traceback` can
                // actually select between them.
                build.file(src.join(format!("{alg}_trace_striped_{infix}_{w}.c")));
                build.file(src.join(format!("{alg}_striped_{infix}_{w}.c")));
            }
        }
        build.compile(libname);
    }

    println!("cargo:rerun-if-changed=third_party/parasail");
    println!("cargo:rerun-if-changed=build.rs");
}

fn configure(build: &mut cc::Build, vendor: &Path) {
    build
        .include(vendor)
        .include(vendor.join("src"))
        .flag_if_supported("-std=gnu99")
        .define("_POSIX_C_SOURCE", "200112L")
        .warnings(false)
        .opt_level(3);
}
