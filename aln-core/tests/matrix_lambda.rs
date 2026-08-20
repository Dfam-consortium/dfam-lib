//! Pin the lambda port against RepeatMasker's `Matrix.pm`.
//!
//! Reference values produced by:
//!
//! ```sh
//! perl -Mlib=/usr/local/RepeatMasker -MMatrix -e '
//!   my $m = Matrix->new(fileName => shift);
//!   printf("%.10f\n", $m->getLambda());' <matrix>
//! ```
//!
//! `_calculateLambda` bisects to a tolerance of 1e-5 on lambda itself, so
//! agreement is asserted to 1e-6 — tight enough to catch a wrong recurrence or a
//! misplaced frequency, loose enough to survive the last bisection step landing
//! on a different side.
//!
//! These run against the installed RepeatMasker tree and are skipped when it is
//! absent, so the suite still passes on a machine without it.

use std::path::Path;

use aln_core::SubstMatrix;

const MATRIX_DIR: &str = "/usr/local/RepeatMasker/Matrices/crossmatch";

/// `(file stem, lambda from Matrix.pm)`
const REFERENCE: &[(&str, f64)] = &[
    ("14p35g", 0.128_395_080_6),
    ("25p41g", 0.123_222_351_1),
    ("20p43g", 0.122_673_034_7),
    ("14p41g", 0.125_328_064_0),
];

#[test]
fn lambda_matches_repeatmasker() {
    let dir = Path::new(MATRIX_DIR);
    if !dir.is_dir() {
        eprintln!("skipping: {MATRIX_DIR} not present");
        return;
    }

    let mut checked = 0;
    for &(stem, expected) in REFERENCE {
        let path = dir.join(format!("{stem}.matrix"));
        if !path.exists() {
            continue;
        }
        let m = SubstMatrix::from_file(&path)
            .unwrap_or_else(|e| panic!("parsing {}: {e}", path.display()));
        assert_eq!(m.name(), stem);

        let lambda = m
            .lambda()
            .unwrap_or_else(|| panic!("{stem} has FREQS but yielded no lambda"));
        assert!(
            (lambda - expected).abs() < 1e-6,
            "{stem}: lambda {lambda:.10} != Matrix.pm's {expected:.10}"
        );
        checked += 1;
    }
    assert!(checked > 0, "no reference matrices found under {MATRIX_DIR}");
}

/// Every shipped crossmatch matrix must parse, and any with `FREQS` must yield a
/// lambda that actually solves the Karlin-Altschul equation.
#[test]
fn all_crossmatch_matrices_parse_and_solve() {
    let dir = Path::new(MATRIX_DIR);
    if !dir.is_dir() {
        eprintln!("skipping: {MATRIX_DIR} not present");
        return;
    }

    let mut seen = 0;
    for entry in std::fs::read_dir(dir).unwrap() {
        let path = entry.unwrap().path();
        if path.extension().and_then(|e| e.to_str()) != Some("matrix") {
            continue;
        }
        let m = SubstMatrix::from_file(&path)
            .unwrap_or_else(|e| panic!("parsing {}: {e}", path.display()));
        seen += 1;

        let Some(lambda) = m.lambda() else { continue };
        let mut sum = 0.0;
        for (i, &fi) in m.freqs().iter().enumerate() {
            for (j, &fj) in m.freqs().iter().enumerate() {
                if fi > 0.0 && fj > 0.0 {
                    sum += fi * fj * (lambda * m.score_idx(i, j) as f64).exp();
                }
            }
        }
        assert!(
            (sum - 1.0).abs() < 1e-3,
            "{}: sum f_i f_j exp(lambda S_ij) = {sum}, expected 1.0",
            path.display()
        );
    }
    assert!(seen > 0, "no .matrix files found under {MATRIX_DIR}");
}
