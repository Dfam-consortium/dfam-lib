//! What does the traceback kernel cost over score-only?
//!
//! Farrar fuses traceback into its scoring pass (Kokhanyy's 2009 modification
//! to `swStripedWord`), so a score cutoff can never save work there. parasail
//! keeps the two kernels separate, so a two-phase "score first, trace only the
//! survivors" is possible. This measures what that would buy.
//!
//! cargo run --release -p aln-parasail --example tracecost -- <fa> <matrix> [n]

use std::time::Instant;

use aln_core::{io, SubstMatrix};
use aln_engine::{AlignMode, AlignParams, PairwiseAligner};
use aln_parasail::ParasailAligner;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut a = std::env::args().skip(1);
    let fa = a.next().expect("fasta");
    let mx = a.next().expect("matrix");
    let n: usize = a.next().map(|v| v.parse().unwrap()).unwrap_or(20);

    let seqs = io::read_fasta(std::io::BufReader::new(std::fs::File::open(&fa)?))?;
    let seqs = &seqs[..n.min(seqs.len())];
    let matrix = SubstMatrix::parse(&std::fs::read_to_string(&mx)?)?;
    println!(
        "{} seqs, mean len {}, {} pairs",
        seqs.len(),
        seqs.iter().map(|s| s.len()).sum::<usize>() / seqs.len(),
        seqs.len() * seqs.len()
    );

    let mut p = AlignParams::from_matrix(&matrix);
    p.mode = AlignMode::Local;
    p.min_score = 1;
    let al = ParasailAligner::new(matrix, p)?;

    // Prepare each subject once, as the real pipeline does.
    let profiles: Vec<_> = seqs
        .iter()
        .map(|s| al.prepare_subject(s))
        .collect::<Result<_, _>>()?;

    let t = Instant::now();
    let mut sum = 0i64;
    for pr in &profiles {
        for q in seqs {
            if let Some(r) = al.align_prepared(pr, q)? {
                sum += r.score as i64;
            }
        }
    }
    let trace = t.elapsed().as_secs_f64();

    let t = Instant::now();
    let mut sum2 = 0i64;
    for pr in &profiles {
        for q in seqs {
            if let Some(s) = al.score_prepared(pr, q)? {
                sum2 += s as i64;
            }
        }
    }
    let score = t.elapsed().as_secs_f64();

    assert_eq!(sum, sum2, "the two kernel families must agree");
    println!("  traceback  {trace:.3}s");
    println!("  score-only {score:.3}s   ({:.2}x faster)", trace / score);
    Ok(())
}
