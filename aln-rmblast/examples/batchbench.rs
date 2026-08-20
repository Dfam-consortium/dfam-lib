//! Where does the batched search time actually go? Compares one all-vs-all call
//! against the pair-at-a-time loop it replaced.
//!
//! cargo run --release -p aln-rmblast --example batchbench -- <fa> <matrix>
use std::time::Instant;

use aln_core::{io, SubstMatrix};
use aln_engine::engine::{SearchEngine, SearchParams, SeqSource};
use aln_rmblast::{RmblastEngine, RmblastOptions};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut a = std::env::args().skip(1);
    let fa = a.next().expect("fasta");
    let mx = a.next().expect("matrix");
    let seqs = io::read_fasta(std::io::BufReader::new(std::fs::File::open(&fa)?))?;
    let matrix = SubstMatrix::parse(&std::fs::read_to_string(&mx)?)?;
    let min_score: i32 = a.next().map(|v| v.parse().unwrap()).unwrap_or(1);
    let n: usize = a.next().map(|v| v.parse().unwrap()).unwrap_or(seqs.len());
    let seqs = &seqs[..n.min(seqs.len())];
    println!("{} seqs, mean len {}", seqs.len(),
        seqs.iter().map(|s| s.len()).sum::<usize>() / seqs.len());

    let p = SearchParams {
        matrix: Some(matrix),
        gap_init: -25,
        ins_gap_ext: -5,
        del_gap_ext: -5,
        min_match: 7,
        min_score,
        mask_level: 80,
        cores: Some(1),
        ..Default::default()
    };
    let e = RmblastEngine::new(p, RmblastOptions::default())?;

    let t = Instant::now();
    let hits = e.all_vs_all(seqs, 101)?;
    println!("all_vs_all: {:.2}s, {} hsps", t.elapsed().as_secs_f64(), hits.len());
    let mut sc: Vec<i32> = hits.iter().map(|(_, _, a)| a.score).collect();
    sc.sort_unstable_by(|a, b| b.cmp(a));
    let pct = |q: f64| sc.get(((sc.len() as f64 * q) as usize).min(sc.len().saturating_sub(1))).copied().unwrap_or(0);
    println!(
        "  n={} max={} p99={} p90={} median={} | >=150: {} ({:.3}%) | >=60: {} ({:.3}%)",
        sc.len(), sc.first().copied().unwrap_or(0), pct(0.01), pct(0.10), pct(0.50),
        sc.iter().filter(|&&v| v >= 150).count(),
        100.0 * sc.iter().filter(|&&v| v >= 150).count() as f64 / sc.len().max(1) as f64,
        sc.iter().filter(|&&v| v >= 60).count(),
        100.0 * sc.iter().filter(|&&v| v >= 60).count() as f64 / sc.len().max(1) as f64);
    if std::env::var("BENCH_ONLY_SCORES").is_ok() { return Ok(()); }

    let t = Instant::now();
    let mut k = 0;
    for r in seqs {
        k += e.one_to_many(r, seqs, None)?.len();
    }
    println!("one_to_many x{}: {:.2}s, {} hsps", seqs.len(), t.elapsed().as_secs_f64(), k);

    println!("--- min_score {min_score} ---");
    let t = Instant::now();
    let mut k = 0;
    for r in seqs.iter().take(4) {
        for q in seqs {
            k += e.search(
                &SeqSource::Memory(vec![q.clone()]),
                &SeqSource::Memory(vec![r.clone()]),
            )?.len();
        }
    }
    let per = t.elapsed().as_secs_f64() / 4.0;
    println!("pair loop (4 refs): {:.2}s/ref -> {:.1}s extrapolated, {} hsps",
        per, per * seqs.len() as f64, k);
    Ok(())
}
