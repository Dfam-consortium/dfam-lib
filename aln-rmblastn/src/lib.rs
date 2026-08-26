//! [`SearchEngine`] backed by the **external `rmblastn` binary**.
//!
//! The sibling to `aln-rmblast`, which is an in-process port of the same
//! algorithm. Both exist on purpose:
//!
//! * **This one** is what RepeatMasker actually runs, so its output is the
//!   reference for anything that must agree with published annotation. It is
//!   also the only path that can consume a prepared BLAST database, which the
//!   in-process port rejects.
//! * **The in-process port** needs no external binary, is far faster on small
//!   inputs (no process spawn, no temp files, no database build), and can be
//!   unit-tested.
//!
//! Selecting between them is a deployment decision, not a correctness one, so
//! both sit behind [`SearchEngine`] and callers choose.
//!
//! # Parameter mapping
//!
//! Flags follow `NCBIBlastSearchEngine.pm`, which is what RepeatMasker and
//! RepeatModeler's `Refiner` emit. Two conversions are easy to get wrong:
//!
//! * **Gap costs.** crossmatch charges `open + (k-1)*ext` for a gap of length
//!   `k`; NCBI charges `open + k*ext`. So NCBI's open must be reduced by one
//!   extension to describe the same scoring system — a crossmatch `gap_init` of
//!   -25 with `ext` -5 becomes `-gapopen 20 -gapextend 5`.
//! * **X-drop.** NCBI derives all three cutoffs from the score floor:
//!   `minScore*2`, `minScore/2`, `minScore`. They are not independent knobs
//!   here, which is why lowering `min_score` also shortens extensions.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use aln_core::{Alignment, Sequence, Strand};
use aln_engine::engine::{ScoreMode, SearchEngine, SearchParams, SeqSource};
use aln_engine::{EngineError, Result};

/// Backend name, for diagnostics and `.out` provenance lines.
pub const NAME: &str = "rmblastn";

/// Tabular fields requested from `rmblastn`, in order.
///
/// `qseq`/`sseq` are the aligned (gapped) strings — without them there is no
/// traceback to recover, only coordinates.
const OUTFMT: &str = "6 score qseqid qstart qend qlen sstrand sseqid sstart send slen qseq sseq";

/// Options with no [`SearchParams`] equivalent.
#[derive(Clone, Debug, Default)]
pub struct RmblastnOptions {
    /// `makeblastdb` executable, used when a subject must be indexed.
    ///
    /// Leave as `None` to take the one sitting **next to `rmblastn`**. That
    /// matters: a machine can carry more than one BLAST+ install (here,
    /// `/usr/local/rmblast/bin` is 2.17.1 while `/usr/bin` is 2.12.0), and
    /// building a database with one version then searching it with another is a
    /// silent-corruption hazard rather than a clean failure.
    pub makeblastdb: Option<PathBuf>,
    /// Low-complexity filtering. RepeatMasker runs `-dust no`.
    pub dust: bool,
    /// Matrix file to pass as `-matrix`, e.g. RepeatMasker's
    /// `Matrices/ncbi/nt/20p43g.matrix`.
    ///
    /// **Required if you want a matrix at all.** `SearchParams::matrix` holds a
    /// *parsed* matrix in crossmatch layout, and this is deliberately not
    /// converted for you: NCBI's matrices are the transpose of crossmatch's, so
    /// synthesising one risks writing a silently wrong scoring system. Point at
    /// a real NCBI-format file instead.
    ///
    /// Two rmblastn requirements this handles: the value must be a **bare
    /// filename** (an absolute path is rejected), resolved against `BLASTMAT`,
    /// which is set from this path's parent.
    pub matrix_path: Option<PathBuf>,
}

/// Locate `name` in the same directory as `exe`, resolving `exe` through `PATH`
/// when it is a bare command. Returns `None` if no such file exists.
fn sibling_of(exe: &Path, name: &str) -> Option<PathBuf> {
    let full = if exe.components().count() > 1 {
        exe.to_path_buf()
    } else {
        std::env::var_os("PATH")
            .map(|paths| std::env::split_paths(&paths).map(|d| d.join(exe)).find(|c| c.is_file()))?
            ?
    };
    let cand = full.parent()?.join(name);
    cand.is_file().then_some(cand)
}

/// A search engine that shells out to `rmblastn`.
pub struct RmblastnEngine {
    params: SearchParams,
    opts: RmblastnOptions,
    exe: PathBuf,
    /// Resolved at construction so it always matches `exe`'s install.
    makeblastdb: PathBuf,
}

impl RmblastnEngine {
    /// Fails immediately on a configuration that cannot work, rather than
    /// after temp files have been written and a database built.
    pub fn new(params: SearchParams, opts: RmblastnOptions) -> Result<Self> {
        if opts.matrix_path.is_none() && params.matrix.is_some() {
            return Err(EngineError::Params(
                "SearchParams::matrix is set but RmblastnOptions::matrix_path is not. \
                 The external rmblastn needs an NCBI-format matrix FILE, and NCBI's \
                 matrices are the transpose of the crossmatch layout SubstMatrix holds, \
                 so this is not converted automatically — writing the wrong orientation \
                 would corrupt every score silently. Point matrix_path at a real file \
                 (e.g. RepeatMasker Matrices/ncbi/nt/20p43g.matrix), or clear \
                 SearchParams::matrix to use rmblastn's default scoring."
                    .to_string(),
            ));
        }
        let exe = params
            .path_to_engine
            .clone()
            .unwrap_or_else(|| PathBuf::from("rmblastn"));
        let makeblastdb = opts
            .makeblastdb
            .clone()
            .or_else(|| sibling_of(&exe, "makeblastdb"))
            .unwrap_or_else(|| PathBuf::from("makeblastdb"));
        Ok(RmblastnEngine { params, opts, exe, makeblastdb })
    }

    /// NCBI's gap-open, converted from the crossmatch convention.
    fn gap_open(&self) -> Result<i32> {
        let ext = self.params.ins_gap_ext.unsigned_abs() as i32;
        (self.params.gap_init.unsigned_abs() as i32)
            .checked_sub(ext)
            .filter(|v| *v >= 0)
            .ok_or_else(|| {
                EngineError::Params(format!(
                    "gap_init ({}) must be at least as large in magnitude as the gap \
                     extension ({}); NCBI's open cost excludes the first position",
                    self.params.gap_init, self.params.ins_gap_ext
                ))
            })
    }

    /// Materialise a source as a FASTA path, writing a temp file if needed.
    ///
    /// The `TempDir` is returned so the caller keeps it alive; dropping it
    /// removes the file.
    fn as_fasta(
        &self,
        src: &SeqSource,
        what: &str,
    ) -> Result<(PathBuf, Option<tempfile::TempDir>)> {
        match src {
            SeqSource::Fasta(p) => Ok((p.clone(), None)),
            SeqSource::Memory(seqs) => {
                let dir = self.tempdir()?;
                let path = dir.path().join(format!("{what}.fa"));
                let mut f = std::fs::File::create(&path)
                    .map_err(|e| EngineError::backend(NAME, format!("creating {what} FASTA: {e}")))?;
                for s in seqs {
                    writeln!(f, ">{}", s.name)
                        .and_then(|_| f.write_all(&s.seq))
                        .and_then(|_| writeln!(f))
                        .map_err(|e| EngineError::backend(NAME, format!("writing {what}: {e}")))?;
                }
                Ok((path, Some(dir)))
            }
            SeqSource::TwoBit(p) => Err(EngineError::unsupported(
                NAME,
                format!(
                    "rmblastn reads FASTA and BLAST databases; convert {} with twoBitToFa first",
                    p.display()
                ),
            )),
            SeqSource::BlastDb(p) => Ok((p.clone(), None)),
        }
    }

    fn tempdir(&self) -> Result<tempfile::TempDir> {
        match &self.params.temp_dir {
            Some(d) => tempfile::TempDir::new_in(d),
            None => tempfile::TempDir::new(),
        }
        .map_err(|e| EngineError::backend(NAME, format!("creating a temp dir: {e}")))
    }

    /// Index a FASTA subject unless it is already a prepared database.
    fn ensure_db(&self, subject: &SeqSource, fasta: &Path) -> Result<PathBuf> {
        if matches!(subject, SeqSource::BlastDb(_)) {
            return Ok(fasta.to_path_buf());
        }
        let out = fasta.to_path_buf();
        let status = Command::new(&self.makeblastdb)
            .args(["-blastdb_version", "4", "-dbtype", "nucl", "-in"])
            .arg(fasta)
            .arg("-out")
            .arg(&out)
            .output()
            .map_err(|e| {
                EngineError::backend(
                    NAME,
                    format!("could not run {}: {e}", self.makeblastdb.display()),
                )
            })?;
        if !status.status.success() {
            return Err(EngineError::backend(
                NAME,
                format!(
                    "makeblastdb failed: {}",
                    String::from_utf8_lossy(&status.stderr).trim()
                ),
            ));
        }
        Ok(out)
    }
}

impl SearchEngine for RmblastnEngine {
    fn name(&self) -> &'static str {
        NAME
    }

    /// Reports both binaries, since they must come from the same install and a
    /// mismatch is otherwise invisible.
    fn version(&self) -> String {
        let ver = |p: &Path| {
            Command::new(p)
                .arg("-version")
                .output()
                .ok()
                .and_then(|o| String::from_utf8(o.stdout).ok())
                .and_then(|s| s.lines().next().map(str::to_string))
                .unwrap_or_else(|| format!("{} (version unavailable)", p.display()))
        };
        format!("{} / {}", ver(&self.exe), ver(&self.makeblastdb))
    }

    fn params(&self) -> &SearchParams {
        &self.params
    }

    fn accepts(&self, source: &SeqSource) -> bool {
        // The one thing this engine can do that the in-process port cannot.
        !matches!(source, SeqSource::TwoBit(_))
    }

    fn search(&self, query: &SeqSource, subject: &SeqSource) -> Result<Vec<Alignment>> {
        let (qpath, _qtmp) = self.as_fasta(query, "query")?;
        let (spath, _stmp) = self.as_fasta(subject, "subject")?;
        let db = self.ensure_db(subject, &spath)?;

        let min_score = self.params.min_score.max(1);
        let word = self.params.word_raw.unwrap_or(self.params.min_match);
        let gap_open = self.gap_open()?;
        let gap_ext = self.params.ins_gap_ext.unsigned_abs();

        let mut cmd = Command::new(&self.exe);
        cmd.arg("-query").arg(&qpath).arg("-db").arg(&db);
        cmd.args(["-outfmt", OUTFMT]);
        cmd.args(["-gapopen", &gap_open.to_string()]);
        cmd.args(["-gapextend", &gap_ext.to_string()]);
        cmd.args(["-word_size", &word.to_string()]);
        // NCBI ties the X-drop budget to the score floor; see the module docs.
        cmd.args(["-xdrop_ungap", &(min_score * 2).to_string()]);
        cmd.args(["-xdrop_gap", &(min_score / 2).to_string()]);
        cmd.args(["-xdrop_gap_final", &min_score.to_string()]);
        cmd.args(["-min_raw_gapped_score", &min_score.to_string()]);
        cmd.args(["-dust", if self.opts.dust { "yes" } else { "no" }]);
        cmd.args(["-num_threads", &self.params.cores.unwrap_or(1).max(1).to_string()]);
        if self.params.mask_level <= 100 {
            cmd.args(["-mask_level", &self.params.mask_level.to_string()]);
        }
        if self.params.score_mode == ScoreMode::ComplexityAdjusted {
            cmd.arg("-complexity_adjust");
        }
        match (&self.opts.matrix_path, &self.params.matrix) {
            (Some(m), _) => {
                // Bare filename plus BLASTMAT; rmblastn rejects an absolute path.
                let dir = m.parent().filter(|d| !d.as_os_str().is_empty());
                let file = m.file_name().ok_or_else(|| {
                    EngineError::Params(format!("matrix_path {} has no filename", m.display()))
                })?;
                if let Some(d) = dir {
                    cmd.env("BLASTMAT", d);
                }
                cmd.arg("-matrix").arg(file);
            }
            // `new` rejects matrix-without-path, so nothing to do here.
            (None, _) => {}
        }

        let out = cmd.output().map_err(|e| {
            EngineError::backend(NAME, format!("could not run {}: {e}", self.exe.display()))
        })?;
        if !out.status.success() {
            return Err(EngineError::backend(
                NAME,
                format!("rmblastn failed: {}", String::from_utf8_lossy(&out.stderr).trim()),
            ));
        }
        parse_tabular(&String::from_utf8_lossy(&out.stdout))
    }
}

/// Parse `-outfmt 6` rows in `OUTFMT` order into [`Alignment`]s.
///
/// BLAST reports 1-based closed coordinates and, on a minus-strand hit, gives
/// the subject range descending. This crate's convention is 0-based half-open
/// on the forward strand with the strand carried separately, so both are
/// normalised here rather than stored inverted.
pub fn parse_tabular(text: &str) -> Result<Vec<Alignment>> {
    let mut out = Vec::new();
    for (n, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            return Err(EngineError::backend(
                NAME,
                format!("line {}: expected 12 tabular fields, found {}", n + 1, f.len()),
            ));
        }
        let num = |i: usize, what: &str| -> Result<i64> {
            f[i].parse::<i64>().map_err(|_| {
                EngineError::backend(NAME, format!("line {}: bad {what} {:?}", n + 1, f[i]))
            })
        };
        let score = num(0, "score")? as i32;
        let qstart = num(2, "qstart")?;
        let sstart = num(7, "sstart")?;
        let send = num(8, "send")?;
        // f[5] is sstrand, f[6] is sseqid. The descending-coordinate test is a
        // fallback for output that omits or garbles the strand column; it must
        // not be the primary signal, since BLAST tabular can report a minus hit
        // with ascending subject coordinates.
        let strand = if f[5].eq_ignore_ascii_case("minus") || send < sstart {
            Strand::Minus
        } else {
            Strand::Plus
        };
        // 1-based closed -> 0-based half-open, and forward-strand for the subject.
        let subj_lo = sstart.min(send) - 1;
        let aln = Alignment::from_gapped(
            f[1],
            f[6],
            (qstart - 1).max(0) as usize,
            subj_lo.max(0) as usize,
            strand,
            score,
            f[10].as_bytes(),
            f[11].as_bytes(),
        )
        .map_err(|e| EngineError::backend(NAME, format!("line {}: {e}", n + 1)))?;
        out.push(aln);
    }
    Ok(out)
}

/// Convenience: run against sequences already in memory.
pub fn search_memory(
    engine: &RmblastnEngine,
    query: &[Sequence],
    subject: &[Sequence],
) -> Result<Vec<Alignment>> {
    engine.search(
        &SeqSource::Memory(query.to_vec()),
        &SeqSource::Memory(subject.to_vec()),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Parsing must not need the binary present — that is the whole point of
    /// keeping it a free function.
    #[test]
    fn tabular_rows_become_alignments() {
        // score q qs qe qlen strand s ss se slen qseq sseq
        let row = "120\tqry\t11\t20\t100\tplus\tsub\t5\t14\t50\tACGT-ACGTA\tACGTAACG-A";
        let a = parse_tabular(row).unwrap();
        assert_eq!(a.len(), 1);
        assert_eq!(a[0].score, 120);
        assert_eq!(a[0].query_name, "qry");
        assert_eq!(a[0].subj_name, "sub");
        assert_eq!(a[0].query_start, 10, "1-based closed becomes 0-based");
        assert_eq!(a[0].subj_start, 4);
        assert_eq!(a[0].strand, Strand::Plus);
    }

    /// A minus-strand hit arrives with the subject range descending; the stored
    /// coordinates must still satisfy start < end, with the strand separate.
    #[test]
    fn minus_strand_coordinates_are_normalised() {
        let row = "90\tqry\t1\t10\t100\tminus\tsub\t40\t31\t50\tACGTACGTAC\tACGTACGTAC";
        let a = parse_tabular(row).unwrap();
        assert_eq!(a[0].strand, Strand::Minus);
        assert_eq!(a[0].subj_start, 30, "the lower coordinate, 0-based");
        assert!(a[0].subj_start < a[0].subj_end);
    }

    /// Strand must come from the `sstrand` column, not be inferred from
    /// coordinate order. This row is minus with ASCENDING subject coordinates,
    /// so the `send < sstart` fallback cannot save a wrong field index.
    #[test]
    fn strand_is_read_from_the_strand_column() {
        let row = "90\tqry\t1\t10\t100\tminus\tsub\t31\t40\t50\tACGTACGTAC\tACGTACGTAC";
        let a = parse_tabular(row).unwrap();
        assert_eq!(a[0].strand, Strand::Minus, "strand must come from f[5]");
        assert_eq!(a[0].subj_name, "sub", "subject name must come from f[6]");
    }

    #[test]
    fn blank_and_comment_lines_are_skipped() {
        assert!(parse_tabular("\n# comment\n\n").unwrap().is_empty());
    }

    #[test]
    fn a_short_row_is_an_error_not_a_panic() {
        assert!(parse_tabular("1\t2\t3").is_err());
    }
}
