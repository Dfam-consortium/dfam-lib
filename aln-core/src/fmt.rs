//! RepeatMasker-compatible output formats.
//!
//! Ports of `SearchResult.pm`'s `_toCrossMatchFormat` and `_toOUTFileFormat`.
//! Both are byte-for-byte reproductions, including the field widths, so output
//! can be diffed directly against RepeatMasker's.
//!
//! # Orientation in the annotation line
//!
//! Both formats print subject coordinates in the subject's own forward frame.
//! On the minus strand they are printed **descending**, with the remaining-base
//! count moved to the front:
//!
//! ```text
//!  plus:   subjName  begin  end   (left)
//!  minus:  C subjName  (left)  end   begin
//! ```
//!
//! [`Alignment`](crate::Alignment) always stores the span ascending, so the
//! reordering happens here and nowhere else.

use std::fmt::Write as _;

use crate::error::{Error, Result};
use crate::result::SearchResult;
use crate::seq::{self, Strand};

/// How much of the alignment to print, and in which orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AlignmentMode {
    /// Annotation line only — `SearchResult::N_NoAlign`.
    #[default]
    NoAlign,
    /// Show the alignment with the **query** read forward.  On the minus strand
    /// the subject is reverse-complemented and its line carries the `C` marker.
    /// `SearchResult::N_AlignWithQuerySeq`.
    WithQuerySeq,
    /// Show the alignment with the **subject** read forward.  On the minus
    /// strand the query is reverse-complemented and carries the `C` marker.
    /// `SearchResult::N_AlignWithSubjSeq`.
    WithSubjSeq,
}

/// Columns of aligned sequence per block.
const LINE_WIDTH: usize = 50;

/// Mutation marker for a mismatched column, from `SearchResult.pm`'s
/// `%mutChar`: `i` transition, `v` transversion, `-` gap, `?` ambiguous.
///
/// The table is keyed `query` then `subject`, and covers only unambiguous
/// pairings plus every gap combination; anything else falls through to `?`.
fn mut_char(q: u8, s: u8) -> Option<char> {
    let q = q.to_ascii_uppercase();
    let s = s.to_ascii_uppercase();
    if seq::is_gap(q) || seq::is_gap(s) {
        // Only pairings where the other side is a known symbol are in the
        // table, but every IUPAC code is listed, so any gap pairing maps to '-'.
        return Some('-');
    }
    match (q, s) {
        (b'C', b'T') | (b'T', b'C') | (b'A', b'G') | (b'G', b'A') => Some('i'),
        (b'G', b'T') | (b'T', b'G') | (b'G', b'C') | (b'C', b'G')
        | (b'C', b'A') | (b'A', b'C') | (b'A', b'T') | (b'T', b'A') => Some('v'),
        _ => None,
    }
}

fn is_iupac_ambiguity(b: u8) -> bool {
    matches!(
        b.to_ascii_uppercase(),
        b'B' | b'D' | b'H' | b'V' | b'R' | b'Y' | b'K' | b'M' | b'S' | b'W' | b'N' | b'X'
    )
}

// ── .out ──────────────────────────────────────────────────────────────────────

/// One `.out` annotation line.
///
/// Port of `_toOUTFileFormat`, field widths included:
/// `%6d %4.1f %4.1f %4.1f %-17s %8d %8d %8s %1s %-15s %-15s %7s %7s %7s %-5s %3s %3s`
pub fn to_out_line(r: &SearchResult) -> String {
    let a = &r.alignment;
    let (q_b, q_e) = a.query_one_based();
    let (s_b, s_e) = a.subject_one_based();
    let s_left = r.subj_left();

    let (orient, c1, c2, c3) = if a.strand == Strand::Minus {
        ("C", format!("({s_left})"), s_e.to_string(), s_b.to_string())
    } else {
        ("+", s_b.to_string(), s_e.to_string(), format!("({s_left})"))
    };

    format!(
        "{:>6} {:>4.1} {:>4.1} {:>4.1} {:<17} {:>8} {:>8} {:>8} {:1} {:<15} {:<15} \
         {:>7} {:>7} {:>7} {:<5} {:>3} {:>3}\n",
        a.score,
        r.pct_diverge,
        r.pct_delete,
        r.pct_insert,
        a.query_name,
        q_b,
        q_e,
        format!("({})", r.query_left()),
        orient,
        a.subj_name,
        r.subj_class.as_deref().unwrap_or(""),
        c1,
        c2,
        c3,
        r.id.map(|v| v.to_string()).unwrap_or_default(),
        r.lineage_id.as_deref().unwrap_or(""),
        r.overlap.map(|c| c.to_string()).unwrap_or_default(),
    )
}

// No `out_header()` is provided on purpose.  RepeatMasker's `.out` header is
// written by `ProcessRepeats`, not by `SearchResult.pm`, and its column widths
// are computed from the widest value actually present in the file
// (`$colWidths{'class'}`, `$colWidths{'Seq2Begin'}`, …).  A fixed header string
// would only line up by coincidence.  [`to_out_line`] is a faithful port of
// `_toOUTFileFormat`, which does use fixed widths; matching ProcessRepeats'
// adaptive table would mean porting its two-pass width calculation as well.

// ── CAF / CIGAR records ───────────────────────────────────────────────────────

/// The comma-separated header shared by the CAF and CIGAR records.
///
/// ```text
/// score,div,del,ins,qryName,qryStart,qryEnd,qryLeft,
/// subjName,subjType,subjStart,subjEnd,subjLeft,orient,overlap,id,<payload>
/// ```
///
/// `orient` is `1` on the minus strand and `0` otherwise. When the subject name
/// contains a `#` it is split into name and type there, and
/// [`SearchResult::subj_class`] is ignored — matching the Perl's regex branch.
fn record_prefix(r: &SearchResult) -> String {
    let a = &r.alignment;
    let (q_b, q_e) = a.query_one_based();
    let (s_b, s_e) = a.subject_one_based();

    let (s_name, s_type) = match a.subj_name.split_once('#') {
        Some((n, t)) => (n.to_string(), t.to_string()),
        None => (
            a.subj_name.clone(),
            r.subj_class.clone().unwrap_or_default(),
        ),
    };

    format!(
        "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},",
        a.score,
        r.pct_diverge,
        r.pct_delete,
        r.pct_insert,
        a.query_name,
        q_b,
        q_e,
        r.query_left(),
        s_name,
        s_type,
        s_b,
        s_e,
        r.subj_left(),
        u8::from(a.strand == Strand::Minus),
        r.overlap.map(|c| c.to_string()).unwrap_or_default(),
        r.id.map(|v| v.to_string()).unwrap_or_default(),
    )
}

/// Compressed Alignment Format — a port of `_toCAF`.
///
/// The alignment is folded into a single string: identical columns are the base
/// itself, a substitution is `query/subject`, a run of subject gaps is bracketed
/// by `-` and carries the query bases, and a run of query gaps is bracketed by
/// `+` and carries the subject bases.
///
/// ```text
/// query   AAGACTT---A
/// subject AAT--CTAATA   ->   AAG/T-AC-T/CT+AAT+A
/// ```
///
/// A gap run still open when the alignment ends is **not** closed, because the
/// Perl's loop simply terminates — `ACGTACGT--` against `ACGTACGTAC` yields
/// `ACGTACGT+AC`, with no trailing `+`. Reproduced exactly.
///
/// Returns no trailing newline.
pub fn to_caf(r: &SearchResult, query: &[u8], subject: &[u8]) -> Result<String> {
    let (q, s) = r.alignment.gapped(query, subject)?;
    let mut payload = String::with_capacity(q.len() * 2);
    let mut in_ins = false; // run of query gaps, bracketed by '+'
    let mut in_del = false; // run of subject gaps, bracketed by '-'

    for (&qc, &sc) in q.iter().zip(&s) {
        let (q_gap, s_gap) = (seq::is_gap(qc), seq::is_gap(sc));
        if q_gap {
            // The Perl closes an open deletion first — it has seen alignments
            // where the two run kinds abut with no aligned column between.
            if in_del {
                payload.push('-');
                in_del = false;
            }
            if !in_ins {
                payload.push('+');
                in_ins = true;
            }
            payload.push(sc as char);
        } else if s_gap {
            if in_ins {
                payload.push('+');
                in_ins = false;
            }
            if !in_del {
                payload.push('-');
                in_del = true;
            }
            payload.push(qc as char);
        } else {
            if in_del {
                payload.push('-');
                in_del = false;
            } else if in_ins {
                payload.push('+');
                in_ins = false;
            }
            if qc == sc {
                payload.push(qc as char);
            } else {
                payload.push(qc as char);
                payload.push('/');
                payload.push(sc as char);
            }
        }
    }

    Ok(record_prefix(r) + &payload)
}

/// CIGAR record — the same header as [`to_caf`] with a CIGAR payload.
///
/// # This deliberately does *not* reproduce `_toCIGAR`'s output
///
/// RepeatMasker's `_toCIGAR` emits each run's length paired with the *next*
/// run's operator, and opens with a spurious zero-length run:
///
/// ```text
/// query   AAGACTT---A
/// subject AAT--CTAATA
///
///   _toCIGAR emits:              0M3D2M2I3M1M
///   its own doc comment says:    3M2D2M3I1M
/// ```
///
/// Ten aligned columns with no gaps come out as `0M10M`. The output cannot be
/// parsed by anything, so reproducing it faithfully would serve no one; this
/// emits what RepeatMasker documents instead. `caf_and_cigar_perl_parity.rs`
/// pins the broken strings so the divergence stays visible.
///
/// # Operator orientation is inverted relative to SAM
///
/// RepeatMasker uses `I` for a gap in the **query** and `D` for a gap in the
/// **subject** — the opposite of SAM, and therefore of
/// [`EditScript::to_cigar`](crate::EditScript::to_cigar). That is the format's
/// definition rather than a defect, so it is preserved here. Use
/// `EditScript::to_cigar` when you want SAM semantics.
///
/// Returns no trailing newline.
pub fn to_cigar_record(r: &SearchResult, query: &[u8], subject: &[u8]) -> Result<String> {
    let (q, s) = r.alignment.gapped(query, subject)?;

    let mut payload = String::new();
    let mut run_op: Option<char> = None;
    let mut run_len = 0u32;

    for (&qc, &sc) in q.iter().zip(&s) {
        let op = if seq::is_gap(qc) {
            'I'
        } else if seq::is_gap(sc) {
            'D'
        } else {
            'M'
        };
        match run_op {
            Some(prev) if prev == op => run_len += 1,
            Some(prev) => {
                let _ = write!(payload, "{run_len}{prev}");
                run_op = Some(op);
                run_len = 1;
            }
            None => {
                run_op = Some(op);
                run_len = 1;
            }
        }
    }
    if let Some(op) = run_op {
        let _ = write!(payload, "{run_len}{op}");
    }

    Ok(record_prefix(r) + &payload)
}

// ── crossmatch ────────────────────────────────────────────────────────────────

/// Full crossmatch-style output: annotation line, and optionally the alignment.
///
/// `query` and `subject` are the **full forward** source sequences; the aligned
/// spans are sliced out here.
pub fn to_crossmatch(
    r: &SearchResult,
    query: &[u8],
    subject: &[u8],
    mode: AlignmentMode,
) -> Result<String> {
    let mut out = annotation_line(r);
    if mode == AlignmentMode::NoAlign {
        return Ok(out);
    }
    out.push('\n');
    alignment_block(&mut out, r, query, subject, mode)?;
    Ok(out)
}

/// The crossmatch annotation line, without a trailing alignment.
pub fn annotation_line(r: &SearchResult) -> String {
    let a = &r.alignment;
    let (q_b, q_e) = a.query_one_based();
    let (s_b, s_e) = a.subject_one_based();

    // Only the first whitespace-delimited token of each name is printed.
    let q_name = a.query_name.split_whitespace().next().unwrap_or("");
    let s_name = a.subj_name.split_whitespace().next().unwrap_or("");

    let mut s = format!(
        "{} {} {} {} {} {} {} ({}) ",
        a.score, r.pct_diverge, r.pct_delete, r.pct_insert, q_name, q_b, q_e,
        r.query_left()
    );
    if a.strand == Strand::Minus {
        let _ = write!(s, "C {} ({}) {} {}", s_name, r.subj_left(), s_e, s_b);
    } else {
        let _ = write!(s, "{} {} {} ({})", s_name, s_b, s_e, r.subj_left());
    }
    if let Some(l) = &r.lineage_id {
        let _ = write!(s, " {l}");
    }
    if let Some(id) = r.id {
        let _ = write!(s, " {id}");
    }
    if let Some(o) = r.overlap {
        let _ = write!(s, " {o}");
    }
    s.push('\n');
    s
}

/// Per-block tallies carried through the alignment walk.
#[derive(Default)]
struct Tally {
    transitions: u32,
    transversions: u32,
    ambiguous: u32,
    gap_columns: u32,
    gap_openings: u32,
}

fn alignment_block(
    out: &mut String,
    r: &SearchResult,
    query: &[u8],
    subject: &[u8],
    mode: AlignmentMode,
) -> Result<()> {
    let a = &r.alignment;
    // `gapped` returns the query forward with the subject reverse-complemented
    // on the minus strand — exactly RepeatMasker's stored orientation.
    let (mut q, mut s) = a.gapped(query, subject)?;
    if q.len() != s.len() {
        return Err(Error::Alignment("gapped rows differ in length".into()));
    }

    let minus = a.strand == Strand::Minus;
    let flip = minus && mode == AlignmentMode::WithSubjSeq;
    if flip {
        q = seq::revcomp(&q);
        s = seq::revcomp(&s);
    }

    let (q_b, q_e) = a.query_one_based();
    let (s_b, s_e) = a.subject_one_based();
    let mut q_start: i64 = if flip { q_e as i64 } else { q_b as i64 };
    let mut s_start: i64 = if minus && mode == AlignmentMode::WithQuerySeq {
        s_e as i64
    } else {
        s_b as i64
    };
    let mut q_end: i64 = 0;
    let mut s_end: i64 = 0;

    let mut t = Tally::default();
    let q_name: String = a.query_name.chars().take(13).collect();
    let s_name: String = a.subj_name.chars().take(13).collect();

    for chunk in 0..q.len().div_ceil(LINE_WIDTH) {
        let lo = chunk * LINE_WIDTH;
        let hi = (lo + LINE_WIDTH).min(q.len());
        let q_seq = &q[lo..hi];
        let s_seq = &s[lo..hi];

        let insertions = q_seq.iter().filter(|&&b| seq::is_gap(b)).count() as i64;
        let deletions = s_seq.iter().filter(|&&b| seq::is_gap(b)).count() as i64;
        let len = q_seq.len() as i64;

        if chunk > 0 {
            // A block advances the coordinate only if it contained a real base.
            let q_incr = i64::from(len > insertions);
            let s_incr = i64::from(s_seq.len() as i64 > deletions);
            if minus {
                if mode == AlignmentMode::WithSubjSeq {
                    q_start = q_end - q_incr;
                    s_start = s_end + s_incr;
                } else {
                    q_start = q_end + q_incr;
                    s_start = s_end - s_incr;
                }
            } else {
                q_start = q_end + q_incr;
                s_start = s_end + s_incr;
            }
        }

        // ── query line ────────────────────────────────────────────────────
        if flip {
            q_end = q_start - len + 1 + insertions;
            out.push_str("C ");
        } else {
            q_end = q_start + len - 1 - insertions;
            out.push_str("  ");
        }
        if len == insertions {
            q_end = q_start;
        }
        write_seq_line(out, &q_name, q_start, q_seq, q_end);

        // ── mutation line ─────────────────────────────────────────────────
        out.push_str(&" ".repeat(27));
        for j in 0..q_seq.len() {
            let (qc, sc) = (q_seq[j], s_seq[j]);
            if qc.eq_ignore_ascii_case(&sc) {
                out.push(' ');
                continue;
            }
            match mut_char(qc, sc) {
                Some(mc) => {
                    out.push(mc);
                    match mc {
                        'i' => t.transitions += 1,
                        'v' => t.transversions += 1,
                        '-' => {
                            t.gap_columns += 1;
                            // Faithful to the Perl's `$qChars[($j - 1) | 0]`:
                            // at j == 0 that is index -1, which in Perl is the
                            // *last* character of this block, not the previous
                            // one.  Reproduced so gap counts match exactly.
                            let prev = if j == 0 {
                                *q_seq.last().unwrap()
                            } else {
                                q_seq[j - 1]
                            };
                            if !seq::is_gap(prev) {
                                t.gap_openings += 1;
                            }
                        }
                        _ => {}
                    }
                }
                None if is_iupac_ambiguity(qc) || is_iupac_ambiguity(sc) => {
                    out.push('?');
                    t.ambiguous += 1;
                }
                None => out.push(' '),
            }
        }
        out.push('\n');

        // ── subject line ──────────────────────────────────────────────────
        if minus && mode == AlignmentMode::WithQuerySeq {
            out.push_str("C ");
            s_end = s_start - s_seq.len() as i64 + 1 + deletions;
        } else {
            out.push_str("  ");
            s_end = s_start + s_seq.len() as i64 - 1 - deletions;
        }
        write_seq_line(out, &s_name, s_start, s_seq, s_end);
        out.push('\n');
    }

    write_footer(out, r, &t);
    Ok(())
}

fn write_seq_line(out: &mut String, name: &str, start: i64, seq_bytes: &[u8], end: i64) {
    let _ = write!(out, "{name:<13} ");
    let _ = writeln!(
        out,
        "{:>10} {} {}",
        start,
        String::from_utf8_lossy(seq_bytes),
        end
    );
}

fn write_footer(out: &mut String, r: &SearchResult, t: &Tally) {
    let _ = writeln!(
        out,
        "Matrix = {}",
        r.matrix_name.as_deref().unwrap_or("Unknown")
    );
    if let Some(k) = r.kimura_divergence {
        let _ = writeln!(out, "Kimura (with divCpGMod) = {k}");
    }
    if let (Some(c), Some(k)) = (r.cpg_sites, r.raw_kimura_divergence) {
        let _ = writeln!(out, "CpG sites = {c}, Kimura (unadjusted) = {k}");
    }

    out.push_str("Transitions / transversions = ");
    if t.transversions > 0 {
        let _ = write!(
            out,
            "{:.2}",
            t.transitions as f64 / t.transversions as f64
        );
    } else {
        out.push_str("1.00");
    }
    let _ = writeln!(out, " ({}/{})", t.transitions, t.transversions);

    let a = &r.alignment;
    let (q_b, q_e) = a.query_one_based();
    // The Perl divides by getQueryEnd() - getQueryStart(), without the +1.
    let span = q_e as i64 - q_b as i64;
    if span > 0 {
        let _ = write!(
            out,
            "Gap_init rate = {:.2} ({} / {})",
            t.gap_openings as f64 / span as f64,
            t.gap_openings,
            span
        );
    } else {
        let _ = write!(out, "Gap_init rate = 0.0 ( {} / 0 )", t.gap_openings);
    }
    if t.gap_openings > 0 {
        let _ = write!(
            out,
            ", avg. gap size = {:.2} ({} / {})\n\n",
            t.gap_columns as f64 / t.gap_openings as f64,
            t.gap_columns,
            t.gap_openings
        );
    } else {
        out.push_str(", avg. gap size = 0.0 (0 / 0)\n\n");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::Alignment;
    use crate::Strand;

    fn plus_hit() -> (SearchResult, Vec<u8>, Vec<u8>) {
        let query = b"ACGTACGTACGTACGTACGT".to_vec();
        let subject = b"ACGTACGAACGTACGTACGT".to_vec();
        let a = Alignment::from_gapped(
            "chr1", "AluY", 0, 0, Strand::Plus, 500, &query, &subject,
        )
        .unwrap();
        let mut a = a;
        a.query_len = Some(query.len());
        a.subj_len = Some(subject.len());
        let mut r = SearchResult::new(a);
        r.pct_diverge = 5.0;
        r.subj_class = Some("SINE/Alu".into());
        r.id = Some(1);
        (r, query, subject)
    }

    #[test]
    fn out_line_plus_strand_puts_left_last() {
        let (r, _, _) = plus_hit();
        let line = to_out_line(&r);
        assert!(line.contains(" + "), "{line}");
        // subject: begin end (left)
        assert!(line.contains("      1      20     (0)"), "{line}");
    }

    #[test]
    fn out_line_minus_strand_moves_left_to_the_front_and_descends() {
        let (mut r, _, _) = plus_hit();
        r.alignment.strand = Strand::Minus;
        r.alignment.subj_len = Some(30);
        let line = to_out_line(&r);
        assert!(line.contains(" C "), "{line}");
        // subject: (left) end begin  — descending
        let cols: Vec<&str> = line.split_whitespace().collect();
        let left_idx = cols.iter().position(|c| *c == "(10)").expect("no (left)");
        assert_eq!(cols[left_idx + 1], "20", "end should follow (left)");
        assert_eq!(cols[left_idx + 2], "1", "begin should follow end");
    }

    #[test]
    fn annotation_line_truncates_names_at_whitespace() {
        let (mut r, _, _) = plus_hit();
        r.alignment.query_name = "chr1 some description".into();
        let line = annotation_line(&r);
        assert!(line.contains(" chr1 "), "{line}");
        assert!(!line.contains("description"), "{line}");
    }

    #[test]
    fn no_align_mode_emits_only_the_annotation() {
        let (r, q, s) = plus_hit();
        let out = to_crossmatch(&r, &q, &s, AlignmentMode::NoAlign).unwrap();
        assert_eq!(out.lines().count(), 1);
    }

    #[test]
    fn alignment_block_marks_transitions_and_transversions() {
        let (r, q, s) = plus_hit();
        let out = to_crossmatch(&r, &q, &s, AlignmentMode::WithQuerySeq).unwrap();
        // Query T vs subject A at column 7 is a transversion.
        assert!(out.contains('v'), "expected a transversion marker:\n{out}");
        assert!(out.contains("Transitions / transversions = "), "{out}");
        assert!(out.contains("Matrix = Unknown"), "{out}");
    }

    #[test]
    fn sequence_lines_use_the_repeatmasker_column_layout() {
        let (r, q, s) = plus_hit();
        let out = to_crossmatch(&r, &q, &s, AlignmentMode::WithQuerySeq).unwrap();
        let qline = out
            .lines()
            .find(|l| l.starts_with("  chr1"))
            .expect("no query line");
        // Two-space marker, 13-char name field, space, then a 10-wide position.
        assert_eq!(&qline[..2], "  ");
        assert_eq!(&qline[2..15], "chr1         ");
        assert_eq!(&qline[15..16], " ");
        assert_eq!(&qline[16..26], "         1");
    }

    #[test]
    fn long_alignments_wrap_at_fifty_columns() {
        let query: Vec<u8> = b"ACGT".iter().cycle().take(120).copied().collect();
        let a = Alignment::from_gapped(
            "q", "s", 0, 0, Strand::Plus, 1000, &query, &query,
        )
        .unwrap();
        let mut a = a;
        a.query_len = Some(query.len());
        a.subj_len = Some(query.len());
        let r = SearchResult::new(a);
        let out = to_crossmatch(&r, &query, &query, AlignmentMode::WithQuerySeq).unwrap();

        let q_lines: Vec<&str> = out.lines().filter(|l| l.starts_with("  q ")).collect();
        assert_eq!(q_lines.len(), 3, "120 columns should wrap into 3 blocks");
        // First block covers 1..50, second 51..100, third 101..120.
        assert!(q_lines[0].contains(" 1 "), "{}", q_lines[0]);
        assert!(q_lines[0].trim_end().ends_with(" 50"), "{}", q_lines[0]);
        assert!(q_lines[2].trim_end().ends_with(" 120"), "{}", q_lines[2]);
    }

    #[test]
    fn minus_strand_subject_line_carries_the_c_marker_and_descends() {
        let query = b"ACGTACGTACGTACGT".to_vec();
        // Subject whose reverse complement equals the query.
        let subject = seq::revcomp(&query);
        let mut a = Alignment::from_gapped(
            "q", "s", 0, 0, Strand::Minus, 400, &query, &query,
        )
        .unwrap();
        a.query_len = Some(query.len());
        a.subj_len = Some(subject.len());
        let r = SearchResult::new(a);

        let out = to_crossmatch(&r, &query, &subject, AlignmentMode::WithQuerySeq).unwrap();
        let s_line = out
            .lines()
            .find(|l| l.starts_with("C s"))
            .expect("subject line should carry the C marker");
        let nums: Vec<i64> = s_line
            .split_whitespace()
            .filter_map(|t| t.parse().ok())
            .collect();
        assert!(
            nums[0] > nums[1],
            "minus-strand subject coordinates should descend: {s_line}"
        );
    }

    #[test]
    fn with_subj_seq_mode_moves_the_marker_to_the_query() {
        let query = b"ACGTACGTACGTACGT".to_vec();
        let subject = seq::revcomp(&query);
        let mut a = Alignment::from_gapped(
            "q", "s", 0, 0, Strand::Minus, 400, &query, &query,
        )
        .unwrap();
        a.query_len = Some(query.len());
        a.subj_len = Some(subject.len());
        let r = SearchResult::new(a);

        let out = to_crossmatch(&r, &query, &subject, AlignmentMode::WithSubjSeq).unwrap();
        assert!(
            out.lines().any(|l| l.starts_with("C q")),
            "query line should carry the C marker in WithSubjSeq mode:\n{out}"
        );
        assert!(
            !out.lines().any(|l| l.starts_with("C s")),
            "subject line should not:\n{out}"
        );
    }

    #[test]
    fn gaps_are_marked_and_counted() {
        let gq = b"ACGT--ACGT".to_vec();
        let gs = b"ACGTTTACGT".to_vec();
        let query = b"ACGTACGT".to_vec();
        let subject = b"ACGTTTACGT".to_vec();
        let mut a =
            Alignment::from_gapped("q", "s", 0, 0, Strand::Plus, 200, &gq, &gs).unwrap();
        a.query_len = Some(query.len());
        a.subj_len = Some(subject.len());
        let r = SearchResult::new(a);

        let out = to_crossmatch(&r, &query, &subject, AlignmentMode::WithQuerySeq).unwrap();
        assert!(out.contains("--"), "gap should appear in the query row:\n{out}");
        assert!(out.contains("avg. gap size = "), "{out}");
    }

    #[test]
    fn footer_reports_kimura_when_present() {
        let (mut r, q, s) = plus_hit();
        r.matrix_name = Some("14p35g".into());
        r.kimura_divergence = Some(4.25);
        r.raw_kimura_divergence = Some(5.5);
        r.cpg_sites = Some(3);
        let out = to_crossmatch(&r, &q, &s, AlignmentMode::WithQuerySeq).unwrap();
        assert!(out.contains("Matrix = 14p35g"), "{out}");
        assert!(out.contains("Kimura (with divCpGMod) = 4.25"), "{out}");
        assert!(out.contains("CpG sites = 3, Kimura (unadjusted) = 5.5"), "{out}");
    }

    #[test]
    fn mut_char_is_keyed_query_then_subject() {
        assert_eq!(mut_char(b'C', b'T'), Some('i'));
        assert_eq!(mut_char(b'A', b'C'), Some('v'));
        assert_eq!(mut_char(b'A', b'-'), Some('-'));
        assert_eq!(mut_char(b'-', b'A'), Some('-'));
        // Ambiguity codes are not in the table; the caller emits '?'.
        assert_eq!(mut_char(b'A', b'N'), None);
    }
}
