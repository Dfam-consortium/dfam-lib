//! Minimal sequence I/O.
//!
//! Only what the tools need: FASTA in and out, plus the GIRI `.ig` header form
//! `autocons` writes by default. Anything richer (2bit, Stockholm, indexed
//! databases) belongs in the crate that needs it — `aln-core` stays free of
//! heavyweight dependencies.

use std::io::{BufRead, Write};

use crate::error::{Error, Result};
use crate::seq::Sequence;

/// Read every record from a FASTA stream.
///
/// Blank lines are skipped and sequence lines are concatenated with whitespace
/// removed. The record name is the first whitespace-delimited token of the
/// header; the remainder is discarded — callers that need the full defline
/// should use [`read_fasta_with_deflines`].
pub fn read_fasta<R: BufRead>(reader: R) -> Result<Vec<Sequence>> {
    Ok(read_fasta_with_deflines(reader)?
        .into_iter()
        .map(|(s, _)| s)
        .collect())
}

/// Like [`read_fasta`], but also returns each record's full header line.
pub fn read_fasta_with_deflines<R: BufRead>(reader: R) -> Result<Vec<(Sequence, String)>> {
    let mut out: Vec<(Sequence, String)> = Vec::new();
    let mut current: Option<(String, String, Vec<u8>)> = None;

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(header) = trimmed.strip_prefix('>') {
            if let Some((name, defline, seq)) = current.take() {
                out.push((Sequence::new(name, seq), defline));
            }
            let defline = header.to_string();
            let name = header.split_whitespace().next().unwrap_or("").to_string();
            current = Some((name, defline, Vec::new()));
        } else {
            match current.as_mut() {
                Some((_, _, seq)) => {
                    seq.extend(trimmed.bytes().filter(|b| !b.is_ascii_whitespace()));
                }
                None => {
                    return Err(Error::Io(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        "FASTA sequence data before the first '>' header",
                    )))
                }
            }
        }
    }
    if let Some((name, defline, seq)) = current.take() {
        out.push((Sequence::new(name, seq), defline));
    }
    Ok(out)
}

/// Read FASTA from a file path.
pub fn read_fasta_file(path: impl AsRef<std::path::Path>) -> Result<Vec<Sequence>> {
    let f = std::fs::File::open(path.as_ref())?;
    read_fasta(std::io::BufReader::new(f))
}

/// Write one FASTA record, wrapping at `width` columns (0 disables wrapping).
///
/// `comment` is appended to the header after a space, which is how `autocons`
/// carries its `SCORE=` annotation.
pub fn write_fasta<W: Write>(
    w: &mut W,
    name: &str,
    comment: Option<&str>,
    seq: &[u8],
    width: usize,
) -> Result<()> {
    match comment {
        Some(c) => writeln!(w, ">{name} {c}")?,
        None => writeln!(w, ">{name}")?,
    }
    write_body(w, seq, width)
}

/// Write one record in GIRI `.ig` form: two `;` comment lines, the name on its
/// own line, then the sequence.
///
/// This is `autocons`'s default output when `--fa` is not given.
pub fn write_ig<W: Write>(
    w: &mut W,
    name: &str,
    comment: Option<&str>,
    seq: &[u8],
    width: usize,
) -> Result<()> {
    // Matched byte-for-byte against GIRI's autocons: a space after the
    // semicolon, a second bare `; ` line, and NO `1`/`2` terminator -- strict
    // IG has one, GIRI does not emit it.
    writeln!(w, "; {}", comment.unwrap_or(""))?;
    writeln!(w, "; ")?;
    writeln!(w, "{name}")?;
    write_body(w, seq, width)?;
    Ok(())
}

fn write_body<W: Write>(w: &mut W, seq: &[u8], width: usize) -> Result<()> {
    if width == 0 || seq.len() <= width {
        w.write_all(seq)?;
        writeln!(w)?;
        return Ok(());
    }
    for chunk in seq.chunks(width) {
        w.write_all(chunk)?;
        writeln!(w)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_multiple_records() {
        let text = ">one first record\nACGT\nACGT\n>two\nTTTT\n";
        let seqs = read_fasta(text.as_bytes()).unwrap();
        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[0].name, "one");
        assert_eq!(seqs[0].seq, b"ACGTACGT");
        assert_eq!(seqs[1].name, "two");
        assert_eq!(seqs[1].seq, b"TTTT");
    }

    #[test]
    fn keeps_the_full_defline_when_asked() {
        let text = ">one first record\nACGT\n";
        let recs = read_fasta_with_deflines(text.as_bytes()).unwrap();
        assert_eq!(recs[0].0.name, "one");
        assert_eq!(recs[0].1, "one first record");
    }

    #[test]
    fn blank_lines_and_internal_whitespace_are_ignored() {
        let text = ">a\n\nAC GT\n\nAC\tGT\n";
        let seqs = read_fasta(text.as_bytes()).unwrap();
        assert_eq!(seqs[0].seq, b"ACGTACGT");
    }

    #[test]
    fn data_before_a_header_is_an_error() {
        assert!(read_fasta("ACGT\n>a\nACGT\n".as_bytes()).is_err());
    }

    #[test]
    fn empty_input_yields_no_records() {
        assert!(read_fasta("".as_bytes()).unwrap().is_empty());
    }

    #[test]
    fn fasta_round_trips() {
        let mut buf = Vec::new();
        write_fasta(&mut buf, "x", None, b"ACGTACGTAC", 4).unwrap();
        assert_eq!(String::from_utf8(buf.clone()).unwrap(), ">x\nACGT\nACGT\nAC\n");
        let back = read_fasta(buf.as_slice()).unwrap();
        assert_eq!(back[0].seq, b"ACGTACGTAC");
    }

    #[test]
    fn a_comment_lands_on_the_header() {
        let mut buf = Vec::new();
        write_fasta(&mut buf, "cons", Some("SCORE=12.50"), b"ACGT", 0).unwrap();
        assert_eq!(String::from_utf8(buf).unwrap(), ">cons SCORE=12.50\nACGT\n");
    }

    #[test]
    fn ig_format_matches_giri_autocons() {
        // Pinned against the real thing: `bin/autocons fam.fa --orig --sse2`
        // emits "; <comment>", a bare "; ", the name, then the sequence — with
        // no `1`/`2` terminator, even though strict IG has one.  autocons is a
        // drop-in replacement, so it follows GIRI here rather than the spec.
        let mut buf = Vec::new();
        write_ig(&mut buf, "cons", Some("SCORE=1.00"), b"ACGT", 0).unwrap();
        assert_eq!(
            String::from_utf8(buf).unwrap(),
            "; SCORE=1.00\n; \ncons\nACGT\n"
        );
    }
}
