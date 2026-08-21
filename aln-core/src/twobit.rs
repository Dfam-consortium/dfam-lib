//! Minimal random-access reader for UCSC 2bit genome files.
//!
//! Consolidated here from three copies — dfam-curator, RepeatAfterMe, and this
//! crate — so that one implementation serves every consumer. This is
//! RepeatAfterMe's variant, which was the most careful of the three: it rejects
//! inverted ranges and short-circuits empty ones, where the others underflow
//! `end - start` on `u64` and panic in `Vec::with_capacity`.
//!
//! The index
//! (names → offsets, N-block tables) is built on open; `fetch` then reads only
//! the packed bytes covering the requested range via `pread`, so the struct is
//! `Send + Sync` and `fetch` takes `&self`. Returns uppercase ASCII `ACGTN` —
//! mask blocks are ignored, which matches the C loader's unconditional
//! `toUpperN`. Handles both endiannesses and versions 0 (32-bit offsets) and
//! 1 (64-bit). Unix-only (`FileExt::read_exact_at`).

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::os::unix::fs::FileExt;
use std::path::Path;

struct SeqEntry {
    dna_size: u64,
    n_starts: Vec<u32>,
    n_sizes: Vec<u32>,
    /// Absolute byte offset where the packed 2-bit DNA begins.
    packed_dna_offset: u64,
}

pub struct TwoBitReader {
    file: File,
    seqs: HashMap<String, SeqEntry>,
    /// Sequence names in file index order (the order kentsrc tools use).
    order: Vec<String>,
}

impl TwoBitReader {
    /// Open a 2bit file and build the index. Reads no packed DNA.
    pub fn open(path: &Path) -> io::Result<Self> {
        let mut f = File::open(path)?;

        let mut buf4 = [0u8; 4];
        f.read_exact(&mut buf4)?;
        let is_le = match u32::from_be_bytes(buf4) {
            0x1A412743 => false,
            0x4327411A => true,
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid 2bit signature",
                ))
            }
        };

        macro_rules! ru32 {
            () => {{
                let mut b = [0u8; 4];
                f.read_exact(&mut b)?;
                if is_le {
                    u32::from_le_bytes(b)
                } else {
                    u32::from_be_bytes(b)
                }
            }};
        }
        macro_rules! ru64 {
            () => {{
                let mut b = [0u8; 8];
                f.read_exact(&mut b)?;
                if is_le {
                    u64::from_le_bytes(b)
                } else {
                    u64::from_be_bytes(b)
                }
            }};
        }

        let version = ru32!();
        if version > 1 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "unsupported 2bit version",
            ));
        }
        let seq_count = ru32!() as usize;
        ru32!(); // reserved

        let mut index: Vec<(String, u64)> = Vec::with_capacity(seq_count);
        for _ in 0..seq_count {
            let mut len_buf = [0u8; 1];
            f.read_exact(&mut len_buf)?;
            let mut name_buf = vec![0u8; len_buf[0] as usize];
            f.read_exact(&mut name_buf)?;
            let name = String::from_utf8(name_buf)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            let offset = if version == 0 { ru32!() as u64 } else { ru64!() };
            index.push((name, offset));
        }

        let order: Vec<String> = index.iter().map(|(n, _)| n.clone()).collect();

        // Visit sequence records in record-offset order so header reads are
        // sequential.
        index.sort_unstable_by_key(|(_, off)| *off);

        let mut seqs = HashMap::with_capacity(seq_count);
        for (name, seq_off) in index {
            f.seek(SeekFrom::Start(seq_off))?;

            let dna_size = ru32!() as u64;
            let n_count = ru32!() as usize;

            let mut n_starts = Vec::with_capacity(n_count);
            let mut n_sizes = Vec::with_capacity(n_count);
            for _ in 0..n_count {
                n_starts.push(ru32!());
            }
            for _ in 0..n_count {
                n_sizes.push(ru32!());
            }

            let mask_count = ru32!() as usize;
            // Skip mask-block starts/sizes and the reserved word.
            f.seek(SeekFrom::Current((mask_count as i64) * 8 + 4))?;

            let packed_dna_offset = f.stream_position()?;
            seqs.insert(
                name,
                SeqEntry {
                    dna_size,
                    n_starts,
                    n_sizes,
                    packed_dna_offset,
                },
            );
        }

        Ok(TwoBitReader { file: f, seqs, order })
    }

    /// Fetch bases in the 0-based half-open range `[start, end)`.
    /// Returns uppercase ASCII `ACGTN`.
    pub fn fetch(&self, chrom: &str, start: u64, end: u64) -> io::Result<Vec<u8>> {
        let entry = self.seqs.get(chrom).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("chromosome {chrom:?} not found in 2bit file"),
            )
        })?;

        if end > entry.dna_size || start > end {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "range {start}..{end} exceeds {chrom} length {}",
                    entry.dna_size
                ),
            ));
        }
        if start == end {
            return Ok(Vec::new());
        }

        // 2-bit packing: T=0, C=1, A=2, G=3, MSB first within each byte.
        const BASES: [u8; 4] = [b'T', b'C', b'A', b'G'];

        let first_byte = start / 4;
        let last_byte = (end - 1) / 4;
        let byte_count = (last_byte - first_byte + 1) as usize;

        let mut packed = vec![0u8; byte_count];
        self.file
            .read_exact_at(&mut packed, entry.packed_dna_offset + first_byte)?;

        let length = (end - start) as usize;
        let mut result = Vec::with_capacity(length);
        for pos in start..end {
            let byte_idx = ((pos / 4) - first_byte) as usize;
            let shift = (3 - (pos % 4)) * 2;
            result.push(BASES[((packed[byte_idx] >> shift) & 0x03) as usize]);
        }

        // Overlay N-blocks.
        for (&ns, &nz) in entry.n_starts.iter().zip(entry.n_sizes.iter()) {
            let ns = ns as u64;
            let ne = ns + nz as u64;
            let ov_s = ns.max(start);
            let ov_e = ne.min(end);
            for i in ov_s..ov_e.max(ov_s) {
                result[(i - start) as usize] = b'N';
            }
        }

        Ok(result)
    }

    pub fn contains(&self, chrom: &str) -> bool {
        self.seqs.contains_key(chrom)
    }

    pub fn seq_len(&self, chrom: &str) -> Option<u64> {
        self.seqs.get(chrom).map(|e| e.dna_size)
    }

    /// Iterate over (name, length) pairs, unordered.
    /// Sequences in file index order (deterministic, matches kentsrc tools).
    pub fn sequences(&self) -> impl Iterator<Item = (&str, u64)> + '_ {
        self.order
            .iter()
            .map(|n| (n.as_str(), self.seqs[n].dna_size))
    }
}
