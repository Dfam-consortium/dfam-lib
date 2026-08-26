# Third-party code redistributed in this repository

Everything written for `dfam-lib` is CC0-1.0 (see `LICENSE`). Some crates
additionally redistribute upstream source, which carries its own terms. This
file is the index of that; each entry maps a directory to exactly one upstream
project and one licence.

**The rule:** one vendored upstream = one directory under
`<crate>/third_party/<name>/` = one licence file inside it. Do not flatten a
second upstream into an existing directory — a scanner (and a downstream
packager) has to be able to attribute every file by its path alone.

## Inventory

| path | upstream | version | licence | licence file |
|---|---|---|---|---|
| `aln-parasail/third_party/parasail/` | [parasail](https://github.com/jeffdaily/parasail) | 2.6.2 | BSD-3-Clause (see note) | `LICENSE.md` |

Nothing else in the workspace vendors third-party source. `aln-rmblastn` shells
out to an external `rmblastn` binary but redistributes none of it; `aln-rmblast`
takes a git dependency on [RMBlast](https://github.com/Dfam-consortium/RMBlast)
(CC0-1.0) and likewise vendors nothing: cargo fetches it at build time, and the
ALP 1.98 sources RMBlast vendors for its `alp-fit` feature are attributed in
that repository's own `LICENSE`.

## Notes

### parasail

Copyright (c) 2015, Battelle Memorial Institute.

`aln-parasail/Cargo.toml` declares `license = "BSD-3-Clause AND CC0-1.0"`,
because Cargo accepts only SPDX identifiers and BSD-3-Clause is the closest one
— it is also how upstream describes itself. The actual grant is *not* verbatim
BSD-3-Clause. Two additions matter:

- **Name use.** "Other than as used herein, neither the name Battelle Memorial
  Institute or Battelle may be used in any form whatsoever without the express
  written consent of Battelle." This is broader than the usual BSD
  no-endorsement clause.
- **Citation.** "Redistributions of the software in any form, and publications
  based on work performed using the software should include the following
  citation as a reference: Daily, Jeff. (2016). Parasail: SIMD C library for
  global, semi-global, and local pairwise sequence alignments. *BMC
  Bioinformatics*, 17(1), 1-11. doi:10.1186/s12859-016-0930-z" — a "should", so
  a request rather than a condition, but cite it anyway in anything published
  from output this backend produced.

The binary-redistribution clause is why `LICENSE.md` must ship inside
`aln-parasail/third_party/parasail/` and be reproduced in the documentation of
anything that links this crate.

Only a curated subset of parasail is vendored (~30 of 595 generated `.c` files);
`aln-parasail/third_party/parasail/README.md` records exactly what and why.

### Adding an upstream

1. `mkdir aln-<crate>/third_party/<name>/` and copy the source in, with its
   unmodified licence file.
2. Write a `README.md` beside it: upstream URL, version, what was taken, what
   was left out, and how to re-vendor.
3. Extend the inventory table above, and any note the licence needs.
4. Update the crate's `license` field to the full SPDX expression — e.g. adding
   an MIT-licensed dependency makes `aln-parasail` `"CC0-1.0 AND BSD-3-Clause
   AND MIT"`.

The ARM/NEON work will exercise this: if parasail's NEON path needs a
compatibility shim (SIMDe and sais are both MIT), it lands as its own directory
alongside `parasail/`, not inside it.
