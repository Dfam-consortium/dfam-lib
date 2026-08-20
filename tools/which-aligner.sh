#!/bin/sh
# which-aligner.sh <path-to-autocons> [test.fa]
#
# Reports whether an autocons binary defaults to Farrar (SSE2) or Monardo.
#
# Why behavioural and not `nm`/`strings`: both aligner classes are linked into
# every build (they are reachable via --mon / --sse2), so symbol presence tells
# you nothing. But the *default* must match exactly one of them, and Monardo and
# Farrar produce different alignments on any non-trivial input.
set -e
BIN=${1:?usage: which-aligner.sh <autocons> [test.fa]}
FA=${2:-$(dirname "$0")/aligner-probe.fa}
[ -r "$FA" ] || { echo "no test FASTA at $FA" >&2; exit 2; }

run() { "$BIN" "$FA" --fa "$@" 2>/dev/null | grep -v '^>' | tr -d '\n'; }
def=$(run); mon=$(run --mon); sse=$(run --sse2)

[ -n "$def" ] || { echo "INCONCLUSIVE: binary produced no output" >&2; exit 2; }
if [ "$mon" = "$sse" ]; then
  echo "INCONCLUSIVE: --mon and --sse2 agree on this input; use a more divergent family" >&2
  exit 2
fi
if [ "$def" = "$sse" ]; then
  echo "Farrar (SSE2)  <- intended default; HAVE_CONFIG_H was defined at build time"
elif [ "$def" = "$mon" ]; then
  echo "Monardo        <- BUILD REGRESSION: config.h not included, so HAVE_SSE2 undefined"
  echo "                  rebuild with -DHAVE_CONFIG_H (see build-acons.sh)"
  exit 1
else
  echo "UNKNOWN: default matches neither --mon nor --sse2" >&2
  exit 2
fi
