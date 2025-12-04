#!/usr/bin/env bash
set -euo pipefail

# gtf_to_tss.sh - extract TSS (gene-level) from a GTF (or GTF.gz) to a sorted BED
#
# Usage:
#   ./gtf_to_tss.sh INPUT.gtf [OUTPUT.bed]
#
# Notes:
# - Expects gene features in column 3 ($3 == "gene").
# - BED fields: # Chromosome, Start, End, Gene, Score, Strand, Transcript_type
# - Output is sorted by chrom then start.

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 INPUT.gtf[.gz] [OUTPUT.bed]" >&2
  exit 1
fi

INPUT="$1"
if [[ $# -eq 2 ]]; then
  OUTPUT="$2"
else
  base="$(basename "${INPUT%.*}")"
  OUTPUT="${base}_tss.bed"
fi

# check input exists
if [[ ! -e "$INPUT" ]]; then
  echo "ERROR: input file not found: $INPUT" >&2
  exit 2
fi

# choose decompressor (support gzipped GTF)
if [[ "$INPUT" == *.gz ]]; then
  DECOMP="zcat --silent"
else
  DECOMP="cat"
fi

# temp file
tmp="$(mktemp --suffix=.tss.tmp)"
trap 'rm -f "$tmp"' EXIT

# extract gene TSS lines and write unsorted tmp
$DECOMP "$INPUT" | awk 'BEGIN{FS="\t";OFS="\t"}
  $0 ~ /^#/ { next }           # skip headers
  $3 == "gene" {
    # determine TSS depending on strand
    if ($7 == "+") tss = $4 + 0
    else tss = $5 + 0

    start = tss - 1
    end = tss

    name = "."
    if (match($0, /gene_name "([^"]+)"/, a)) {
      name = a[1]
    } else if (match($0, /gene_id "([^"]+)"/, b)) {
      name = b[1]
    }

    # Add dummy Transcript_type column as NA
    transcript_type = "NA"

    # print headered line for Polars-compatible BED
    print $1, start, end, name, 0, $7, transcript_type
  }' > "$tmp"

# sort and uniq (keep unique lines); write to final output
{
  # print standard header
  echo -e "# Chromosome\tStart\tEnd\tGene\tScore\tStrand\tTranscript_type"
  sort -k1,1 -k2,2n "$tmp" | uniq
} > "$OUTPUT"

echo "Wrote TSS BED to: $OUTPUT"

