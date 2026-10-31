#!/usr/bin/env bash
#
# Apply a Benjamini-Hochberg (BH) FDR threshold to MAGMA gene-level results.
#
# The script:
#   1. Reads a MAGMA results table (tab-separated, with a header).
#   2. Sorts the gene rows by their p-value (column 10) in ascending order.
#   3. Assigns each gene a rank and calculates its BH critical value:
#          BHvalue = (rank / number_of_genes) * FDR
#   4. Keeps rows whose p-value is strictly smaller than their BH critical value.
#   5. Writes the retained rows, with rank and BHvalue columns, to the output TSV.
#
# Author: Giann Karlo Aguirre-Samboní
#
# Usage:
#   ./run_magma_bh_correction.sh [input.tsv] [output.tsv] [FDR]
#
# Defaults:
#   input.tsv  = pso.scores.genes.out_converted.tsv
#   output.tsv = pso.scores.genes.out_converted_bhcorrection.tsv
#   FDR        = 0.05
#
# Notes:
#   - The input is expected to be tab-separated and to contain a header row.
#   - By default, the p-value is expected in column 10 of the input file.
#   - The denominator is the total number of data rows (excluding the header).
#   - The original command used a strict '<' comparison; this script preserves it.
#   - LC_ALL=C ensures that decimal points and scientific notation are parsed
#     consistently, making the comma-conversion workaround unnecessary.

set -Eeuo pipefail
export LC_ALL=C

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

INPUT_FILE="${1:-pso.scores.genes.out_converted.tsv}"
OUTPUT_FILE="${2:-pso.scores.genes.out_converted_bhcorrection.tsv}"
FDR="${3:-0.05}"
P_VALUE_COLUMN=10

# -----------------------------------------------------------------------------
# Checks
# -----------------------------------------------------------------------------

if [[ ! -f "$INPUT_FILE" ]]; then
    printf 'Error: input file not found: %s\n' "$INPUT_FILE" >&2
    exit 1
fi

if [[ ! -s "$INPUT_FILE" ]]; then
    printf 'Error: input file is empty: %s\n' "$INPUT_FILE" >&2
    exit 1
fi

# ensure the FDR value is numeric and lies in the interval (0, 1].
if ! awk -v fdr="$FDR" 'BEGIN {
    if (fdr ~ /^[+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/ &&
        fdr > 0 && fdr <= 1) {
        exit 0
    }
    exit 1
}'; then
    printf 'Error: FDR must be a number greater than 0 and no greater than 1 (received: %s).\n' "$FDR" >&2
    exit 1
fi

# check that the input has the expected p-value column and that every data-row p-value is numeric.
if ! awk -F '\t' -v pcol="$P_VALUE_COLUMN" '
    NR == 1 {
        if (NF < pcol) {
            printf "Error: the input has %d columns; p-values are expected in column %d.\n", NF, pcol > "/dev/stderr"
            exit 1
        }
        next
    }
    {
        p = $pcol
        if (p !~ /^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/) {
            printf "Error: non-numeric p-value on input line %d: %s\n", NR, p > "/dev/stderr"
            exit 1
        }
    }
' "$INPUT_FILE"; then
    exit 1
fi

# count genes (data rows), excluding the header.
GENE_COUNT="$(awk 'END { print NR - 1 }' "$INPUT_FILE")"

if (( GENE_COUNT < 1 )); then
    printf 'Error: no gene rows found in %s.\n' "$INPUT_FILE" >&2
    exit 1
fi

# create the output directory if the caller supplied a path containing one.
OUTPUT_DIR="$(dirname -- "$OUTPUT_FILE")"
mkdir -p "$OUTPUT_DIR"

# use temporary files so an interrupted run does not leave a partial result.
SORTED_TMP="$(mktemp)"
OUTPUT_TMP="$(mktemp "${OUTPUT_DIR}/.magma_bh_output.XXXXXX")"
trap 'rm -f "$SORTED_TMP" "$OUTPUT_TMP"' EXIT

# -----------------------------------------------------------------------------
# Sort the gene rows by p-value, preserving the header as the first line
# -----------------------------------------------------------------------------

TAB="$(printf '\t')"
{
    head -n 1 "$INPUT_FILE"
    tail -n +2 "$INPUT_FILE" |
        sort -t "$TAB" -g -k "${P_VALUE_COLUMN}"
} > "$SORTED_TMP"

# -----------------------------------------------------------------------------
# Calculate BH critical values and retain rows passing the threshold
# -----------------------------------------------------------------------------

awk -F '\t' \
    -v OFS='\t' \
    -v fdr="$FDR" \
    -v gene_count="$GENE_COUNT" \
    -v pcol="$P_VALUE_COLUMN" '
    NR == 1 {
        # add two new column names to the original header.
        print "rank", $0, "BHvalue"
        next
    }
    {
        rank = NR - 1
        p_value = $pcol + 0
        bh_value = (rank / gene_count) * fdr

        # strict inequality between p-values and bh-values.
        if (p_value < bh_value) {
            # print the BH threshold at high precision with scientific notation.
            printf "%d\t%s\t%e\n", rank, $0, bh_value
        }
    }
' "$SORTED_TMP" > "$OUTPUT_TMP"

# move the completed result into place only after the calculation succeeds.
mv -f "$OUTPUT_TMP" "$OUTPUT_FILE"

printf 'BH correction complete.\n'
printf '  Input:      %s\n' "$INPUT_FILE"
printf '  Number of genes: %s\n' "$GENE_COUNT"
printf '  FDR:        %s\n' "$FDR"
printf '  Output:     %s\n' "$OUTPUT_FILE"
