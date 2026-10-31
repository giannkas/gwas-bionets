#!/usr/bin/env bash

# -----------------------------------------------------------------------------
# Psoriasis MAGMA workflow
#
# This script:
#   1. Keeps SNP associations with 0 < P <= 5e-8 and sorts them by P-value.
#   2. Converts the significant PLINK association file to MAGMA SNP-location
#      input format (SNP, CHR, BP, P).
#   3. Runs MAGMA gene annotation using a 50 kb window.
#   4. Extracts the annotated gene IDs, maps them to gene names, and writes a
#      file containing only the gene names.
#
# Author: Giann Karlo Aguirre-Samboní
#
# Run this script from the project root directory. You can redefine all input 
# and output variables as needed; the values shown here are for illustrative purposes only.
#
# Required commands: bash, awk, sort, cut, Rscript, magma
# Required R packages: readr, dplyr, gprofiler2
# -----------------------------------------------------------------------------

set -Eeuo pipefail

trap 'echo "ERROR: command failed on line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

INPUT_DIR="results/psoriasis/magma_scores/magma_output_50kb_window_snp-array"
INPUT_ASSOC="${INPUT_DIR}/snp_association_plink.assoc"

# The output directory is the input directory name with this suffix appended.
OUTPUT_DIR="${INPUT_DIR}_significant_snps"

# MAGMA gene-coordinate file (GRCh38 / NCBI38).
MAGMA_GENE_LOC="${MAGMA_GENE_LOC:-magma/genome_data/NCBI38/NCBI38.gene.loc}"

# MAGMA annotation window
WINDOW=50

# P-value threshold
PVALUE=5e-8

# Output files
SIGNIFICANT_ASSOC="${OUTPUT_DIR}/snp_association_plink.assoc"
MAGMA_SNP_INPUT="${OUTPUT_DIR}/snp_association_plink.tsv"
MAGMA_OUTPUT_PREFIX="${OUTPUT_DIR}/pso"
GENES_ANNOT="${MAGMA_OUTPUT_PREFIX}.genes.annot"
GENE_IDS_FILE="${OUTPUT_DIR}/pso.genes.annot_only_gene_ids.tsv"
GENE_NAMES_WITH_IDS="${OUTPUT_DIR}/pso.genes.annot_gene_ids_names.tsv"
GENE_NAMES_FILE="${OUTPUT_DIR}/pso.genes.annot_only_gene_names.tsv"

# -----------------------------------------------------------------------------
# Checks
# -----------------------------------------------------------------------------

for command in awk sort cut Rscript magma; do
    if ! command -v "$command" >/dev/null 2>&1; then
        echo "ERROR: required command not found: $command" >&2
        exit 1
    fi
done

if [[ ! -f "$INPUT_ASSOC" ]]; then
    echo "ERROR: PLINK association file not found: $INPUT_ASSOC" >&2
    exit 1
fi

if [[ ! -f "$MAGMA_GENE_LOC" ]]; then
    echo "ERROR: MAGMA gene-location file not found: $MAGMA_GENE_LOC" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# -----------------------------------------------------------------------------
# 1. Filter significant SNP associations
# -----------------------------------------------------------------------------
# Keep the original header, retain rows with 0 < P <= 5e-8 (P is column 9),
# and sort the retained rows by ascending P-value. You can modify the threshold
# in the code below.

echo "[1/6] Filtering significant SNP associations..."
{
    head -n 1 "$INPUT_ASSOC"
    tail -n +2 "$INPUT_ASSOC" \
        | LC_ALL=C awk -v pval="$PVALUE" '$9 > 0 && $9 <= pval' \
        | LC_ALL=C sort -g -k9,9
} > "$SIGNIFICANT_ASSOC"

if [[ ! -s "$SIGNIFICANT_ASSOC" ]]; then
    echo "WARNING: No significant associations found (empty file)" >&2
fi

# -----------------------------------------------------------------------------
# 2. Convert the significant association file to MAGMA SNP-location format
# -----------------------------------------------------------------------------
# MAGMA input is written as a tab-separated file with the columns:
# SNP, CHR, BP, P.

echo "[2/6] Formatting significant associations for MAGMA..."
Rscript - "$SIGNIFICANT_ASSOC" "$MAGMA_SNP_INPUT" <<'RSCRIPT'
suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
input_file  <- args[[1]]
output_file <- args[[2]]

required_columns <- c("SNP", "CHR", "BP", "P")
association_data <- readr::read_table(
  input_file,
  show_col_types = FALSE
)

missing_columns <- setdiff(required_columns, names(association_data))
if (length(missing_columns) > 0) {
  stop(
    "Missing required column(s) in the PLINK association file: ",
    paste(missing_columns, collapse = ", ")
  )
}

association_data %>%
  dplyr::select(dplyr::all_of(required_columns)) %>%
  readr::write_tsv(output_file)
RSCRIPT

# -----------------------------------------------------------------------------
# 3. Run MAGMA gene annotation
# -----------------------------------------------------------------------------
# Annotate SNPs to genes using a 50 kb window around each gene.

echo "[3/6] Running MAGMA gene annotation..."
magma \
    --annotate window="$WINDOW" \
    --snp-loc "$MAGMA_SNP_INPUT" \
    --gene-loc "$MAGMA_GENE_LOC" \
    --out "$MAGMA_OUTPUT_PREFIX"

# -----------------------------------------------------------------------------
# 4. Extract gene IDs from the MAGMA annotation output
# -----------------------------------------------------------------------------
# Keep the first field from each data row, skip the file's first line, and
# replace any leading '#' marker (MAGMA format) with the expected GENE column label.

echo "[4/6] Extracting gene IDs from the MAGMA annotation..."
awk 'NR > 1 { gsub(/#/, "GENE"); print $1 }' "$GENES_ANNOT" > "$GENE_IDS_FILE"

# -----------------------------------------------------------------------------
# 5. Convert gene IDs to gene names
# -----------------------------------------------------------------------------
# gconvert maps Entrez Gene accession IDs to gene names. The result is joined
# back to the extracted MAGMA gene IDs and duplicate IDs are removed, keeping
# the first occurrence.

echo "[5/6] Mapping gene IDs to gene names..."
Rscript - "$GENE_IDS_FILE" "$GENE_NAMES_WITH_IDS" <<'RSCRIPT' || exit 1
suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
input_file  <- args[[1]]
output_file <- args[[2]]

magma_genes <- readr::read_table(
  input_file,
  show_col_types = FALSE
)

if (!"GENE" %in% names(magma_genes)) {
  stop("The extracted gene-ID file does not contain a GENE column.")
}

converted_genes <- gprofiler2::gconvert(
  magma_genes$GENE,
  numeric_ns = "ENTREZGENE_ACC",
  mthreshold = 1,
  filter_na = TRUE
)

gene_names <- converted_genes %>%
  dplyr::mutate(GENE = as.numeric(input)) %>%
  dplyr::select(GENE, name) %>%
  dplyr::left_join(magma_genes, by = "GENE") %>%
  dplyr::rename(Gene = name)

readr::write_tsv(gene_names, output_file)
RSCRIPT

if [[ ! -s "$GENE_NAMES_WITH_IDS" ]]; then
    echo "WARNING: empty file after accession conversion." >&2
fi

# Remove duplicate gene IDs while preserving their first occurrence.
TEMP_GENE_NAMES="$(mktemp "${GENE_NAMES_WITH_IDS}.XXXXXX")"
awk -F '\t' '!seen[$1]++' "$GENE_NAMES_WITH_IDS" > "$TEMP_GENE_NAMES" || exit 1
mv "$TEMP_GENE_NAMES" "$GENE_NAMES_WITH_IDS" || exit 1

# -----------------------------------------------------------------------------
# 6. Write a file containing only gene names
# -----------------------------------------------------------------------------
# The second column contains the mapped gene names; the header is retained.

echo "[6/6] Writing the gene-name-only file..."
cut -f 2 "$GENE_NAMES_WITH_IDS" > "$GENE_NAMES_FILE"

# -----------------------------------------------------------------------------
# Completion summary
# -----------------------------------------------------------------------------

echo
echo "Significant genes workflow completed successfully."
echo "Output directory:            $OUTPUT_DIR"
echo "Significant associations:    $SIGNIFICANT_ASSOC"
echo "MAGMA SNP input:              $MAGMA_SNP_INPUT"
echo "MAGMA annotation:             $GENES_ANNOT"
echo "Gene IDs:                     $GENE_IDS_FILE"
echo "Gene IDs with gene names:     $GENE_NAMES_WITH_IDS"
echo "Gene names only:              $GENE_NAMES_FILE"
