#!/bin/bash
# --- Purpose ---
# Confidence CLI: prepare known samples, assess support bins, and export confidence estimates.

#################################################################
# GeneParliamentID confidence estimation                            #
# Benedikt Kuhnhaeuser                                          #
# Royal Botanic Gardens, Kew                                    #
# 2026                                                          #
#################################################################

set -euo pipefail

# --- Installation paths, workflow defaults, and shared state ---
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SCRIPTS_DIR="$SCRIPT_DIR"
PROJECT_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
VERSION_FILE="${GPID_VERSION_FILE:-$PROJECT_DIR/VERSION}"
if [ -z "${GPID_VERSION:-}" ] && [ -f "$VERSION_FILE" ]; then
    IFS= read -r GPID_VERSION < "$VERSION_FILE" || GPID_VERSION=""
    GPID_VERSION=${GPID_VERSION%$'\r'}
fi
GPID_VERSION="${GPID_VERSION:-unknown}"
OUTPUT_DIR="confidence"
PREPARATIONS_DIR="$OUTPUT_DIR/preparations"
TESTS_DIR="$OUTPUT_DIR/tests"
DEFAULT_PREPARED_FILE="$PREPARATIONS_DIR/confidence_prepared.rds"
DEFAULT_TOP_IDS_FILE="$TESTS_DIR/confidence_top_ids.rds"

declare -A REFERENCE_GENE_FILES=()
declare -A CONFIDENCE_GENE_FILES=()

# usage(): Print command usage, supported options, and output locations.
usage() {
    printf 'GPID version: %s\n\n' "$GPID_VERSION"
    cat <<'EOF'
Usage: gpid confidence <command> [arguments]

Commands:
  prepare       Prepare BLAST and R input data for confidence
  estimate      Estimate confidence for different support bins
  bins          Save confidence support probabilities for selected bins
  help          Show this help message

Use gpid confidence <command> -h for command-specific help.

Examples:
  gpid confidence prepare -r reference -i confidence_samples
  gpid confidence estimate -g gene_performance.csv -t thresholds_filtering.csv
  gpid confidence bins -b 5
EOF
}

# usage_prepare(): Print confidence preparation inputs, grouping options, and outputs.
usage_prepare() {
    printf 'GPID version: %s\n\n' "$GPID_VERSION"
    cat <<'EOF'
Usage: gpid confidence prepare -r <reference directory> -i <confidence dataset directory> [-s <species groups file>]

Required:
  -r  Reference dataset directory containing one FASTA file per gene
      Missing BLAST databases are built automatically.
  -i  Confidence dataset directory containing one FASTA file per gene

Optional:
  -s  Species groups CSV with header genus_species,species_group
      If omitted, species groups are derived from genus names.
  -h  Show this help message

Outputs (paths relative to the working directory):
  confidence/preparations/confidence_blast.tsv
  confidence/preparations/confidence_prepared.rds
EOF
}

# usage_estimate(): Print confidence estimation inputs and output locations.
usage_estimate() {
    printf 'GPID version: %s\n\n' "$GPID_VERSION"
    cat <<'EOF'
Usage: gpid confidence estimate [-i <prepared confidence RDS>] -g <gene performance CSV> -t <filtering thresholds CSV>

Required:
  -g  Gene performance CSV produced by gpid calibrate genes
      Normally saved as: gene_performance.csv
      Unique gene and performance columns are required; extra columns are allowed.
      NA performance is treated as 0 for filtering, with a warning.
  -t  Filtering thresholds CSV from gpid calibrate combine (parameter,value or legacy single-row format)
      Normally saved as: thresholds_filtering.csv

Optional:
  -i  Intermediate RDS produced by gpid confidence prepare
      Default: confidence/preparations/confidence_prepared.rds
  -h  Show this help message

Outputs (paths relative to the working directory):
  confidence/tests/confidence_estimate.pdf
  confidence/tests/confidence_top_ids.rds
EOF
}

# usage_bins(): Print confidence bin options and explain CSV versus plot semantics.
usage_bins() {
    printf 'GPID version: %s\n\n' "$GPID_VERSION"
    cat <<'EOF'
Usage: gpid confidence bins -b <number of bins> [-i <confidence top IDs RDS>]

Required:
  -b  Number of confidence support bins (integer from 1 to 100)

Optional:
  -i  Top IDs RDS produced by gpid confidence estimate
      Default: confidence/tests/confidence_top_ids.rds
  -h  Show this help message

Outputs (paths relative to the working directory):
  confidence_support.csv
  confidence/confidence_support.pdf

In the CSV, probability_close includes correct and close identifications.
The plot shows correct, close and wrong as separate categories.
EOF
}

# --- Shared logging and input helpers ---
# log(): Write a progress or result message to standard output.
log() {
    printf '%s\n' "$1"
}

# warn(): Write a nonfatal warning to standard error.
warn() {
    printf 'Warning: %s\n' "$1" >&2
}

# die(): Report a fatal error and terminate the shell workflow.
die() {
    printf 'Error: %s\n' "$1" >&2
    exit 1
}

# trim_cr(): Remove a trailing carriage return from a Windows-format input line.
trim_cr() {
    local value="$1"
    value=${value%$'\r'}
    printf '%s' "$value"
}

# gene_name_from_path(): Derive the gene key from a FASTA basename, ignoring extension case.
gene_name_from_path() {
    local file_name
    file_name=$(basename "$1")
    case "${file_name,,}" in
        *.fasta) printf '%s\n' "${file_name:0:${#file_name}-6}" ;;
        *.fna) printf '%s\n' "${file_name:0:${#file_name}-4}" ;;
        *.fa) printf '%s\n' "${file_name:0:${#file_name}-3}" ;;
        *) return 1 ;;
    esac
}

# collect_gene_files(): Populate the named associative array with gene-to-FASTA paths; reject duplicate gene keys.
collect_gene_files() {
    local dir="$1"
    local target="$2"
    local -n gene_map="$target"
    local files=()
    local file=""
    local gene=""

    gene_map=()

    shopt -s nullglob nocaseglob
    files=( "$dir"/*.fna "$dir"/*.fasta "$dir"/*.fa )
    shopt -u nullglob nocaseglob

    if [ "${#files[@]}" -eq 0 ]; then
        return 1
    fi

    while IFS= read -r file; do
        gene=$(gene_name_from_path "$file") || die "Unsupported FASTA filename encountered: $file"

        if [[ -z "$gene" ]]; then
            die "Encountered a FASTA file without a gene name before the suffix: $file"
        fi

        if [[ "$gene" == *.* ]]; then
            die "FASTA filenames must only contain the gene name before the FASTA suffix: $file"
        fi

        if [[ -n "${gene_map[$gene]+x}" ]]; then
            die "Gene '$gene' occurs more than once in directory '$dir'. Use one FASTA file per gene."
        fi

        gene_map["$gene"]="$file"
    done < <(printf '%s\n' "${files[@]}" | awk '!seen[$0]++' | sort)

    return 0
}

# validate_multi_sequence_fasta(): Check multi-sample FASTA records for sequence data, unique names, and species-formatted headers.
validate_multi_sequence_fasta() {
    local fasta_file="$1"
    local context="$2"
    local line=""
    local line_number=0
    local header=""
    local file_failed=0
    local header_seen=0
    local sequence_seen=0
    local -a special_headers=()
    declare -A seen_headers=()

    while IFS= read -r line || [ -n "$line" ]; do
        line_number=$((line_number + 1))
        line=$(trim_cr "$line")

        [ -n "$line" ] || continue

        if [[ "$line" == ">"* ]]; then
            header="${line#>}"
            header_seen=1

            if [ -z "$header" ]; then
                printf 'Error: %s:%s contains an empty FASTA header.\n' "$fasta_file" "$line_number" >&2
                file_failed=1
                continue
            fi

            if [ -n "${seen_headers[$header]+x}" ]; then
                printf 'Error: Duplicate sequence name in %s: %s\n' "$fasta_file" "$header" >&2
                file_failed=1
            else
                seen_headers["$header"]=1
            fi

            if [[ "$header" =~ [[:space:]] ]]; then
                printf 'Error: Sequence name in %s contains whitespace and should use underscores as separators: %s\n' "$fasta_file" "$header" >&2
                file_failed=1
            fi

            if [[ ! "$header" =~ ^[A-Z][A-Za-z]+_[a-z][A-Za-z]+ ]]; then
                printf 'Error: Sequence name in %s does not start with <Genus>_<species>: %s\n' "$fasta_file" "$header" >&2
                file_failed=1
            fi

            if [[ "$header" != *_*_* ]]; then
                warn "Sequence name in $fasta_file does not include a sample-specific identifier after <Genus>_<species>: $header"
            fi

            if [[ ! "$header" =~ ^[A-Za-z0-9_]+$ ]]; then
                special_headers+=( "$header" )
            fi
        elif [ "$header_seen" -eq 0 ]; then
            printf 'Error: %s:%s contains sequence data before the first FASTA header.\n' "$fasta_file" "$line_number" >&2
            file_failed=1
        else
            sequence_seen=1
        fi
    done < "$fasta_file"

    if [ "$header_seen" -eq 0 ]; then
        printf 'Error: No FASTA headers found in %s.\n' "$fasta_file" >&2
        file_failed=1
    fi

    if [ "$sequence_seen" -eq 0 ]; then
        printf 'Error: No FASTA sequence data found in %s.\n' "$fasta_file" >&2
        file_failed=1
    fi

    if [ "${#special_headers[@]}" -gt 0 ]; then
        warn "Special characters were found in sequence names in $context file $fasta_file. This could cause issues with downstream data processing, and sequence names should only include underscores as separators."
        printf '%s\n' "${special_headers[@]}" >&2
    fi

    return "$file_failed"
}

# require_file(): Fail early if the required input file does not exist.
require_file() {
    local file="$1"
    [ -f "$file" ] || die "File not found: $file"
}

# require_csv_extension(): Reject input paths without the expected CSV filename extension.
require_csv_extension() {
    local file="$1"
    if [[ "${file##*.}" != "csv" && "${file##*.}" != "CSV" ]]; then
        die "Expected a comma-separated .csv file: $file"
    fi
}

# validate_csv_has_commas(): Reject files whose header does not look comma-separated.
validate_csv_has_commas() {
    local file="$1"
    local header=""
    header=$(awk 'NF { gsub(/\r$/, "", $0); print; exit }' "$file")
    [ -n "$header" ] || die "CSV file is empty: $file"
    [[ "$header" == *,* ]] || die "CSV file does not appear to be comma-separated: $file"
}

# validate_gene_performance_file(): Run the shared AWK CSV checker, preserving NA warnings and detailed errors.
validate_gene_performance_file() {
    local file="$1"
    [ -r "$file" ] || die "Gene performance file is missing or unreadable: $file"
    require_csv_extension "$file"
    if LC_ALL=C awk -f "$SCRIPT_DIR/gene_performance.awk" "$file"; then
        log "Gene performance calibration file format check passed."
    else
        die "Gene performance file check failed: $file"
    fi
}

# validate_thresholds_file(): Validate the eight-parameter filtering CSV before downstream filtering.
validate_thresholds_file() {
    local file="$1"
    local status=0

    [ -f "$file" ] || die "File not found: $file"
    require_csv_extension "$file"
    if LC_ALL=C awk -v enforce_ranges=1 -f "$SCRIPT_DIR/filtering_thresholds.awk" "$file"; then
        status=0
    else
        status=$?
    fi
    case $status in
        0) log "Filtering thresholds calibration file format check passed." ;;
        10) die "Filtering thresholds file must use parameter,value columns or the eight parameter names followed by one value row: $file" ;;
        11) die "Filtering thresholds file row width does not match its header: $file" ;;
        12) die "Filtering thresholds file contains an unexpected parameter: $file" ;;
        13) die "Filtering thresholds file contains a duplicated parameter: $file" ;;
        14) die "Filtering thresholds file contains a missing or non-numeric threshold value: $file" ;;
        15) die "Filtering thresholds file must contain each of the eight required parameters exactly once: $file" ;;
        16) die "Filtering thresholds file contains a value outside the allowed range: $file" ;;
        17) die "Legacy filtering thresholds file must contain exactly one value row: $file" ;;
        *) die "Unable to validate filtering thresholds file: $file" ;;
    esac
}

# validate_bins_value(): Require an integer confidence bin count between 1 and 100.
validate_bins_value() {
    local bins="$1"

    [[ "$bins" =~ ^[0-9]+$ ]] || die "Bins must be a positive integer."
    [ "$bins" -ge 1 ] || die "Bins must be at least 1."
    [ "$bins" -le 100 ] || die "Bins must be at most 100."
}

# --- Confidence workflow commands ---
# run_prepare(): Prepare references and confidence BLAST matches, then label matches using species groups or genera.
run_prepare() {
    local reference_dir=""
    local confidence_dir=""
    local species_groups_file=""
    local input_checks_failed=0
    local matched_genes=()
    local missing_reference_genes=()
    local gene=""
    local sample_count=0
    local blast_file="$PREPARATIONS_DIR/confidence_blast.tsv"
    local prepared_file="$PREPARATIONS_DIR/confidence_prepared.rds"

    if [ "$#" -eq 0 ]; then
        usage_prepare
        exit 1
    fi

    while getopts ":r:i:s:h" opt; do
        case "$opt" in
            r) reference_dir="$OPTARG" ;;
            i) confidence_dir="$OPTARG" ;;
            s) species_groups_file="$OPTARG" ;;
            h)
                usage_prepare
                exit 0
                ;;
            :)
                die "Option -$OPTARG requires an argument. Use -h for help."
                ;;
            \?)
                die "Unknown option: -$OPTARG. Use -h for help."
                ;;
        esac
    done

    [ -n "$reference_dir" ] || die "Reference directory is required. Use -r <reference directory>."
    [ -n "$confidence_dir" ] || die "Confidence dataset directory is required. Use -i <confidence dataset directory>."
    [ -d "$reference_dir" ] || die "Reference directory not found: $reference_dir"
    [ -d "$confidence_dir" ] || die "Confidence dataset directory not found: $confidence_dir"
    [ -z "$species_groups_file" ] || require_file "$species_groups_file"

    reference_dir=${reference_dir%/}
    confidence_dir=${confidence_dir%/}

    log "Checking and preparing reference dataset..."
    bash "$SCRIPTS_DIR/reference.sh" -r "$reference_dir"

    collect_gene_files "$reference_dir" REFERENCE_GENE_FILES || die "No FASTA gene files found in reference directory. Expected files ending in .FNA, .fasta or .fa."
    collect_gene_files "$confidence_dir" CONFIDENCE_GENE_FILES || die "No FASTA gene files found in confidence dataset directory. Expected files ending in .FNA, .fasta or .fa."

    log "Checking confidence dataset..."
    for gene in "${!CONFIDENCE_GENE_FILES[@]}"; do
        log "Checking $(basename "${CONFIDENCE_GENE_FILES[$gene]}")"
        if ! validate_multi_sequence_fasta "${CONFIDENCE_GENE_FILES[$gene]}" "confidence"; then
            input_checks_failed=1
        fi
    done

    if [ "$input_checks_failed" -ne 0 ]; then
        die "Confidence dataset checks failed. Please fix the FASTA files and run the script again."
    fi

    for gene in "${!CONFIDENCE_GENE_FILES[@]}"; do
        if [ -n "${REFERENCE_GENE_FILES[$gene]+x}" ]; then
            matched_genes+=( "$gene" )
        else
            missing_reference_genes+=( "$gene" )
            warn "Gene found in the confidence dataset but not in the reference dataset: $gene"
        fi
    done

    if [ "${#matched_genes[@]}" -eq 0 ]; then
        die "None of the genes in the confidence dataset were found in the reference dataset."
    fi

    if [ "${#missing_reference_genes[@]}" -eq 0 ]; then
        log "All confidence gene names were found in the reference dataset."
    fi

    sample_count=$(
        for gene in "${matched_genes[@]}"; do
            awk '/^>/ { sub(/^>/, "", $0); gsub(/\r$/, "", $0); print }' "${CONFIDENCE_GENE_FILES[$gene]}"
        done | sort -u | awk 'NF { count++ } END { print count + 0 }'
    )

    log "Number of samples being used for confidence: $sample_count"

    mkdir -p "$PREPARATIONS_DIR"
    command -v blastn >/dev/null 2>&1 || die "blastn command not found in PATH."
    command -v Rscript >/dev/null 2>&1 || die "Rscript command not found in PATH."

    log "Matching confidence genes against reference databases..."
    {
        printf 'gene\tquery\ttarget\tpident\tlength\tmismatch\tgapopen\tevalue\tbitscore\n'
        while IFS= read -r gene; do
            [ -n "$gene" ] || continue

            blastn \
                -query "${CONFIDENCE_GENE_FILES[$gene]}" \
                -db "${REFERENCE_GENE_FILES[$gene]}" \
                -task megablast \
                -outfmt "6 qseqid sseqid pident length mismatch gapopen evalue bitscore" \
                -max_target_seqs 1000000 |
                sort -t $'\t' -k1,1 -k8,8rn -k2,2R |
                awk '!seen[$1]++' |
                awk -v genename="$gene" '{print genename "\t" $0}'
        done < <(printf '%s\n' "${matched_genes[@]}" | sort)
    } > "$blast_file"

    log "BLAST confidence file written:"
    log "$blast_file"

    log "Preparing confidence data for downstream R analyses..."
    if [ -n "$species_groups_file" ]; then
        Rscript "$SCRIPTS_DIR/confidence_preparations.R" "$blast_file" "$prepared_file" "$species_groups_file" "$reference_dir"
    else
        Rscript "$SCRIPTS_DIR/confidence_preparations.R" "$blast_file" "$prepared_file"
    fi
}

# run_estimate(): Validate calibrated inputs and launch the confidence-bin comparison in R.
run_estimate() {
    local input_file="$DEFAULT_PREPARED_FILE"
    local gene_performance_file=""
    local thresholds_file=""

    if [ "$#" -eq 0 ]; then
        usage_estimate
        exit 1
    fi

    while getopts ":i:g:t:h" opt; do
        case "$opt" in
            i) input_file="$OPTARG" ;;
            g) gene_performance_file="$OPTARG" ;;
            t) thresholds_file="$OPTARG" ;;
            h)
                usage_estimate
                exit 0
                ;;
            :)
                die "Option -$OPTARG requires an argument. Use -h for help."
                ;;
            \?)
                die "Unknown option: -$OPTARG. Use -h for help."
                ;;
        esac
    done

    [ -n "$gene_performance_file" ] || die "Gene performance CSV is required. Use -g <gene performance CSV>."
    [ -n "$thresholds_file" ] || die "Filtering thresholds CSV is required. Use -t <filtering thresholds CSV>."

    log "Checking confidence input files..."
    require_file "$input_file"
    validate_gene_performance_file "$gene_performance_file"
    validate_thresholds_file "$thresholds_file"
    command -v Rscript >/dev/null 2>&1 || die "Rscript command not found in PATH."
    log "Done."
    log ""

    mkdir -p "$TESTS_DIR"

    log "Estimating confidence for different support bins..."
    Rscript "$SCRIPTS_DIR/confidence_estimate.R" "$input_file" "$gene_performance_file" "$thresholds_file" "$TESTS_DIR"
    log ""
    log "Inspect the confidence plot above, then select the number of bins with gpid confidence bins -b <number>."
}

# run_bins(): Check the chosen bin count and export confidence probabilities and the selected-bin plot.
run_bins() {
    local input_file="$DEFAULT_TOP_IDS_FILE"
    local bins=""

    if [ "$#" -eq 0 ]; then
        usage_bins
        exit 1
    fi

    while getopts ":i:b:h" opt; do
        case "$opt" in
            i) input_file="$OPTARG" ;;
            b) bins="$OPTARG" ;;
            h)
                usage_bins
                exit 0
                ;;
            :)
                die "Option -$OPTARG requires an argument. Use -h for help."
                ;;
            \?)
                die "Unknown option: -$OPTARG. Use -h for help."
                ;;
        esac
    done

    [ -n "$bins" ] || die "Number of bins is required. Use -b <number of bins>."
    validate_bins_value "$bins"
    require_file "$input_file"
    command -v Rscript >/dev/null 2>&1 || die "Rscript command not found in PATH."

    mkdir -p "$OUTPUT_DIR"
    Rscript "$SCRIPTS_DIR/confidence_bins.R" "$input_file" "$bins" "$OUTPUT_DIR"
}

# --- Dispatch the confidence subcommand ---
if [ "$#" -eq 0 ]; then
    usage
    exit 1
fi

command_name="$1"
shift

case "$command_name" in
    prepare)
        run_prepare "$@"
        ;;
    estimate)
        run_estimate "$@"
        ;;
    bins)
        run_bins "$@"
        ;;
    help|-h|--help)
        usage
        ;;
    *)
        printf 'Error: Unknown confidence command: %s\n' "$command_name" >&2
        usage >&2
        exit 1
        ;;
esac
