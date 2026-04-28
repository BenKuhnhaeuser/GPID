#!/bin/bash

#################################################################
# GeneParliamentID method validation                            #
# Benedikt Kuhnhaeuser                                          #
# Royal Botanic Gardens, Kew                                    #
# 2026                                                          #
#################################################################

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SCRIPTS_DIR="$SCRIPT_DIR"
PROJECT_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
VERSION_FILE="${GPID_VERSION_FILE:-$PROJECT_DIR/VERSION}"
if [ -z "${GPID_VERSION:-}" ] && [ -f "$VERSION_FILE" ]; then
    IFS= read -r GPID_VERSION < "$VERSION_FILE" || GPID_VERSION=""
    GPID_VERSION=${GPID_VERSION%$'\r'}
fi
GPID_VERSION="${GPID_VERSION:-unknown}"
OUTPUT_DIR="validation"
PREPARATIONS_DIR="$OUTPUT_DIR/preparations"
TESTS_DIR="$OUTPUT_DIR/tests"
DEFAULT_PREPARED_FILE="$PREPARATIONS_DIR/validation_prepared.rds"
DEFAULT_TOP_IDS_FILE="$TESTS_DIR/validate_top_ids.rds"

declare -A REFERENCE_GENE_FILES=()
declare -A VALIDATION_GENE_FILES=()

usage() {
    printf 'GPID version: %s\n\n' "$GPID_VERSION"
    cat <<'EOF'
Usage: gpid validate <command> [arguments]

Commands:
  prepare       Prepare BLAST and R input data for validation
  confidence    Estimate validation confidence for different support bins
  bins          Save validation confidence support probabilities for selected bins

Examples:
  gpid validate prepare -r reference -i validation_samples
  gpid validate confidence -i validation/preparations/validation_prepared.rds -g calibration/calibration_gene_performance.csv -t calibration/calibration_filtering_thresholds.csv
  gpid validate bins -b 5

Calibration inputs normally produced before validation:
  -g  calibration/calibration_gene_performance.csv
      Gene performance file with columns gene,performance
  -t  calibration/calibration_filtering_thresholds.csv
      Filtering thresholds file produced by gpid calibrate combine
EOF
}

usage_prepare() {
    printf 'GPID version: %s\n\n' "$GPID_VERSION"
    cat <<'EOF'
Usage: gpid validate prepare -r <reference directory> -i <validation dataset directory> [-s <species groups file>]

Required:
  -r  Reference dataset directory containing one FASTA file per gene and BLAST databases
  -i  Validation dataset directory containing one FASTA file per gene

Optional:
  -s  Species groups CSV with header genus_species,species_group
      If omitted, species groups are derived from genus names.

Default outputs:
  validation/preparations/validation_blast.tsv
  validation/preparations/validation_prepared.rds
EOF
}

usage_confidence() {
    printf 'GPID version: %s\n\n' "$GPID_VERSION"
    cat <<'EOF'
Usage: gpid validate confidence [-i <prepared validation RDS>] -g <gene performance CSV> -t <filtering thresholds CSV>

Required:
  -g  Gene performance CSV produced by calibration
      Normally saved as: calibration/calibration_gene_performance.csv
  -t  Filtering thresholds CSV produced by gpid calibrate combine
      Normally saved as: calibration/calibration_filtering_thresholds.csv

Optional:
  -i  Intermediate RDS produced by gpid validate prepare
      Default: validation/preparations/validation_prepared.rds

Default outputs:
  validation/tests/validate_confidence.pdf
  validation/tests/validate_top_ids.rds
EOF
}

usage_bins() {
    printf 'GPID version: %s\n\n' "$GPID_VERSION"
    cat <<'EOF'
Usage: gpid validate bins -b <number of bins> [-i <validate top IDs RDS>]

Required:
  -b  Number of confidence support bins to use

Optional:
  -i  Top IDs RDS produced by gpid validate confidence
      Default: validation/tests/validate_top_ids.rds

Default output:
  validation/validation_confidence_support.csv
  validation/validation_confidence_support.pdf
EOF
}

log() {
    printf '%s\n' "$1"
}

warn() {
    printf 'Warning: %s\n' "$1" >&2
}

die() {
    printf 'Error: %s\n' "$1" >&2
    exit 1
}

trim_cr() {
    local value="$1"
    value=${value%$'\r'}
    printf '%s' "$value"
}

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

require_file() {
    local file="$1"
    [ -f "$file" ] || die "File not found: $file"
}

require_csv_extension() {
    local file="$1"
    if [[ "${file##*.}" != "csv" && "${file##*.}" != "CSV" ]]; then
        die "Expected a comma-separated .csv file: $file"
    fi
}

validate_csv_has_commas() {
    local file="$1"
    local header=""
    header=$(awk 'NF { gsub(/\r$/, "", $0); print; exit }' "$file")
    [ -n "$header" ] || die "CSV file is empty: $file"
    [[ "$header" == *,* ]] || die "CSV file does not appear to be comma-separated: $file"
}

validate_gene_performance_file() {
    local file="$1"
    local status=0

    require_file "$file"
    require_csv_extension "$file"
    validate_csv_has_commas "$file"

    if awk -F',' '
        BEGIN { found = 0 }
        NF {
            gsub(/\r$/, "", $0)
            if (found == 0) {
                if (NF != 2 || $1 != "gene" || $2 != "performance") {
                    exit 10
                }
                found = 1
                next
            }

            if (NF != 2 || $1 == "" || $2 == "" || $2 !~ /^[0-9]+([.][0-9]+)?$/) {
                exit 11
            }

            if (($2 + 0) < 0 || ($2 + 0) > 100) {
                exit 12
            }
        }
        END {
            if (found == 0) {
                exit 13
            }
            if (found == 1 && NR == 1) {
                exit 14
            }
        }
    ' "$file"; then
        status=0
    else
        status=$?
    fi

    case $status in
        0) log "Gene performance calibration file format check passed." ;;
        10) die "Gene performance file must be a comma-separated CSV with header 'gene,performance': $file" ;;
        11) die "Gene performance file must contain two columns: gene name and numeric percentage performance: $file" ;;
        12) die "Gene performance values must be between 0 and 100: $file" ;;
        13) die "Gene performance file is empty: $file" ;;
        14) die "Gene performance file must contain at least one gene performance row: $file" ;;
        *) die "Unable to validate gene performance file: $file" ;;
    esac
}

validate_thresholds_file() {
    local file="$1"
    local expected_header="min_similarity,min_length,max_gapopens,max_mismatches,max_evalue,min_bitscore,min_gene_performance,min_parliament_size"
    local status=0

    require_file "$file"
    require_csv_extension "$file"
    validate_csv_has_commas "$file"

    if awk -F',' -v expected_header="$expected_header" '
        BEGIN {
            row = 0
            split("0 0 0 0 0 0 0 0", mins, " ")
            split("100 99999 99999 99999 100 99999 100 99999", maxs, " ")
        }
        NF {
            gsub(/\r$/, "", $0)
            row++

            if (row == 1) {
                if ($0 != expected_header) {
                    exit 10
                }
                next
            }

            if (row == 2) {
                if (NF != 8) {
                    exit 11
                }

                for (i = 1; i <= NF; i++) {
                    if ($i == "" || $i == "NA" || $i !~ /^[0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?$/) {
                        exit 12
                    }
                    if (($i + 0) < mins[i] || ($i + 0) > maxs[i]) {
                        exit 13
                    }
                }
                next
            }

            exit 14
        }
        END {
            if (row == 0) {
                exit 15
            }
            if (row == 1) {
                exit 16
            }
        }
    ' "$file"; then
        status=0
    else
        status=$?
    fi

    case $status in
        0) log "Filtering thresholds calibration file format check passed." ;;
        10) die "Filtering thresholds file must use the template header '$expected_header': $file" ;;
        11) die "Filtering thresholds file must contain exactly eight threshold values in the second row: $file" ;;
        12) die "All filtering threshold values need to be specified as numeric values before this file can be used: $file" ;;
        13) die "Filtering thresholds file contains a value outside the allowed range: $file" ;;
        14) die "Filtering thresholds file must contain only a header row and one row of thresholds: $file" ;;
        15) die "Filtering thresholds file is empty: $file" ;;
        16) die "Filtering thresholds file is missing the thresholds row: $file" ;;
        *) die "Unable to validate filtering thresholds file: $file" ;;
    esac
}

validate_bins_value() {
    local bins="$1"

    [[ "$bins" =~ ^[0-9]+$ ]] || die "Bins must be a positive integer."
    [ "$bins" -ge 1 ] || die "Bins must be at least 1."
    [ "$bins" -le 100 ] || die "Bins must be at most 100."
}

run_prepare() {
    local reference_dir=""
    local validation_dir=""
    local species_groups_file=""
    local validation_failed=0
    local matched_genes=()
    local missing_reference_genes=()
    local gene=""
    local sample_count=0
    local blast_file="$PREPARATIONS_DIR/validation_blast.tsv"
    local prepared_file="$PREPARATIONS_DIR/validation_prepared.rds"

    if [ "$#" -eq 0 ]; then
        usage_prepare
        exit 1
    fi

    while getopts ":r:i:s:h" opt; do
        case "$opt" in
            r) reference_dir="$OPTARG" ;;
            i) validation_dir="$OPTARG" ;;
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
    [ -n "$validation_dir" ] || die "Validation dataset directory is required. Use -i <validation dataset directory>."
    [ -d "$reference_dir" ] || die "Reference directory not found: $reference_dir"
    [ -d "$validation_dir" ] || die "Validation dataset directory not found: $validation_dir"
    [ -z "$species_groups_file" ] || require_file "$species_groups_file"

    reference_dir=${reference_dir%/}
    validation_dir=${validation_dir%/}

    log "Checking and preparing reference dataset..."
    bash "$SCRIPTS_DIR/reference.sh" -r "$reference_dir"

    collect_gene_files "$reference_dir" REFERENCE_GENE_FILES || die "No FASTA gene files found in reference directory. Expected files ending in .FNA, .fasta or .fa."
    collect_gene_files "$validation_dir" VALIDATION_GENE_FILES || die "No FASTA gene files found in validation dataset directory. Expected files ending in .FNA, .fasta or .fa."

    log "Checking validation dataset..."
    for gene in "${!VALIDATION_GENE_FILES[@]}"; do
        log "Checking $(basename "${VALIDATION_GENE_FILES[$gene]}")"
        if ! validate_multi_sequence_fasta "${VALIDATION_GENE_FILES[$gene]}" "validation"; then
            validation_failed=1
        fi
    done

    if [ "$validation_failed" -ne 0 ]; then
        die "Validation dataset validation failed. Please fix the FASTA files and run the script again."
    fi

    for gene in "${!VALIDATION_GENE_FILES[@]}"; do
        if [ -n "${REFERENCE_GENE_FILES[$gene]+x}" ]; then
            matched_genes+=( "$gene" )
        else
            missing_reference_genes+=( "$gene" )
            warn "Gene found in the validation dataset but not in the reference dataset: $gene"
        fi
    done

    if [ "${#matched_genes[@]}" -eq 0 ]; then
        die "None of the genes in the validation dataset were found in the reference dataset."
    fi

    if [ "${#missing_reference_genes[@]}" -eq 0 ]; then
        log "All validation gene names were found in the reference dataset."
    fi

    sample_count=$(
        for gene in "${matched_genes[@]}"; do
            awk '/^>/ { sub(/^>/, "", $0); gsub(/\r$/, "", $0); print }' "${VALIDATION_GENE_FILES[$gene]}"
        done | sort -u | awk 'NF { count++ } END { print count + 0 }'
    )

    log "Number of samples being used for validation: $sample_count"

    mkdir -p "$PREPARATIONS_DIR"
    command -v blastn >/dev/null 2>&1 || die "blastn command not found in PATH."
    command -v Rscript >/dev/null 2>&1 || die "Rscript command not found in PATH."

    log "Matching validation genes against reference databases..."
    {
        printf 'gene\tquery\ttarget\tpident\tlength\tmismatch\tgapopen\tevalue\tbitscore\n'
        while IFS= read -r gene; do
            [ -n "$gene" ] || continue

            blastn \
                -query "${VALIDATION_GENE_FILES[$gene]}" \
                -db "${REFERENCE_GENE_FILES[$gene]}" \
                -task megablast \
                -outfmt "6 qseqid sseqid pident length mismatch gapopen evalue bitscore" \
                -max_target_seqs 1000000 |
                sort -t $'\t' -k1,1 -k8,8rn -k2,2R |
                awk '!seen[$1]++' |
                awk -v genename="$gene" '{print genename "\t" $0}'
        done < <(printf '%s\n' "${matched_genes[@]}" | sort)
    } > "$blast_file"

    log "BLAST validation file written:"
    log "$blast_file"

    log "Preparing validation data for downstream R analyses..."
    if [ -n "$species_groups_file" ]; then
        Rscript "$SCRIPTS_DIR/validation_preparations.R" "$blast_file" "$prepared_file" "$species_groups_file"
    else
        Rscript "$SCRIPTS_DIR/validation_preparations.R" "$blast_file" "$prepared_file"
    fi

    log "Prepared validation data written:"
    log "$prepared_file"
}

run_confidence() {
    local input_file="$DEFAULT_PREPARED_FILE"
    local gene_performance_file=""
    local thresholds_file=""

    if [ "$#" -eq 0 ]; then
        usage_confidence
        exit 1
    fi

    while getopts ":i:g:t:h" opt; do
        case "$opt" in
            i) input_file="$OPTARG" ;;
            g) gene_performance_file="$OPTARG" ;;
            t) thresholds_file="$OPTARG" ;;
            h)
                usage_confidence
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

    log "Checking validation confidence input files..."
    require_file "$input_file"
    validate_gene_performance_file "$gene_performance_file"
    validate_thresholds_file "$thresholds_file"
    command -v Rscript >/dev/null 2>&1 || die "Rscript command not found in PATH."
    log "Done."
    log ""

    mkdir -p "$TESTS_DIR"

    log "Estimating validation confidence for different support bins..."
    Rscript "$SCRIPTS_DIR/validation_confidence.R" "$input_file" "$gene_performance_file" "$thresholds_file" "$TESTS_DIR"
    log "Done."
    log "Validation confidence files written:"
    log "$TESTS_DIR/validate_confidence.pdf"
    log "$TESTS_DIR/validate_top_ids.rds"
    log ""
    log "Inspect the file validate_confidence.pdf to decide on the optimal number of bins. Specify the number of bins using gpid validate bins."
}

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
    Rscript "$SCRIPTS_DIR/validation_bins.R" "$input_file" "$bins" "$OUTPUT_DIR"
}

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
    confidence)
        run_confidence "$@"
        ;;
    bins)
        run_bins "$@"
        ;;
    help|-h|--help)
        usage
        ;;
    *)
        printf 'Error: Unknown validation command: %s\n' "$command_name" >&2
        usage >&2
        exit 1
        ;;
esac
