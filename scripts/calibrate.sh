#!/bin/bash

#################################################################
# GeneParliamentID method calibration                           #
# Benedikt Kuhnhaeuser                                          #
# Royal Botanic Gardens, Kew                                    #
# 2026                                                          #
#################################################################

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
SCRIPTS_DIR="$SCRIPT_DIR"
VERSION_FILE="${GPID_VERSION_FILE:-$PROJECT_DIR/VERSION}"
if [ -z "${GPID_VERSION:-}" ] && [ -f "$VERSION_FILE" ]; then
    IFS= read -r GPID_VERSION < "$VERSION_FILE" || GPID_VERSION=""
    GPID_VERSION=${GPID_VERSION%$'\r'}
fi
GPID_VERSION="${GPID_VERSION:-unknown}"
TEMPLATE_FILE="$PROJECT_DIR/templates/calibration_filtering_thresholds_template.csv"
OUTPUT_DIR="calibration"
PREPARATIONS_OUTPUT_DIR="$OUTPUT_DIR/preparations"
TEST_OUTPUT_DIR="$OUTPUT_DIR/tests"
ALIGNMENTS_TEST_OUTPUT_DIR="$TEST_OUTPUT_DIR/alignments"
GENES_TEST_OUTPUT_DIR="$TEST_OUTPUT_DIR/genes"
PARLIAMENT_TEST_OUTPUT_DIR="$TEST_OUTPUT_DIR/parliament"
DEFAULT_PREPARED_FILE="$PREPARATIONS_OUTPUT_DIR/calibration_prepared.rds"

BLAST_SUFFIXES=(
    ".ndb"
    ".nhr"
    ".nin"
    ".njs"
    ".nog"
    ".nos"
    ".not"
    ".nsq"
    ".ntf"
    ".nto"
)

declare -A REFERENCE_GENE_FILES=()
declare -A CALIBRATION_GENE_FILES=()

usage() {
    printf 'GPID version: %s\n\n' "$GPID_VERSION"
    cat <<'EOF'
Usage: gpid calibrate <command> [arguments]

Commands:
  prepare       Prepare BLAST and R input data for calibration
  alignments    Identify optimal alignment filtering thresholds
  genes         Estimate gene performance and identify gene threshold
  parliament    Identify minimum parliament size threshold
  combine       Combine manually selected thresholds

Examples:
  gpid calibrate prepare -r reference -i calibration_samples
  gpid calibrate alignments
  gpid calibrate genes -a calibration/manual_input_needed/calibration_alignments.tsv
  gpid calibrate parliament -a calibration/manual_input_needed/calibration_alignments.tsv -g calibration/manual_input_needed/calibration_genes.tsv
  gpid calibrate combine -a calibration/manual_input_needed/calibration_alignments.tsv -g calibration/manual_input_needed/calibration_genes.tsv -p calibration/manual_input_needed/calibration_parliament.tsv
EOF
}

usage_preparations() {
    printf 'GPID version: %s\n\n' "$GPID_VERSION"
    cat <<'EOF'
Usage: gpid calibrate prepare -r <reference directory> -i <calibration dataset directory>

Required:
  -r  Reference dataset directory containing one FASTA file per gene and BLAST databases
      This is usually prepared with gpid reference.
  -i  Calibration dataset directory containing one FASTA file per gene

Default outputs created for later calibration steps:
  calibration/preparations/calibration_blast.tsv
  calibration/preparations/calibration_prepared.rds
EOF
}

usage_alignments() {
    printf 'GPID version: %s\n\n' "$GPID_VERSION"
    cat <<'EOF'
Usage: gpid calibrate alignments [-i <prepared calibration RDS>]

Optional:
  -i  Intermediate RDS produced by gpid calibrate prepare
      Default: calibration/preparations/calibration_prepared.rds

The default input file is created by gpid calibrate prepare and stored as:
  calibration/preparations/calibration_prepared.rds

Default output created for the next calibration step:
  calibration/manual_input_needed/calibration_alignments.tsv

Calibration test CSV/PDF outputs are written to:
  calibration/tests/alignments/
EOF
}

usage_genes() {
    printf 'GPID version: %s\n\n' "$GPID_VERSION"
    cat <<'EOF'
Usage: gpid calibrate genes [-i <prepared calibration RDS>] -a <alignment thresholds TSV>

Required:
  -a  Alignment thresholds TSV produced and manually edited after gpid calibrate alignments
      Default output path from gpid calibrate alignments:
        calibration/manual_input_needed/calibration_alignments.tsv

Optional:
  -i  Intermediate RDS produced by gpid calibrate prepare
      Default: calibration/preparations/calibration_prepared.rds

The default input file is created by gpid calibrate prepare and stored as:
  calibration/preparations/calibration_prepared.rds

Default output created for the next calibration step:
  calibration/manual_input_needed/calibration_genes.tsv

Additional output created by this command:
  calibration/calibration_gene_performance.csv

Calibration test CSV/PDF outputs are written to:
  calibration/tests/genes/
EOF
}

usage_parliament() {
    printf 'GPID version: %s\n\n' "$GPID_VERSION"
    cat <<'EOF'
Usage: gpid calibrate parliament [-i <prepared calibration RDS>] -a <alignment thresholds TSV> -g <gene threshold TSV>

Required:
  -a  Alignment thresholds TSV produced and manually edited after gpid calibrate alignments
      Default output path from gpid calibrate alignments:
        calibration/manual_input_needed/calibration_alignments.tsv
  -g  Gene threshold TSV produced and manually edited after gpid calibrate genes
      Default output path from gpid calibrate genes:
        calibration/manual_input_needed/calibration_genes.tsv

Optional:
  -i  Intermediate RDS produced by gpid calibrate prepare
      Default: calibration/preparations/calibration_prepared.rds

The default input file is created by gpid calibrate prepare and stored as:
  calibration/preparations/calibration_prepared.rds

Default output created for the combine step:
  calibration/manual_input_needed/calibration_parliament.tsv

Calibration test CSV/PDF outputs are written to:
  calibration/tests/parliament/
EOF
}

usage_combine() {
    printf 'GPID version: %s\n\n' "$GPID_VERSION"
    cat <<'EOF'
Usage: gpid calibrate combine -a <alignment thresholds TSV> -g <gene threshold TSV> -p <parliament threshold TSV>

Required:
  -a  Alignment thresholds TSV produced and manually edited after gpid calibrate alignments
      Default output path from gpid calibrate alignments:
        calibration/manual_input_needed/calibration_alignments.tsv
  -g  Gene threshold TSV produced and manually edited after gpid calibrate genes
      Default output path from gpid calibrate genes:
        calibration/manual_input_needed/calibration_genes.tsv
  -p  Parliament threshold TSV produced and manually edited after gpid calibrate parliament
      Default output path from gpid calibrate parliament:
        calibration/manual_input_needed/calibration_parliament.tsv

Default output created by this command:
  calibration/calibration_filtering_thresholds.csv
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

blast_db_complete() {
    local reference_fasta="$1"
    local suffix=""

    for suffix in "${BLAST_SUFFIXES[@]}"; do
        if [ ! -f "${reference_fasta}${suffix}" ]; then
            return 1
        fi
    done

    return 0
}

csv_header() {
    local file="$1"
    awk 'NF { gsub(/\r$/, "", $0); print; exit }' "$file"
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

validate_single_threshold_file() {
    local file="$1"
    local expected_parameter="$2"
    local value_label="$3"
    local min_value="$4"
    local max_value="$5"
    local status=0

    require_file "$file"

    if awk -F'\t' -v expected_parameter="$expected_parameter" -v min="$min_value" -v max="$max_value" '
        function trim(value) {
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
            return value
        }
        BEGIN { row = 0 }
        NF {
            gsub(/\r$/, "", $0)
            row++

            if (row == 1) {
                if (NF != 2 || trim($1) != "parameter" || trim($2) != "value") {
                    exit 10
                }
                next
            }

            if (row == 2) {
                if (NF != 2) {
                    exit 11
                }

                parameter = trim($1)
                value = trim($2)

                if (parameter != expected_parameter) {
                    exit 18
                }

                if (value == "NA" || value == "") {
                    exit 12
                }

                if (value !~ /^[0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?$/) {
                    exit 13
                }

                if ((value + 0) < min || (value + 0) > max) {
                    exit 14
                }
                next
            }

            exit 15
        }
        END {
            if (row == 0) {
                exit 16
            }
            if (row == 1) {
                exit 17
            }
        }
    ' "$file"; then
        status=0
    else
        status=$?
    fi

    case $status in
        0) ;;
        10) die "$value_label file must use a tab-separated header: parameter<TAB>value: $file" ;;
        11) die "$value_label file must contain exactly one threshold value in the second row: $file" ;;
        12) die "All values need to be specified before this file can be used: $file" ;;
        13) die "$value_label file contains a non-numeric threshold value: $file" ;;
        14) die "$value_label must be between $min_value and $max_value: $file" ;;
        15) die "$value_label file must contain only a header row and one row of thresholds: $file" ;;
        16) die "$value_label file is empty: $file" ;;
        17) die "$value_label file is missing the thresholds row: $file" ;;
        18) die "$value_label file must contain parameter '$expected_parameter': $file" ;;
        *) die "Unable to validate $value_label file: $file" ;;
    esac
}

threshold_value() {
    local file="$1"
    awk -F'\t' '
        function trim(value) {
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
            return value
        }
        NF {
            gsub(/\r$/, "", $0)
            if (trim($1) == "parameter") {
                next
            }
            print trim($2)
            exit
        }
    ' "$file"
}

validate_alignment_thresholds_for_combine() {
    local file="$1"
    local expected_parameters="min_similarity min_length max_gapopens max_mismatches max_evalue min_bitscore"
    local status=0

    require_file "$file"

    if awk -F'\t' -v expected_parameters="$expected_parameters" '
        function trim(value) {
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
            return value
        }
        BEGIN {
            row = 0
            split(expected_parameters, expected, " ")
            split("0 0 0 0 0 0", mins, " ")
            split("100 99999 99999 99999 100 99999", maxs, " ")
        }
        NF {
            gsub(/\r$/, "", $0)
            row++

            if (row == 1) {
                if (NF != 2 || trim($1) != "parameter" || trim($2) != "value") {
                    exit 10
                }
                next
            }

            if (NF != 2) {
                exit 11
            }

            parameter = trim($1)
            value = trim($2)
            found = 0

            for (i = 1; i <= 6; i++) {
                if (parameter == expected[i]) {
                    found = i
                }
            }

            if (found == 0) {
                exit 18
            }

            if (seen[parameter] == 1) {
                exit 19
            }
            seen[parameter] = 1

            if (value == "NA" || value == "") {
                exit 12
            }
            if (value !~ /^[0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?$/) {
                exit 13
            }
            if ((value + 0) < mins[found] || (value + 0) > maxs[found]) {
                exit 14
            }
        }
        END {
            if (row == 0) {
                exit 16
            }
            if (row == 1) {
                exit 17
            }
            for (i = 1; i <= 6; i++) {
                if (seen[expected[i]] != 1) {
                    exit 20
                }
            }
        }
    ' "$file"; then
        status=0
    else
        status=$?
    fi

    case $status in
        0) ;;
        10) die "Alignment thresholds file must use a tab-separated header: parameter<TAB>value: $file" ;;
        11) die "Alignment thresholds file must contain two tab-separated columns per row: $file" ;;
        12) die "All values need to be specified before this file can be used: $file" ;;
        13) die "Alignment thresholds file contains a non-numeric threshold value: $file" ;;
        14) die "Alignment thresholds file contains a value outside the allowed range: $file" ;;
        16) die "Alignment thresholds file is empty: $file" ;;
        17) die "Alignment thresholds file is missing threshold rows: $file" ;;
        18) die "Alignment thresholds file contains an unexpected parameter name: $file" ;;
        19) die "Alignment thresholds file contains a duplicated parameter name: $file" ;;
        20) die "Alignment thresholds file is missing one or more required parameters: $file" ;;
        *) die "Unable to validate alignment thresholds file: $file" ;;
    esac
}

print_manual_threshold_instructions() {
    cat <<'EOF'
MANUAL INPUT NEEDED: Specify the optimal alignment filtering thresholds in the file calibration/manual_input_needed/calibration_alignments.tsv based on the calibration test files in directory calibration/tests/alignments.
EOF
}

run_preparations() {
    local reference_dir=""
    local calibration_dir=""
    local validation_failed=0
    local matched_genes=()
    local missing_reference_genes=()
    local missing_databases=()
    local gene=""
    local sample_count=0
    local blast_file="$PREPARATIONS_OUTPUT_DIR/calibration_blast.tsv"
    local prepared_file="$PREPARATIONS_OUTPUT_DIR/calibration_prepared.rds"

    if [ "$#" -eq 0 ]; then
        usage_preparations
        exit 1
    fi

    while getopts ":r:i:h" opt; do
        case "$opt" in
            r) reference_dir="$OPTARG" ;;
            i) calibration_dir="$OPTARG" ;;
            h)
                usage_preparations
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
    [ -n "$calibration_dir" ] || die "Calibration dataset directory is required. Use -i <calibration dataset directory>."
    [ -d "$reference_dir" ] || die "Reference directory not found: $reference_dir"
    [ -d "$calibration_dir" ] || die "Calibration dataset directory not found: $calibration_dir"

    reference_dir=${reference_dir%/}
    calibration_dir=${calibration_dir%/}

    log "Checking and preparing reference dataset..."
    bash "$SCRIPTS_DIR/reference.sh" -r "$reference_dir"

    collect_gene_files "$reference_dir" REFERENCE_GENE_FILES || die "No FASTA gene files found in reference directory. Expected files ending in .FNA, .fasta or .fa."
    collect_gene_files "$calibration_dir" CALIBRATION_GENE_FILES || die "No FASTA gene files found in calibration dataset directory. Expected files ending in .FNA, .fasta or .fa."

    log "Checking calibration dataset..."
    for gene in "${!CALIBRATION_GENE_FILES[@]}"; do
        log "Checking $(basename "${CALIBRATION_GENE_FILES[$gene]}")"
        if ! validate_multi_sequence_fasta "${CALIBRATION_GENE_FILES[$gene]}" "calibration"; then
            validation_failed=1
        fi
    done

    if [ "$validation_failed" -ne 0 ]; then
        die "Calibration dataset validation failed. Please fix the FASTA files and run the script again."
    fi

    for gene in "${!CALIBRATION_GENE_FILES[@]}"; do
        if [ -n "${REFERENCE_GENE_FILES[$gene]+x}" ]; then
            matched_genes+=( "$gene" )
        else
            missing_reference_genes+=( "$gene" )
            warn "Gene found in the calibration dataset but not in the reference dataset: $gene"
        fi
    done

    if [ "${#matched_genes[@]}" -eq 0 ]; then
        die "None of the genes in the calibration dataset were found in the reference dataset."
    fi

    if [ "${#missing_reference_genes[@]}" -eq 0 ]; then
        log "All calibration gene names were found in the reference dataset."
    fi

    for gene in "${matched_genes[@]}"; do
        if ! blast_db_complete "${REFERENCE_GENE_FILES[$gene]}"; then
            missing_databases+=( "$gene" )
        fi
    done

    if [ "${#missing_databases[@]}" -gt 0 ]; then
        while IFS= read -r gene; do
            [ -n "$gene" ] || continue
            warn "Missing BLAST database files for reference gene: $gene"
        done < <(printf '%s\n' "${missing_databases[@]}" | sort)

        die "Reference BLAST databases are missing. Please run gpid reference for the reference directory first."
    fi

    sample_count=$(
        for gene in "${matched_genes[@]}"; do
            awk '/^>/ { sub(/^>/, "", $0); gsub(/\r$/, "", $0); print }' "${CALIBRATION_GENE_FILES[$gene]}"
        done | sort -u | awk 'NF { count++ } END { print count + 0 }'
    )

    log "Number of samples being used for calibration: $sample_count"

    mkdir -p "$PREPARATIONS_OUTPUT_DIR"
    command -v blastn >/dev/null 2>&1 || die "blastn command not found in PATH."
    command -v Rscript >/dev/null 2>&1 || die "Rscript command not found in PATH."

    log "Matching calibration genes against reference databases..."
    {
        printf 'gene\tquery\ttarget\tpident\tlength\tmismatch\tgapopen\tevalue\tbitscore\n'
        while IFS= read -r gene; do
            [ -n "$gene" ] || continue

            blastn \
                -query "${CALIBRATION_GENE_FILES[$gene]}" \
                -db "${REFERENCE_GENE_FILES[$gene]}" \
                -task megablast \
                -outfmt "6 qseqid sseqid pident length mismatch gapopen evalue bitscore" \
                -max_target_seqs 1000000 |
                sort -t $'\t' -k1,1 -k8,8rn -k2,2R |
                awk '!seen[$1]++' |
                awk -v genename="$gene" '{print genename "\t" $0}'
        done < <(printf '%s\n' "${matched_genes[@]}" | sort)
    } > "$blast_file"

    log "BLAST calibration file written:"
    log "$blast_file"

    log "Preparing calibration data for downstream R analyses..."
    Rscript "$SCRIPTS_DIR/calibration_preparations.R" "$blast_file" "$prepared_file"

    log "Prepared calibration data written:"
    log "$prepared_file"
}

run_alignments() {
    local input_file="$DEFAULT_PREPARED_FILE"

    while getopts ":i:h" opt; do
        case "$opt" in
            i) input_file="$OPTARG" ;;
            h)
                usage_alignments
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

    require_file "$input_file"
    command -v Rscript >/dev/null 2>&1 || die "Rscript command not found in PATH."
    mkdir -p "$OUTPUT_DIR" "$ALIGNMENTS_TEST_OUTPUT_DIR"

    Rscript "$SCRIPTS_DIR/calibration_alignments.R" "$input_file" "$OUTPUT_DIR" "$ALIGNMENTS_TEST_OUTPUT_DIR" "$TEMPLATE_FILE"
    print_manual_threshold_instructions
}

run_genes() {
    local input_file="$DEFAULT_PREPARED_FILE"
    local alignments_file=""

    if [ "$#" -eq 0 ]; then
        usage_genes
        exit 1
    fi

    while getopts ":i:a:h" opt; do
        case "$opt" in
            i) input_file="$OPTARG" ;;
            a) alignments_file="$OPTARG" ;;
            h)
                usage_genes
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

    [ -n "$alignments_file" ] || die "Alignment thresholds TSV is required. Use -a <alignment thresholds TSV>."
    require_file "$input_file"
    require_file "$alignments_file"
    command -v Rscript >/dev/null 2>&1 || die "Rscript command not found in PATH."
    mkdir -p "$OUTPUT_DIR" "$GENES_TEST_OUTPUT_DIR"

    Rscript "$SCRIPTS_DIR/calibration_genes.R" "$input_file" "$alignments_file" "$OUTPUT_DIR" "$GENES_TEST_OUTPUT_DIR" "$TEMPLATE_FILE"
    log "MANUAL INPUT NEEDED: Specify the optimal gene performance threshold in the file calibration/manual_input_needed/calibration_genes.tsv based on the calibration test files in directory calibration/tests/genes."
}

run_parliament() {
    local input_file="$DEFAULT_PREPARED_FILE"
    local alignments_file=""
    local genes_file=""

    if [ "$#" -eq 0 ]; then
        usage_parliament
        exit 1
    fi

    while getopts ":i:a:g:h" opt; do
        case "$opt" in
            i) input_file="$OPTARG" ;;
            a) alignments_file="$OPTARG" ;;
            g) genes_file="$OPTARG" ;;
            h)
                usage_parliament
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

    [ -n "$alignments_file" ] || die "Alignment thresholds TSV is required. Use -a <alignment thresholds TSV>."
    [ -n "$genes_file" ] || die "Gene threshold TSV is required. Use -g <gene threshold TSV>."
    require_file "$input_file"
    require_file "$alignments_file"
    require_file "$genes_file"
    command -v Rscript >/dev/null 2>&1 || die "Rscript command not found in PATH."
    mkdir -p "$OUTPUT_DIR" "$PARLIAMENT_TEST_OUTPUT_DIR"

    Rscript "$SCRIPTS_DIR/calibration_parliament.R" "$input_file" "$alignments_file" "$genes_file" "$OUTPUT_DIR" "$PARLIAMENT_TEST_OUTPUT_DIR" "$TEMPLATE_FILE"
    log "MANUAL INPUT NEEDED: Specify the optimal minimum parliament size threshold in the file calibration/manual_input_needed/calibration_parliament.tsv based on the calibration test files in directory calibration/tests/parliament."
}

run_combine() {
    local alignments_file=""
    local genes_file=""
    local parliament_file=""
    local output_file="$OUTPUT_DIR/calibration_filtering_thresholds.csv"
    local alignment_values=""
    local gene_value=""
    local parliament_value=""

    if [ "$#" -eq 0 ]; then
        usage_combine
        exit 1
    fi

    while getopts ":a:g:p:h" opt; do
        case "$opt" in
            a) alignments_file="$OPTARG" ;;
            g) genes_file="$OPTARG" ;;
            p) parliament_file="$OPTARG" ;;
            h)
                usage_combine
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

    [ -n "$alignments_file" ] || die "Alignment thresholds TSV is required. Use -a <alignment thresholds TSV>."
    [ -n "$genes_file" ] || die "Gene threshold TSV is required. Use -g <gene threshold TSV>."
    [ -n "$parliament_file" ] || die "Parliament threshold TSV is required. Use -p <parliament threshold TSV>."

    log "Step 1/5: Checking alignment filtering thresholds."
    validate_alignment_thresholds_for_combine "$alignments_file"

    log "Step 2/5: Checking gene performance threshold."
    validate_single_threshold_file "$genes_file" "min_gene_performance" "Gene performance threshold" 0 100

    log "Step 3/5: Checking parliament size threshold."
    validate_single_threshold_file "$parliament_file" "min_parliament_size" "Parliament size threshold" 0 99999

    log "Step 4/5: Reading specified calibration thresholds."
    alignment_values=$(awk -F'\t' '
        function trim(value) {
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
            return value
        }
        NF {
            gsub(/\r$/, "", $0)
            parameter = trim($1)
            if (parameter == "parameter") {
                next
            }
            values[parameter] = trim($2)
        }
        END {
            print values["min_similarity"] "," values["min_length"] "," values["max_gapopens"] "," values["max_mismatches"] "," values["max_evalue"] "," values["min_bitscore"]
        }
    ' "$alignments_file")
    gene_value=$(threshold_value "$genes_file")
    parliament_value=$(threshold_value "$parliament_file")

    log "Step 5/5: Writing combined filtering thresholds CSV."
    mkdir -p "$OUTPUT_DIR"
    {
        csv_header "$TEMPLATE_FILE"
        printf '%s,%s,%s\n' "$alignment_values" "$gene_value" "$parliament_value"
    } > "$output_file"

    IFS=',' read -r min_similarity min_length max_gapopens max_mismatches max_evalue min_bitscore <<< "$alignment_values"

    log "Specified calibration thresholds:"
    log "Minimum alignment similarity: $min_similarity"
    log "Minimum alignment length: $min_length"
    log "Maximum number of gap openings: $max_gapopens"
    log "Maximum number of alignment mismatches: $max_mismatches"
    log "Maximum E-value: $max_evalue"
    log "Minimum Bit-score: $min_bitscore"
    log "Minimum gene performance: $gene_value"
    log "Minimum gene parliament size: $parliament_value"
    log "Combined filtering thresholds written:"
    log "$output_file"
}

if [ "$#" -eq 0 ]; then
    usage
    exit 1
fi

command_name="$1"
shift

case "$command_name" in
    prepare)
        run_preparations "$@"
        ;;
    alignments)
        run_alignments "$@"
        ;;
    genes)
        run_genes "$@"
        ;;
    parliament)
        run_parliament "$@"
        ;;
    combine)
        run_combine "$@"
        ;;
    help|-h|--help)
        usage
        ;;
    *)
        printf 'Error: Unknown calibration command: %s\n' "$command_name" >&2
        usage >&2
        exit 1
        ;;
esac
