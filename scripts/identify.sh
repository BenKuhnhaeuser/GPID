#!/bin/bash

#################################################################
# GeneParliamentID pipeline                                     #
# Benedikt Kuhnhaeuser                                          #
# Royal Botanic Gardens, Kew                                    #
# 2026                                                          #
#################################################################

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
VERSION_FILE="${GPID_VERSION_FILE:-$PROJECT_DIR/VERSION}"
if [ -z "${GPID_VERSION:-}" ] && [ -f "$VERSION_FILE" ]; then
    IFS= read -r GPID_VERSION < "$VERSION_FILE" || GPID_VERSION=""
    GPID_VERSION=${GPID_VERSION%$'\r'}
fi
GPID_VERSION="${GPID_VERSION:-unknown}"

DEFAULT_GENE_PERFORMANCE_FILE="calibration/calibration_gene_performance.csv"
DEFAULT_THRESHOLDS_FILE="calibration/calibration_filtering_thresholds.csv"
DEFAULT_CONFIDENCE_SUPPORT_FILE="validation/validation_confidence_support.csv"

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

sample_dir=""
reference_dir=""
gene_performance_file=""
thresholds_file=""
confidence_support_file=""
species_groups_file=""
output_dir=""
output_formats_raw="csv,pdf"
bypass_calibration=0
overwrite_output_dir=0
remove_intermediates=0
temp_dir=""
gene_performance_file_specified=0
thresholds_file_specified=0
confidence_support_file_specified=0

declare -A SAMPLE_GENE_FILES=()
declare -A REFERENCE_GENE_FILES=()
declare -A REFERENCE_SPECIES=()
declare -A OUTPUT_FORMAT_SET=()

usage() {
    printf 'GPID version: %s\n\n' "$GPID_VERSION"
    cat <<'EOF'
Usage: gpid identify -i <sample directory> -r <reference directory> [-g <gene performance file> -t <thresholds file> -c <confidence support file>] [-s <species groups file>] [-o <output directory>] [-f <output formats>] [--bypass_calibration] [--overwrite_outputs] [--remove_intermediates]

Required:
  -i  Sample directory containing one FASTA file per gene for the sample to identify
  -r  Reference directory containing one FASTA file per gene and the corresponding BLAST databases

Calibration and validation:
  -g  Gene performance CSV
      If omitted, use: calibration/calibration_gene_performance.csv
      Default file produced by: gpid calibrate genes
  -t  Filtering thresholds CSV
      If omitted, use: calibration/calibration_filtering_thresholds.csv
      Default file produced by: gpid calibrate combine
  -c  Confidence support CSV
      If omitted, use: validation/validation_confidence_support.csv
      Default file produced by: gpid validate bins

      During a run, gpid identify prints the specified filtering thresholds
      and the confidence support bins read from these files.

Optional:
  -s  Species groups CSV
  -o  Output directory (default: identification/<sample_name>)
  -f  Comma-separated output formats: csv, jpg, svg, pdf, or all (default: csv,pdf)
  --bypass_calibration  Ignore -g/-t/-c and create dummy calibration files instead
  --overwrite_outputs  Allow existing output files in the output directory to be overwritten
  --remove_intermediates  Do not save intermediate files; only keep final output files
  -h, --help  Show this help message
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

cleanup() {
    if [ -n "$temp_dir" ] && [ -d "$temp_dir" ]; then
        rm -rf "$temp_dir"
    fi
}

trap cleanup EXIT

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
            die "Sample/reference FASTA filenames must only contain the gene name before the FASTA suffix: $file"
        fi

        if [[ -n "${gene_map[$gene]+x}" ]]; then
            die "Gene '$gene' occurs more than once in directory '$dir'. Use one FASTA file per gene."
        fi

        gene_map["$gene"]="$file"
    done < <(printf '%s\n' "${files[@]}" | awk '!seen[$0]++' | sort)

    return 0
}

validate_single_sequence_fasta() {
    local fasta_file="$1"
    local context="$2"
    local line=""
    local header_count=0
    local sequence_seen=0

    while IFS= read -r line || [ -n "$line" ]; do
        line=$(trim_cr "$line")

        [ -n "$line" ] || continue

        if [[ "$line" == ">"* ]]; then
            header_count=$((header_count + 1))
            if [ "$header_count" -gt 1 ]; then
                die "$context file '$fasta_file' contains more than one FASTA record. Each gene file must contain a single sequence."
            fi
            if [ "${#line}" -le 1 ]; then
                die "$context file '$fasta_file' contains an empty FASTA header."
            fi
        else
            if [ "$header_count" -eq 0 ]; then
                die "$context file '$fasta_file' contains sequence data before the FASTA header."
            fi
            sequence_seen=1
        fi
    done < "$fasta_file"

    if [ "$header_count" -eq 0 ]; then
        die "$context file '$fasta_file' does not contain a FASTA header."
    fi

    if [ "$sequence_seen" -eq 0 ]; then
        die "$context file '$fasta_file' does not contain sequence data."
    fi
}

extract_reference_species() {
    local fasta_file="$1"
    local line=""
    local header=""
    local species=""

    while IFS= read -r line || [ -n "$line" ]; do
        line=$(trim_cr "$line")

        [[ "$line" == ">"* ]] || continue
        header=${line#>}

        if [[ "$header" =~ ^([^_]+_[^_]+) ]]; then
            species="${BASH_REMATCH[1]}"
            REFERENCE_SPECIES["$species"]=1
        fi
    done < "$fasta_file"
}

csv_header() {
    local file="$1"
    awk 'NF { gsub(/\r$/, "", $0); print; exit }' "$file"
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
    header=$(csv_header "$file")
    [ -n "$header" ] || die "CSV file is empty: $file"
    [[ "$header" == *,* ]] || die "CSV file does not appear to be comma-separated: $file"
}

validate_gene_performance_file() {
    local file="$1"
    local status=0

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
        }
        END {
            if (found == 0) {
                exit 12
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
        12) die "Gene performance file is empty: $file" ;;
        *) die "Unable to validate gene performance file: $file" ;;
    esac
}

validate_thresholds_file() {
    local file="$1"
    local expected_header="min_similarity,min_length,max_gapopens,max_mismatches,max_evalue,min_bitscore,min_gene_performance,min_parliament_size"
    local status=0

    require_csv_extension "$file"
    validate_csv_has_commas "$file"

    if awk -F',' -v expected_header="$expected_header" '
        BEGIN { row = 0 }
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
                    if ($i == "" || $i !~ /^[0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?$/) {
                        exit 12
                    }
                }
                next
            }

            exit 13
        }
        END {
            if (row == 0) {
                exit 14
            }
            if (row == 1) {
                exit 15
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
        12) die "Filtering thresholds file contains a non-numeric threshold value: $file" ;;
        13) die "Filtering thresholds file must contain only a header row and one row of thresholds: $file" ;;
        14) die "Filtering thresholds file is empty: $file" ;;
        15) die "Filtering thresholds file is missing the thresholds row: $file" ;;
        *) die "Unable to validate filtering thresholds file: $file" ;;
    esac
}

validate_confidence_support_file() {
    local file="$1"
    local status=0

    require_csv_extension "$file"
    validate_csv_has_commas "$file"

    if awk '
        BEGIN {
            row = 0
            status = 0
            expected_header = "range_support,probability_correct,probability_close,probability_wrong"
        }
        NF {
            gsub(/\r$/, "", $0)
            row++

            if (row == 1) {
                header = $0
                gsub(/"/, "", header)
                if (header != expected_header) {
                    status = 10
                }
                next
            }

            if (status != 0) {
                next
            }

            if ($0 !~ /^"[[(][[:space:]]*[0-9]+([.][0-9]+)?[[:space:]]*,[[:space:]]*[0-9]+([.][0-9]+)?[[:space:]]*[])]",(NA|[0-9]+([.][0-9]+)?),(NA|[0-9]+([.][0-9]+)?),(NA|[0-9]+([.][0-9]+)?)$/) {
                status = 11
                next
            }
        }
        END {
            if (status == 0) {
                if (row == 0) {
                    status = 12
                } else if (row == 1) {
                    status = 13
                }
            }
            exit status
        }
    ' "$file"; then
        status=0
    else
        status=$?
    fi

    case $status in
        0) log "Confidence support calibration file format check passed." ;;
        10) die "Confidence support file must use the header 'range_support,probability_correct,probability_close,probability_wrong': $file" ;;
        11) die "Confidence support file must contain rows like \"[0,20]\",NA,NA,NA with a quoted support range and numeric or NA probabilities: $file" ;;
        12) die "Confidence support file is empty: $file" ;;
        13) die "Confidence support file must contain at least one confidence bin: $file" ;;
        *) die "Unable to validate confidence support file: $file" ;;
    esac
}

validate_species_groups_file() {
    local file="$1"
    local status=0

    require_csv_extension "$file"
    validate_csv_has_commas "$file"

    if awk -F',' '
        BEGIN { row = 0 }
        NF {
            gsub(/\r$/, "", $0)
            row++

            if (row == 1) {
                if (NF != 2 || $1 != "genus_species" || $2 != "species_group") {
                    exit 10
                }
                next
            }

            if (NF != 2 || $1 == "" || $2 == "") {
                exit 11
            }

            if ($1 !~ /^[A-Z][A-Za-z]+_[a-z][A-Za-z]+$/) {
                exit 12
            }

            if ($2 !~ /^[A-Za-z0-9_]+$/) {
                exit 13
            }
        }
        END {
            if (row == 0) {
                exit 14
            }
            if (row == 1) {
                exit 15
            }
        }
    ' "$file"; then
        status=0
    else
        status=$?
    fi

    case $status in
        0) log "Species groups file format check passed." ;;
        10) die "Species groups file must use the header 'genus_species,species_group': $file" ;;
        11) die "Species groups file must contain exactly two populated columns in every row: $file" ;;
        12) die "Species groups file must use species names in the format Genus_species in the first column: $file" ;;
        13) die "Species group names may only contain letters, numbers and underscores: $file" ;;
        14) die "Species groups file is empty: $file" ;;
        15) die "Species groups file must contain at least one species-to-group mapping: $file" ;;
        *) die "Unable to validate species groups file: $file" ;;
    esac
}

check_species_groups_match_reference() {
    local file="$1"
    local matched=0
    local species=""
    local missing_species=()

    while IFS=',' read -r species _ || [ -n "$species" ]; do
        species=$(trim_cr "$species")
        [ "$species" != "genus_species" ] || continue
        [ -n "$species" ] || continue

        if [ -n "${REFERENCE_SPECIES[$species]+x}" ]; then
            matched=1
        else
            warn "Species groups file contains species not found in the reference FASTA headers: $species"
        fi
    done < "$file"

    for species in "${!REFERENCE_SPECIES[@]}"; do
        if ! grep -Fqx "$species" < <(tail -n +2 "$file" | cut -d',' -f1 | tr -d '\r'); then
            missing_species+=( "$species" )
        fi
    done

    if [ "$matched" -eq 0 ]; then
        die "None of the species names in the species groups file match the species names found in the reference FASTA headers."
    fi

    if [ "${#missing_species[@]}" -gt 0 ]; then
        while IFS= read -r species; do
            [ -n "$species" ] || continue
            warn "No species group was provided for reference species: $species"
        done < <(printf '%s\n' "${missing_species[@]}" | sort)
    else
        log "Species groups file matches all reference species."
    fi
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

lookup_threshold_value() {
    local file="$1"
    local column_name="$2"

    awk -F',' -v column_name="$column_name" '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                gsub(/\r$/, "", $i)
                if ($i == column_name) {
                    column = i
                }
            }
            next
        }
        NR == 2 {
            gsub(/\r$/, "", $column)
            print $column
            exit
        }
    ' "$file"
}

print_filtering_threshold_summary() {
    local file="$1"
    local values=""
    local min_similarity=""
    local min_length=""
    local max_gapopens=""
    local max_mismatches=""
    local max_evalue=""
    local min_bitscore=""
    local min_gene_performance=""
    local min_parliament_size=""

    values=$(awk -F',' '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                gsub(/\r$/, "", $i)
                columns[$i] = i
            }
            next
        }
        NR == 2 {
            print $columns["min_similarity"] "," \
                $columns["min_length"] "," \
                $columns["max_gapopens"] "," \
                $columns["max_mismatches"] "," \
                $columns["max_evalue"] "," \
                $columns["min_bitscore"] "," \
                $columns["min_gene_performance"] "," \
                $columns["min_parliament_size"]
            exit
        }
    ' "$file")

    IFS=',' read -r min_similarity min_length max_gapopens max_mismatches max_evalue min_bitscore min_gene_performance min_parliament_size <<< "$values"

    log "Specified filtering thresholds:"
    log "Minimum alignment similarity: $min_similarity"
    log "Minimum alignment length: $min_length"
    log "Maximum number of gap openings: $max_gapopens"
    log "Maximum number of alignment mismatches: $max_mismatches"
    log "Maximum E-value: $max_evalue"
    log "Minimum Bit-score: $min_bitscore"
    log "Minimum gene performance: $min_gene_performance"
    log "Minimum gene parliament size: $min_parliament_size"
}

print_confidence_support_summary() {
    local file="$1"

    log "Confidence support bins:"
    awk '
        NR == 1 {
            next
        }
        NF {
            line = $0
            gsub(/\r$/, "", line)
            range_support = line
            sub(/^"/, "", range_support)
            sub(/",.*$/, "", range_support)

            probabilities = line
            sub(/^"[^"]+",/, "", probabilities)
            split(probabilities, values, ",")

            printf "  %s: probability_correct=%s, probability_close=%s, probability_wrong=%s\n", range_support, values[1], values[2], values[3]
        }
    ' "$file"
}

normalise_output_formats() {
    local raw_formats="$1"
    local token=""
    local cleaned=""
    local all_requested=0
    local tokens=()

    OUTPUT_FORMAT_SET=()

    IFS=',' read -r -a tokens <<< "$raw_formats"
    [ "${#tokens[@]}" -gt 0 ] || die "At least one output format must be specified with -f."

    for token in "${tokens[@]}"; do
        cleaned=${token,,}
        cleaned=${cleaned//[[:space:]]/}

        [ -n "$cleaned" ] || die "Output formats in -f must not contain empty entries."

        case "$cleaned" in
            csv|jpg|svg|pdf)
                OUTPUT_FORMAT_SET["$cleaned"]=1
                ;;
            all)
                all_requested=1
                ;;
            *)
                die "Unsupported output format '$token'. Allowed values are csv, jpg, svg, pdf and all."
                ;;
        esac
    done

    if [ "$all_requested" -eq 1 ]; then
        OUTPUT_FORMAT_SET=(
            ["csv"]=1
            ["jpg"]=1
            ["svg"]=1
            ["pdf"]=1
        )
    fi
}

output_file_requested() {
    local format="$1"
    [ -n "${OUTPUT_FORMAT_SET[$format]+x}" ]
}

ensure_temp_dir() {
    if [ -z "$temp_dir" ] || [ ! -d "$temp_dir" ]; then
        temp_dir=$(mktemp -d "${TMPDIR:-/tmp}/gpid-identification.XXXXXX")
    fi
}

planned_output_files() {
    local files=()

    if [ "$remove_intermediates" -eq 0 ]; then
        files+=(
            "$output_dir/genelist_high_performance.txt"
            "$output_dir/${sample}_blast.tsv"
        )
    fi

    if output_file_requested csv; then
        files+=( "$output_dir/${sample}_gpid.csv" )
    fi

    if output_file_requested jpg; then
        files+=( "$output_dir/${sample}_gpid.jpg" )
    fi

    if output_file_requested pdf; then
        files+=( "$output_dir/${sample}_gpid.pdf" )
    fi

    if output_file_requested svg; then
        files+=( "$output_dir/${sample}_gpid.svg" )
    fi

    printf '%s\n' "${files[@]}"
}

check_existing_output_files() {
    local existing_files=()
    local path=""

    while IFS= read -r path; do
        [ -n "$path" ] || continue

        if [ -e "$path" ]; then
            existing_files+=( "$path" )
        fi
    done < <(planned_output_files)

    if [ "${#existing_files[@]}" -eq 0 ]; then
        return 0
    fi

    if [ "$overwrite_output_dir" -eq 1 ]; then
        warn "Existing output files will be overwritten in: $output_dir"
        while IFS= read -r path; do
            [ -n "$path" ] || continue
            warn "  $path"
        done < <(printf '%s\n' "${existing_files[@]}" | sort)
        return 0
    fi

    while IFS= read -r path; do
        [ -n "$path" ] || continue
        warn "Existing output file detected: $path"
    done < <(printf '%s\n' "${existing_files[@]}" | sort)

    die "Output files already exist in '$output_dir'. Use --overwrite_outputs to overwrite the existing outputs."
}

create_dummy_calibration_files() {
    local gene=""

    ensure_temp_dir

    gene_performance_file="$temp_dir/dummy_gene_performance.csv"
    thresholds_file="$temp_dir/dummy_filtering_thresholds.csv"
    confidence_support_file="$temp_dir/dummy_confidence_support.csv"

    {
        printf 'gene,performance\n'
        while IFS= read -r gene; do
            printf '%s,100\n' "$gene"
        done < <(printf '%s\n' "${!REFERENCE_GENE_FILES[@]}" | sort)
    } > "$gene_performance_file"

    cat > "$thresholds_file" <<'EOF'
min_similarity,min_length,max_gapopens,max_mismatches,max_evalue,min_bitscore,min_gene_performance,min_parliament_size
0,0,99999,99999,99999,0,0,0
EOF

    cat > "$confidence_support_file" <<'EOF'
range_support,probability_correct,probability_close,probability_wrong
"[0,100]",NA,NA,NA
EOF

    warn "Bypassing method calibration will likely reduce identification accuracy considerably. This may still be justified when no suitable test dataset is available or for an initial exploratory analysis."
    log "Dummy calibration files created for calibration bypass:"
    log "  $gene_performance_file"
    log "  $thresholds_file"
    log "  $confidence_support_file"
}

while [ "$#" -gt 0 ]; do
    case "$1" in
        -i)
            [ "$#" -ge 2 ] || die "Option -i requires an argument."
            sample_dir="$2"
            shift 2
            ;;
        -r)
            [ "$#" -ge 2 ] || die "Option -r requires an argument."
            reference_dir="$2"
            shift 2
            ;;
        -g)
            [ "$#" -ge 2 ] || die "Option -g requires an argument."
            gene_performance_file="$2"
            gene_performance_file_specified=1
            shift 2
            ;;
        -t)
            [ "$#" -ge 2 ] || die "Option -t requires an argument."
            thresholds_file="$2"
            thresholds_file_specified=1
            shift 2
            ;;
        -c)
            [ "$#" -ge 2 ] || die "Option -c requires an argument."
            confidence_support_file="$2"
            confidence_support_file_specified=1
            shift 2
            ;;
        -s)
            [ "$#" -ge 2 ] || die "Option -s requires an argument."
            species_groups_file="$2"
            shift 2
            ;;
        -o)
            [ "$#" -ge 2 ] || die "Option -o requires an argument."
            output_dir="$2"
            shift 2
            ;;
        -f)
            [ "$#" -ge 2 ] || die "Option -f requires an argument."
            output_formats_raw="$2"
            shift 2
            ;;
        --bypass_calibration)
            bypass_calibration=1
            shift
            ;;
        --overwrite_outputs)
            overwrite_output_dir=1
            shift
            ;;
        --remove_intermediates)
            remove_intermediates=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "Unknown option: $1. Use -h for help."
            ;;
    esac
done

if [ -z "$sample_dir" ] && [ -z "$reference_dir" ] && [ -z "$gene_performance_file" ] && [ -z "$thresholds_file" ] && [ -z "$confidence_support_file" ] && [ -z "$species_groups_file" ] && [ -z "$output_dir" ] && [ "$output_formats_raw" = "csv,pdf" ] && [ "$bypass_calibration" -eq 0 ] && [ "$overwrite_output_dir" -eq 0 ] && [ "$remove_intermediates" -eq 0 ]; then
    usage
    exit 1
fi

if [ "$bypass_calibration" -eq 0 ]; then
    if [ -z "$gene_performance_file" ]; then
        gene_performance_file="$DEFAULT_GENE_PERFORMANCE_FILE"
    fi

    if [ -z "$thresholds_file" ]; then
        thresholds_file="$DEFAULT_THRESHOLDS_FILE"
    fi

    if [ -z "$confidence_support_file" ]; then
        confidence_support_file="$DEFAULT_CONFIDENCE_SUPPORT_FILE"
    fi
fi

input_path_errors=0

if [ -n "$sample_dir" ] && [ ! -d "$sample_dir" ]; then
    printf 'Error: Sample directory not found: %s\n' "$sample_dir" >&2
    input_path_errors=1
fi

if [ -n "$reference_dir" ] && [ ! -d "$reference_dir" ]; then
    printf 'Error: Reference directory not found: %s\n' "$reference_dir" >&2
    input_path_errors=1
fi

if [ -n "$gene_performance_file" ] && [ "$bypass_calibration" -eq 0 ] && [ ! -f "$gene_performance_file" ]; then
    printf 'Error: Gene performance file not found: %s\n' "$gene_performance_file" >&2
    input_path_errors=1
fi

if [ -n "$thresholds_file" ] && [ "$bypass_calibration" -eq 0 ] && [ ! -f "$thresholds_file" ]; then
    printf 'Error: Filtering thresholds file not found: %s\n' "$thresholds_file" >&2
    input_path_errors=1
fi

if [ -n "$confidence_support_file" ] && [ "$bypass_calibration" -eq 0 ] && [ ! -f "$confidence_support_file" ]; then
    printf 'Error: Confidence support file not found: %s\n' "$confidence_support_file" >&2
    input_path_errors=1
fi

if [ -n "$species_groups_file" ] && [ ! -f "$species_groups_file" ]; then
    printf 'Error: Species groups file not found: %s\n' "$species_groups_file" >&2
    input_path_errors=1
fi

if [ -z "$sample_dir" ] || [ -z "$reference_dir" ]; then
    usage
    die "Both -i <sample directory> and -r <reference directory> are required."
fi

if [ "$input_path_errors" -ne 0 ]; then
    exit 1
fi

sample_dir=${sample_dir%/}
reference_dir=${reference_dir%/}

if [ "$bypass_calibration" -eq 1 ]; then
    [ "$gene_performance_file_specified" -eq 0 ] || warn "Ignoring gene performance file because --bypass_calibration was specified: $gene_performance_file"
    [ "$thresholds_file_specified" -eq 0 ] || warn "Ignoring filtering thresholds file because --bypass_calibration was specified: $thresholds_file"
    [ "$confidence_support_file_specified" -eq 0 ] || warn "Ignoring confidence support file because --bypass_calibration was specified: $confidence_support_file"
fi

sample=$(basename "$sample_dir")

normalise_output_formats "$output_formats_raw"

if [ -z "$output_dir" ]; then
    output_dir="identification/$sample"
fi

output_dir=${output_dir%/}

if [ -e "$output_dir" ] && [ ! -d "$output_dir" ]; then
    die "Output path exists but is not a directory: $output_dir"
fi

check_existing_output_files

mkdir -p "$output_dir" || die "Could not create output directory: $output_dir"

log "Step 1/6: Checking input directories and output settings."
log "Output directory: $output_dir"
if [ "$remove_intermediates" -eq 1 ]; then
    log "Intermediate files will be removed after the run."
fi

collect_gene_files "$sample_dir" SAMPLE_GENE_FILES || die "No FASTA gene files found in sample directory. Expected files ending in .FNA, .fasta or .fa."
collect_gene_files "$reference_dir" REFERENCE_GENE_FILES || die "No FASTA gene files found in reference directory. Expected files ending in .FNA, .fasta or .fa."
log "Genes in sample directory before filtering: ${#SAMPLE_GENE_FILES[@]}"

for gene in "${!SAMPLE_GENE_FILES[@]}"; do
    validate_single_sequence_fasta "${SAMPLE_GENE_FILES[$gene]}" "Sample"
done

for gene in "${!REFERENCE_GENE_FILES[@]}"; do
    extract_reference_species "${REFERENCE_GENE_FILES[$gene]}"
done

if [ "${#REFERENCE_SPECIES[@]}" -eq 0 ]; then
    die "No species names could be extracted from the FASTA headers in the reference directory."
fi

log "Step 2/6: Checking sample genes against the reference dataset."

matched_genes=()
missing_reference_genes=()

for gene in "${!SAMPLE_GENE_FILES[@]}"; do
    if [ -n "${REFERENCE_GENE_FILES[$gene]+x}" ]; then
        matched_genes+=( "$gene" )
    else
        missing_reference_genes+=( "$gene" )
        warn "Gene found in the sample directory but not in the reference directory: $gene"
    fi
done

if [ "${#matched_genes[@]}" -eq 0 ]; then
    die "None of the genes in the sample directory were found in the reference directory."
fi

if [ "${#missing_reference_genes[@]}" -eq 0 ]; then
    log "All sample gene names were found in the reference directory."
fi

missing_databases=()
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

log "Reference directory contains BLAST databases for all sample genes found in the reference."

if [ "$bypass_calibration" -eq 1 ]; then
    log "Step 3/6: Bypassing calibration and validation inputs."
    create_dummy_calibration_files
else
    log "Step 3/6: Checking calibration and validation inputs."
    validate_gene_performance_file "$gene_performance_file"
    validate_thresholds_file "$thresholds_file"
    validate_confidence_support_file "$confidence_support_file"
fi

print_filtering_threshold_summary "$thresholds_file"
print_confidence_support_summary "$confidence_support_file"

if [ -n "$species_groups_file" ]; then
    log "Checking species groups file against reference species."
    validate_species_groups_file "$species_groups_file"
    check_species_groups_match_reference "$species_groups_file"
else
    log "No species groups file provided. Species groups will be derived from the genus names in the BLAST results."
fi

gene_performance_threshold=$(lookup_threshold_value "$thresholds_file" "min_gene_performance")
[ -n "$gene_performance_threshold" ] || die "Unable to read min_gene_performance from $thresholds_file"

if [ "$remove_intermediates" -eq 1 ]; then
    ensure_temp_dir
    genelist_file="$temp_dir/genelist_high_performance.txt"
    blast_file="$temp_dir/${sample}_blast.tsv"
else
    genelist_file="$output_dir/genelist_high_performance.txt"
    blast_file="$output_dir/${sample}_blast.tsv"
fi

log "Step 4/6: Selecting high-performance genes for identification."

genelist_high_performance=$(
    awk -F',' -v threshold="$gene_performance_threshold" '
        NR > 1 {
            gsub(/\r$/, "", $1)
            gsub(/\r$/, "", $2)
            if ($2 + 0 >= threshold) {
                print $1
            }
        }
    ' "$gene_performance_file"
)

{
    comm -12 \
        <(printf '%s\n' "$genelist_high_performance" | awk 'NF' | sort -u) \
        <(printf '%s\n' "${matched_genes[@]}" | sort -u)
} > "$genelist_file"

high_performance_gene_count=$(awk 'NF { count++ } END { print count + 0 }' "$genelist_file")
log "Genes after gene performance filtering: $high_performance_gene_count"

if [ ! -s "$genelist_file" ]; then
    warn "No gene passed gene performance threshold."

    if output_file_requested csv; then
        cat > "$output_dir/${sample}_gpid.csv" <<EOF
Sample,Rank,Identification,Species_group,Support_pct,Support_count,Parliament_size,Data_checks,ID_correct_pct,ID_close_pct,ID_wrong_pct
${sample},NA,NA,NA,NA,NA,0,FAILED_1,,,
EOF
        log "Minimal gene parliament table written:"
        log "$output_dir/${sample}_gpid.csv"
    else
        log "CSV output was not requested, so no gene parliament table was written."
    fi

    log "No gene parliament figure produced."
    log "Terminating."
    exit 1
fi

if [ "$remove_intermediates" -eq 0 ]; then
    log "High-performance gene list written:"
    log "$genelist_file"
else
    log "High-performance gene list created temporarily for internal use only."
fi
log ""

command -v blastn >/dev/null 2>&1 || die "blastn command not found in PATH."
command -v Rscript >/dev/null 2>&1 || die "Rscript command not found in PATH."

log "Step 5/6: Matching sample genes against reference BLAST databases."

while IFS= read -r gene; do
    [ -n "$gene" ] || continue

    blastn \
        -query "${SAMPLE_GENE_FILES[$gene]}" \
        -db "${REFERENCE_GENE_FILES[$gene]}" \
        -task megablast \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen evalue bitscore" \
        -max_target_seqs 1000000 |
        sort -t $'\t' -k8,8rn -k2,2R |
        head -1 |
        awk -v genename="$gene" '{print genename "\t" $0}'
done < "$genelist_file" > "$blast_file"

sed -i '1i gene\tquery\ttarget\tpident\tlength\tmismatch\tgapopen\tevalue\tbitscore' "$blast_file"

if [ "$remove_intermediates" -eq 0 ]; then
    log "Top-match BLAST table written:"
    log "$blast_file"
else
    log "Top-match BLAST table created temporarily for internal use only."
fi
log ""

log "Step 6/6: Computing Gene Parliament and writing requested outputs."

Rscript \
    "$SCRIPT_DIR/parliament.R" \
    "$sample_dir" \
    "$blast_file" \
    "$gene_performance_file" \
    "$thresholds_file" \
    "$confidence_support_file" \
    "$species_groups_file" \
    "$output_dir" \
    "$output_formats_raw" \
    "$overwrite_output_dir"

if output_file_requested csv; then
    log "Gene Parliament table written:"
    log "$output_dir/${sample}_gpid.csv"
fi

if output_file_requested jpg; then
    log "Gene Parliament JPG figure written:"
    log "$output_dir/${sample}_gpid.jpg"
fi

if output_file_requested pdf; then
    log "Gene Parliament PDF figure written:"
    log "$output_dir/${sample}_gpid.pdf"
fi

if output_file_requested svg; then
    log "Gene Parliament SVG figure written:"
    log "$output_dir/${sample}_gpid.svg"
fi

log ""
log "Identification complete."
