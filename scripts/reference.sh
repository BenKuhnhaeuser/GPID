#!/bin/bash
# --- Purpose ---
# Reference workflow: check FASTA headers and build any missing BLAST databases.

#################################################################
# Reference directory preparation                               #
# Benedikt Kuhnhaeuser                                          #
# Royal Botanic Gardens, Kew                                    #
# 2026                                                          #
#################################################################

set -euo pipefail

# --- Installation paths and version ---
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
VERSION_FILE="${GPID_VERSION_FILE:-$PROJECT_DIR/VERSION}"
if [ -z "${GPID_VERSION:-}" ] && [ -f "$VERSION_FILE" ]; then
    IFS= read -r GPID_VERSION < "$VERSION_FILE" || GPID_VERSION=""
    GPID_VERSION=${GPID_VERSION%$'\r'}
fi
GPID_VERSION="${GPID_VERSION:-unknown}"

# usage(): Print command usage, supported options, and output locations.
usage() {
    printf 'GPID version: %s\n\n' "$GPID_VERSION"
    cat <<'EOF'
Usage: gpid reference -r <reference directory>

Check reference FASTA files and build missing BLAST databases with makeblastdb.

Required:
  -r  Reference directory containing FASTA files (.fna, .fasta, .fa; case-insensitive)
      Headers must start with Genus_species, contain no whitespace, and be unique
      within each file. Maximum header length: 50 characters, excluding >.

Optional:
  -h  Show this help message

Outputs:
  BLAST database files beside each reference FASTA file; existing databases are reused.
EOF
}

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

# --- Parse CLI arguments and check the reference directory ---
reference_dir=""

if [ "$#" -eq 0 ]; then
    usage
    exit 1
fi

while getopts ":r:h" opt; do
    case "$opt" in
        r) reference_dir="$OPTARG" ;;
        h)
            usage
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
[ -d "$reference_dir" ] || die "Reference directory not found: $reference_dir"

reference_dir=${reference_dir%/}

# --- Expected BLAST database components ---
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

FASTA_FILES=()

# collect_fasta_files(): Collect supported FASTA extensions case-insensitively without duplicate paths.
collect_fasta_files() {
    local dir="$1"
    local matches=()

    shopt -s nullglob nocaseglob
    matches=( "$dir"/*.fna "$dir"/*.fasta "$dir"/*.fa )
    shopt -u nullglob nocaseglob

    if [ "${#matches[@]}" -eq 0 ]; then
        return 1
    fi

    mapfile -t FASTA_FILES < <(printf '%s\n' "${matches[@]}" | awk '!seen[$0]++' | sort)
}

# blast_db_complete(): Return success only when every expected BLAST database component exists.
blast_db_complete() {
    local fasta_file="$1"
    local suffix=""

    for suffix in "${BLAST_SUFFIXES[@]}"; do
        if [ ! -f "${fasta_file}${suffix}" ]; then
            return 1
        fi
    done

    return 0
}

# validate_fasta_file(): Check reference FASTA records, unique species-prefixed names, and the 50-character header limit.
validate_fasta_file() {
    local fasta_file="$1"
    local line=""
    local line_number=0
    local header=""
    local file_failed=0
    local header_seen=0
    local -a special_headers=()
    declare -A seen_headers=()

    while IFS= read -r line || [ -n "$line" ]; do
        line_number=$((line_number + 1))
        line=${line%$'\r'}

        [ -n "$line" ] || continue

        if [[ "$line" == ">"* ]]; then
            header="${line#>}"
            header_seen=1

            if [ -z "$header" ]; then
                printf 'Error: %s:%s contains an empty FASTA header.\n' "$fasta_file" "$line_number" >&2
                file_failed=1
                continue
            fi

            if [ "${#header}" -gt 50 ]; then
                printf 'Error: FASTA header in %s:%s is %s characters long; the maximum for BLAST reference databases is 50 characters (excluding >). Shorten this header: %s\n' "$fasta_file" "$line_number" "${#header}" "$header" >&2
                file_failed=1
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
        fi
    done < "$fasta_file"

    if [ "$header_seen" -eq 0 ]; then
        printf 'Error: No FASTA headers found in %s.\n' "$fasta_file" >&2
        file_failed=1
    fi

    if [ "${#special_headers[@]}" -gt 0 ]; then
        warn "Special characters were found in sequence names in $fasta_file. This could cause issues with downstream data processing, and sequence names should only include underscores as separators."
        printf '%s\n' "${special_headers[@]}" >&2
    fi

    return "$file_failed"
}

# --- Run reference checks and build missing databases ---
log "This script prepares a reference directory for GPID by validating reference FASTA files and building any missing BLAST databases."
log "Step 1/4: Checking for FASTA files in $reference_dir"

if ! collect_fasta_files "$reference_dir"; then
    die "No FASTA files found. Expected files ending in .FNA, .fasta or .fa."
fi

log "Found ${#FASTA_FILES[@]} FASTA file(s)."

log "Step 2/4: Checking FASTA header format"

validation_failed=0

for fasta_file in "${FASTA_FILES[@]}"; do
    log "Checking $(basename "$fasta_file")"
    if ! validate_fasta_file "$fasta_file"; then
        validation_failed=1
    fi
done

if [ "$validation_failed" -ne 0 ]; then
    die "Reference directory validation failed. Please fix the FASTA files and run the script again."
fi

log "All FASTA header checks passed."

log "Step 3/4: Checking whether BLAST reference databases already exist"

missing_db_files=()

for fasta_file in "${FASTA_FILES[@]}"; do
    if blast_db_complete "$fasta_file"; then
        log "BLAST database files already present for $(basename "$fasta_file")"
    else
        missing_db_files+=( "$fasta_file" )
        log "BLAST database files missing for $(basename "$fasta_file")"
    fi
done

if [ "${#missing_db_files[@]}" -eq 0 ]; then
    log "Step 4/4: Building BLAST reference databases"
    log "BLAST reference databases already exist for all FASTA files. Nothing to do."
    exit 0
fi

log "Step 4/4: Building BLAST reference databases"

command -v makeblastdb >/dev/null 2>&1 || die "makeblastdb command not found on PATH."

for fasta_file in "${missing_db_files[@]}"; do
    log "Building BLAST database for $(basename "$fasta_file")"
    makeblastdb -in "$fasta_file" -parse_seqids -dbtype nucl

    if blast_db_complete "$fasta_file"; then
        log "Successfully created BLAST database files for $(basename "$fasta_file")"
    else
        die "BLAST database creation appears incomplete for $(basename "$fasta_file")."
    fi
done

log "Reference directory preparation completed successfully."
