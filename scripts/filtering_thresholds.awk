# --- Purpose ---
# Validate current parameter/value and legacy single-row threshold CSVs.
# Optional modes return a named value (lookup) or values in canonical order (summary).

# clean(): Normalize whitespace, CRLF endings, and optional CSV field quotes.
function clean(value) {
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
    if (value ~ /^".*"$/) value = substr(value, 2, length(value) - 2)
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
    return value
}

# fail(): Preserve a validation error when AWK subsequently executes END.
function fail(code) { failed = code; exit code }

# record(): Check a named threshold and retain it for downstream lookup.
function record(name, value, maximum) {
    if (!(name in required)) fail(12)
    if (seen[name]++) fail(13)
    if (value !~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/) fail(14)
    if (enforce_ranges) {
        maximum = (name == "min_similarity" || name == "max_evalue" || name == "min_gene_performance") ? 100 : 99999
        if (value + 0 > maximum) fail(16)
    }
    values[name] = value
}

BEGIN {
    FS = ","
    split("min_similarity min_length max_gapopens max_mismatches max_evalue min_bitscore min_gene_performance min_parliament_size", names, " ")
    for (i = 1; i <= 8; i++) required[names[i]] = 1
}

/^[[:space:]]*$/ { next }
{
    row++
    if (row == 1) {
        if (NF == 2 && clean($1) == "parameter" && clean($2) == "value") {
            format = "long"
        } else if (NF == 8) {
            format = "wide"
            for (i = 1; i <= 8; i++) {
                headers[i] = clean($i)
                if (!(headers[i] in required)) fail(12)
                if (header_seen[headers[i]]++) fail(13)
            }
        } else fail(10)
        next
    }
    if (format == "long") {
        if (NF != 2) fail(11)
        record(clean($1), clean($2))
    } else {
        if (row > 2) fail(17)
        if (NF != 8) fail(11)
        for (i = 1; i <= 8; i++) record(headers[i], clean($i))
    }
}

END {
    if (failed) exit failed
    for (i = 1; i <= 8; i++) if (!seen[names[i]]) exit 15
    if (mode == "lookup") print values[parameter]
    if (mode == "summary") {
        for (i = 1; i <= 8; i++) printf "%s%s", values[names[i]], (i == 8 ? "\n" : ",")
    }
}
