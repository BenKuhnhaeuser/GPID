# --- Purpose ---
# Shared shell CSV reader: validate gene performance, or emit genes meeting the filtering threshold.
# Shared CSV validation and gene filtering. Additional columns may occur in any
# order, including quoted fields containing commas, quotes or line breaks.
# fail(): Report a malformed input with context and stop further processing.
function fail(message) {
    print "Error: Gene performance file " FILENAME ": " message > "/dev/stderr"
    failed = 1
    exit 1
}
# trim(): Remove surrounding whitespace before comparing names or parsing values.
function trim(value) {
    sub(/^[[:space:]]+/, "", value)
    sub(/[[:space:]]+$/, "", value)
    return value
}
# csv_fields(): Parse one CSV record, returning a field count, -1 for unfinished quoting, or -2 for malformed quoting.
function csv_fields(line, fields,    i, c, value, quoted, closed, n) {
    n = 1
    value = ""
    quoted = closed = 0
    for (i = 1; i <= length(line); i++) {
        c = substr(line, i, 1)
        if (quoted) {
            if (c == "\"") {
                if (substr(line, i + 1, 1) == "\"") { value = value "\""; i++ }
                else { quoted = 0; closed = 1 }
            } else value = value c
        } else if (c == ",") {
            fields[n++] = trim(value)
            value = ""
            closed = 0
        } else if (c == "\"") {
            if (closed || trim(value) != "") return -2
            value = ""
            quoted = 1
        } else {
            if (closed && c !~ /[[:space:]]/) return -2
            value = value c
        }
    }
    if (quoted) return -1
    fields[n] = trim(value)
    return n
}
# --- Read complete CSV records, including quoted multiline fields ---
{
    sub(/\r$/, "", $0)
    if (!pending && $0 ~ /^[[:space:]]*$/) next
    # A quoted field may span physical lines; validate only complete CSV records.
    record = pending ? record "\n" $0 : $0
    count = csv_fields(record, fields)
    if (count == -1) { pending = 1; next }
    if (count == -2) fail("malformed CSV quoting near line " NR ".")
    pending = 0
    if (!header_seen) {
        columns = count
        for (i = 1; i <= count; i++) {
            if (fields[i] == "gene") { gene_column = i; gene_headers++ }
            if (fields[i] == "performance") { performance_column = i; performance_headers++ }
        }
        if (gene_headers != 1 || performance_headers != 1)
            fail("must contain exactly one column named 'gene' and one named 'performance'; additional columns are allowed.")
        header_seen = 1
        next
    }
    if (count != columns) fail("row ending at line " NR " has " count " fields, but the header has " columns ".")
    gene = fields[gene_column]
    value = fields[performance_column]
    if (gene == "") fail("contains an empty gene name near line " NR ".")
    if (seen[gene]++) fail("contains duplicated gene name '" gene "'.")
    if (value == "NA") {
        na_genes = na_genes (na_genes == "" ? "" : ", ") gene
        value = 0
    } else if (value !~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/) {
        fail("performance for gene '" gene "' must be numeric or NA; found '" value "'.")
    }
    if (value + 0 < 0 || value + 0 > 100) fail("performance for gene '" gene "' must be between 0 and 100; found '" value "'.")
    genes[++rows] = gene
    performance[gene] = value + 0
}
# --- Finish validation and optionally emit genes above the threshold ---
END {
    # AWK executes END even after exit; preserve prior failures and emit no genes.
    if (failed) exit 1
    if (pending) fail("contains an unterminated quoted field.")
    if (!header_seen || !rows) fail("must contain a header and at least one gene performance row.")
    if (na_genes != "" && mode != "filter")
        print "Warning: NA performance for gene(s): " na_genes ". Treating NA as 0 for filtering (" FILENAME ")." > "/dev/stderr"
    if (mode == "filter") {
        for (i = 1; i <= rows; i++)
            if (performance[genes[i]] >= threshold) print genes[i]
    }
}
