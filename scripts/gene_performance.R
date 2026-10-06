# --- Purpose ---
# Shared R CSV reader: validate gene performance and normalize explicit NA values to zero.
# Base-R reader shared by confidence estimation and identification.
# Preserve the source CSV; NA performance is represented as zero in memory.
# read_gene_performance(): Read required gene/performance columns, reject malformed values, and normalize explicit NA to zero.
read_gene_performance <- function(file, warn_na = TRUE) {
  # fail(): Report a malformed input with context and stop further processing.
  fail <- function(message) stop(sprintf("Gene performance file %s: %s", file, message), call. = FALSE)
  field_counts <- count.fields(file, sep = ",", quote = "\"", comment.char = "", blank.lines.skip = TRUE)
  field_counts <- field_counts[!is.na(field_counts)]
  if (length(field_counts) < 2L) fail("must contain a header and at least one gene performance row.")
  if (any(field_counts != field_counts[[1]])) fail("contains rows with a different number of fields than the header.")
  data <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE,
                   colClasses = "character", na.strings = character(),
                   strip.white = TRUE, fill = FALSE)
  names(data) <- trimws(names(data))
  if (sum(names(data) == "gene") != 1L || sum(names(data) == "performance") != 1L) {
    fail("must contain exactly one column named 'gene' and one named 'performance'; additional columns are allowed.")
  }
  data <- data[, c("gene", "performance"), drop = FALSE]
  data$gene <- trimws(data$gene)
  if (anyNA(data$gene) || any(data$gene == "")) fail("contains an empty gene name.")
  if (anyDuplicated(data$gene)) fail(paste("contains duplicated gene names:", paste(unique(data$gene[duplicated(data$gene)]), collapse = ", ")))
  raw <- trimws(data$performance)
  missing <- !is.na(raw) & raw == "NA"
  numeric_value <- grepl("^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$", raw)
  invalid <- !missing & (is.na(raw) | !numeric_value)
  if (any(invalid)) fail(paste("performance must be numeric or NA for gene(s):", paste(data$gene[invalid], collapse = ", ")))
  raw[missing] <- "0"
  data$performance <- suppressWarnings(as.numeric(raw))
  invalid <- !is.finite(data$performance) | data$performance < 0 | data$performance > 100
  if (any(invalid)) fail(paste("performance must be between 0 and 100 for gene(s):", paste(data$gene[invalid], collapse = ", ")))
  if (warn_na && any(missing)) {
    warning(sprintf("NA performance for gene(s): %s. Treating NA as 0 for filtering (%s).",
                    paste(data$gene[missing], collapse = ", "), file), call. = FALSE)
  }
  data
}
