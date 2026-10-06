# --- Purpose ---
# Embed the calibration parameter table so no external template CSV is required.

# These are the former template contents, retained as reference/example values.
# They are not selected automatically: users choose thresholds from their data.
calibration_parameter_table <- data.frame(
  parameter = c("min_similarity", "min_length", "max_gapopens", "max_mismatches",
                "max_evalue", "min_bitscore", "min_gene_performance", "min_parliament_size"),
  value = c(98, 100, 1, 5, 1e-60, 200, 30, 10),
  stringsAsFactors = FALSE
)

# calibration_threshold_rows(): Create editable rows for one step, leaving values unset for user selection.
calibration_threshold_rows <- function(parameters) {
  indices <- match(parameters, calibration_parameter_table$parameter)
  if (anyNA(indices)) stop("Unknown calibration parameter requested.", call. = FALSE)
  rows <- calibration_parameter_table[indices, , drop = FALSE]
  rows$value <- NA_real_
  rownames(rows) <- NULL
  rows
}
