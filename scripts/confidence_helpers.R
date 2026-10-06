# --- Purpose ---
# Shared confidence label formatting and output-table conversion; plotting data stays unchanged.
# Plotmath supplies a true minus sign even with the standard PDF device's fonts.
# Keep the original interval strings in data files for downstream range parsing.
# confidence_bin_labels(): Convert interval labels to readable percentage expressions without changing stored CSV intervals.
confidence_bin_labels <- function(labels) {
  # Callback: Translate one stored interval into a PDF-compatible percentage label.
  as.expression(lapply(as.character(labels), function(label) {
    if (is.na(label)) return(NA_character_)
    bounds <- strsplit(gsub("[[:space:]]", "", label), ",", fixed = TRUE)[[1]]
    lower <- sub("^\\(", ">", sub("^\\[", "", bounds[[1]]))
    upper <- paste0(sub("\\]$", "", bounds[[2]]), "%")
    bquote(.(lower) - .(upper))
  }))
}

# confidence_output_table(): Copy export columns and make close cumulative from counts while leaving plot probabilities unchanged.
confidence_output_table <- function(support) {
  output <- support[, c("range_support", "probability_correct", "probability_close", "probability_wrong")]
  # The CSV reports correct OR close; the plot retains mutually exclusive classes.
  # Combine counts before rounding, avoiding accumulated rounding error.
  output$probability_close <- round((support$count_correct + support$count_close) / support$count_all * 100, 2)
  output
}
