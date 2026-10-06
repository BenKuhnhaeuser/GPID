#!/usr/bin/env Rscript
# --- Purpose ---
# Rebuild the eight calibration plots from saved CSVs and mark the selected thresholds in red.

# --- Read input paths and shared panel definitions ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: calibration_combine.R <thresholds_filtering.csv> <calibration test directory> <output.pdf>", call. = FALSE)
}

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]])
source(file.path(dirname(script_file), "calibration_plots.R"))

# --- Validate selected thresholds and index them by parameter ---
thresholds <- read.csv(args[[1]], stringsAsFactors = FALSE, check.names = FALSE)
if (!identical(names(thresholds), c("parameter", "value")) ||
    anyDuplicated(thresholds$parameter) ||
    !setequal(thresholds$parameter, calibration_plot_specs$parameter)) {
  stop("Filtering thresholds must contain the eight unique calibration parameters in parameter,value columns.", call. = FALSE)
}
values <- suppressWarnings(as.numeric(thresholds$value))
if (any(!is.finite(values)) || any(values < 0)) {
  stop("Selected thresholds must be finite, nonnegative numbers.", call. = FALSE)
}
selected <- setNames(values, thresholds$parameter)

# --- Load the saved diagnostic curves ---
curve_files <- file.path(args[[2]], calibration_plot_specs$file)
missing_files <- curve_files[!file.exists(curve_files)]
if (length(missing_files)) {
  stop("Calibration plot data are missing. Run gpid calibrate alignments, genes and parliament first. Missing file(s): ",
       paste(missing_files, collapse = ", "), call. = FALSE)
}
# Callback: Read one saved curve and check its threshold and performance columns.
curves <- lapply(curve_files, function(file) {
  data <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("threshold", "proportion_correct", "retrievability_all")
  if (!all(required %in% names(data)) || !nrow(data)) {
    stop("Calibration plot data must contain rows with threshold, proportion_correct and retrievability_all columns: ", file, call. = FALSE)
  }
  for (column in required) {
    raw <- data[[column]]
    data[[column]] <- suppressWarnings(as.numeric(raw))
    if (any(!is.na(raw) & is.na(data[[column]])) || any(is.infinite(data[[column]]))) {
      stop("Invalid numeric values in calibration plot column ", column, ": ", file, call. = FALSE)
    }
  }
  if (anyNA(data$threshold)) {
    stop("Missing thresholds in calibration plot data: ", file, call. = FALSE)
  }
  data
})

# --- Check plotting dependencies before rendering ---
required_packages <- c("ggplot2", "ggpubr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
})

# --- Create panels in row-major order and mark selected thresholds ---
# Callback: Build one panel and add its matching selected-threshold marker.
plots <- lapply(seq_len(nrow(calibration_plot_specs)), function(index) {
  spec <- calibration_plot_specs[index, ]
  value <- selected[[spec$parameter]]
  plot <- calibration_plot(curves[[index]], spec$label, spec$trans) +
    geom_vline(xintercept = value, colour = "red", linewidth = 0.5, show.legend = FALSE)
  # Zero is -Inf on a logarithmic axis, placing its line on the left boundary.
  if (spec$trans == "log10" && value == 0) {
    plot <- plot + labs(caption = "Selected E-value = 0 (left boundary of log axis)")
  }
  plot
})

# --- Arrange the eight panels and save the proportional PDF ---
composite <- ggarrange(plotlist = plots, nrow = 4, ncol = 2, common.legend = TRUE, legend = "bottom")
ggsave(args[[3]], composite, width = calibration_plot_width, height = 4 * calibration_plot_row_height, units = "in")
message("Selected calibration thresholds figure written: ", args[[3]])
