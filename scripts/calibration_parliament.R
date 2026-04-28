#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 6) {
  stop("Usage: calibration_parliament.R <prepared_input.rds> <alignment_thresholds.tsv> <gene_threshold.tsv> <output_dir> <test_output_dir> <threshold_template.csv>", call. = FALSE)
}

prepared_input <- args[[1]]
alignment_file <- args[[2]]
gene_file <- args[[3]]
output_dir <- args[[4]]
test_output_dir <- args[[5]]
template_file <- args[[6]]

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(test_output_dir, recursive = TRUE, showWarnings = FALSE)
manual_input_dir <- file.path(output_dir, "manual_input_needed")
dir.create(manual_input_dir, recursive = TRUE, showWarnings = FALSE)

parliament_size_csv <- file.path(test_output_dir, "calibration_parliament_size.csv")
parliament_size_pdf <- file.path(test_output_dir, "calibration_parliament_size.pdf")
parliament_thresholds_tsv <- file.path(manual_input_dir, "calibration_parliament.tsv")

validate_alignment_thresholds <- function(file) {
  message("Step 1/7: Loading and checking alignment filtering thresholds.")
  expected <- c("min_similarity", "min_length", "max_gapopens", "max_mismatches", "max_evalue", "min_bitscore")
  thresholds <- read.delim(file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE, blank.lines.skip = TRUE)
  names(thresholds) <- trimws(names(thresholds))

  if (!identical(names(thresholds), c("parameter", "value"))) {
    stop("Alignment thresholds file must use tab-separated columns: parameter and value.", call. = FALSE)
  }

  thresholds$parameter <- trimws(thresholds$parameter)
  thresholds$value <- trimws(as.character(thresholds$value))

  if (any(duplicated(thresholds$parameter))) {
    stop("Alignment thresholds file contains duplicated parameter names.", call. = FALSE)
  }

  unexpected <- setdiff(thresholds$parameter, expected)
  if (length(unexpected) > 0) {
    stop("Alignment thresholds file contains unexpected parameter names: ", paste(unexpected, collapse = ", "), call. = FALSE)
  }

  missing <- setdiff(expected, thresholds$parameter)
  if (length(missing) > 0) {
    stop("Alignment thresholds file is missing required parameter names: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  thresholds <- thresholds[match(expected, thresholds$parameter), , drop = FALSE]

  if (any(is.na(thresholds$value)) || any(thresholds$value == "NA") || any(thresholds$value == "")) {
    stop("All values need to be specified before this file can be used.", call. = FALSE)
  }

  values <- setNames(as.numeric(thresholds$value), thresholds$parameter)

  ranges <- list(
    min_similarity = c(0, 100),
    min_length = c(0, 99999),
    max_gapopens = c(0, 99999),
    max_mismatches = c(0, 99999),
    max_evalue = c(0, 100),
    min_bitscore = c(0, 99999)
  )

  for (column in names(ranges)) {
    value <- values[[column]]
    if (is.na(value) || value < ranges[[column]][1] || value > ranges[[column]][2]) {
      stop(column, " must be between ", ranges[[column]][1], " and ", ranges[[column]][2], ".", call. = FALSE)
    }
  }

  message("Minimum alignment similarity: ", values[["min_similarity"]])
  message("Minimum alignment length: ", values[["min_length"]])
  message("Maximum number of gap openings: ", values[["max_gapopens"]])
  message("Maximum number of alignment mismatches: ", values[["max_mismatches"]])
  message("Maximum E-value: ", values[["max_evalue"]])
  message("Minimum Bit-score: ", values[["min_bitscore"]])

  as.list(values)
}

validate_gene_threshold <- function(file) {
  message("Step 2/7: Loading and checking gene performance threshold.")
  thresholds <- read.delim(file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE, blank.lines.skip = TRUE)
  names(thresholds) <- trimws(names(thresholds))

  if (!identical(names(thresholds), c("parameter", "value"))) {
    stop("Gene threshold file must use tab-separated columns: parameter and value.", call. = FALSE)
  }

  thresholds$parameter <- trimws(thresholds$parameter)
  thresholds$value <- trimws(as.character(thresholds$value))

  if (nrow(thresholds) != 1 || thresholds$parameter[1] != "min_gene_performance") {
    stop("Gene threshold file must contain exactly one row of threshold values.", call. = FALSE)
  }

  if (any(is.na(thresholds$value)) || any(thresholds$value == "NA") || any(thresholds$value == "")) {
    stop("All values need to be specified before this file can be used.", call. = FALSE)
  }

  value <- as.numeric(thresholds$value[1])
  if (is.na(value) || value < 0 || value > 100) {
    stop("min_gene_performance must be between 0 and 100.", call. = FALSE)
  }

  message("Minimum gene performance: ", value)
  value
}

thresholds <- validate_alignment_thresholds(alignment_file)
gene_performance_threshold <- validate_gene_threshold(gene_file)
message("Step 3/7: Loading prepared calibration data.")
ids <- readRDS(prepared_input)

message("Step 4/7: Applying alignment filters.")
filtered_ids <- ids %>%
  filter(
    pident >= thresholds$min_similarity,
    length >= thresholds$min_length,
    gapopen <= thresholds$max_gapopens,
    mismatch <= thresholds$max_mismatches,
    evalue <= thresholds$max_evalue,
    bitscore >= thresholds$min_bitscore
  )

summarise_threshold <- function(data, threshold) {
  data %>%
    group_by(query, target_sp, id_correct_close, query_samples) %>%
    reframe(count = n()) %>%
    group_by(query) %>%
    mutate(parliament_size = sum(count)) %>%
    arrange(query, desc(count)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    filter(parliament_size >= threshold) %>%
    reframe(
      count_all = n(),
      count_correct = sum(id_correct_close == "correct"),
      count_close = sum(id_correct_close == "close"),
      count_wrong = sum(id_correct_close == "wrong"),
      proportion_correct = count_correct / count_all,
      proportion_close = count_close / count_all,
      proportion_wrong = count_wrong / count_all,
      retrievability_all = count_all / query_samples
    ) %>%
    distinct() %>%
    mutate(threshold = threshold)
}

message("Step 5/7: Applying gene performance filter.")
ids_with_gene_performance <- filtered_ids %>%
  group_by(gene) %>%
  mutate(
    gene_count_all_postfiltering = n(),
    gene_count_correct_postfiltering = sum(id_correct_close == "correct"),
    gene_proportion_correct_postfiltering = gene_count_correct_postfiltering / gene_count_all_postfiltering * 100
  ) %>%
  ungroup() %>%
  filter(gene_proportion_correct_postfiltering >= gene_performance_threshold)

message("Step 6/7: Assessing minimum parliament size thresholds.")
parliament_size_df <- bind_rows(lapply(seq(0, 1000, 1), function(i) {
  summarise_threshold(ids_with_gene_performance, i)
}))

parliament_size_plot <- ggplot(parliament_size_df) +
  geom_line(aes(x = threshold, y = proportion_correct * 100), linewidth = 0.5) +
  geom_line(aes(x = threshold, y = retrievability_all * 100), linewidth = 0.5, linetype = "dashed") +
  scale_x_continuous(breaks = seq(0, 1000, 50)) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100), name = "Accuracy (%)", sec.axis = dup_axis(name = "Retrievability (%)")) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(), legend.key.width = unit(3, "line")) +
  labs(x = "Min. parliament size (n genes)")

write.csv(parliament_size_df, parliament_size_csv, row.names = FALSE)
ggsave(parliament_size_pdf, parliament_size_plot, width = 8, height = 4)

message("Step 7/7: Writing editable parliament size threshold TSV.")
template <- read.csv(template_file, stringsAsFactors = FALSE, check.names = FALSE)
parliament_template <- data.frame(parameter = "min_parliament_size", value = NA, stringsAsFactors = FALSE)
parliament_template <- parliament_template[parliament_template$parameter %in% names(template), , drop = FALSE]
write.table(parliament_template, parliament_thresholds_tsv, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

message("Output files written:")
message(parliament_size_csv)
message(parliament_size_pdf)
message(parliament_thresholds_tsv)
