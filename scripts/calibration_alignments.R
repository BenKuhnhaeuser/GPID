#!/usr/bin/env Rscript
# --- Purpose ---
# Assess six alignment thresholds and write their diagnostic curves and editable threshold TSV.

# --- Read paths and initialize output locations ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: calibration_alignments.R <prepared_input.rds> <output_dir> <test_output_dir>", call. = FALSE)
}

prepared_input <- args[[1]]
output_dir <- args[[2]]
test_output_dir <- args[[3]]

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
})

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(test_output_dir, recursive = TRUE, showWarnings = FALSE)
manual_input_dir <- file.path(output_dir, "manual_input_needed")
dir.create(manual_input_dir, recursive = TRUE, showWarnings = FALSE)
message("Loading prepared calibration data: ", prepared_input)
ids <- readRDS(prepared_input)

alignment_similarity_csv <- file.path(test_output_dir, "calibration_alignments_similarity.csv")
alignment_similarity_pdf <- file.path(test_output_dir, "calibration_alignments_similarity.pdf")
alignment_length_csv <- file.path(test_output_dir, "calibration_alignments_length.csv")
alignment_length_pdf <- file.path(test_output_dir, "calibration_alignments_length.pdf")
alignment_gapopens_csv <- file.path(test_output_dir, "calibration_alignments_gapopens.csv")
alignment_gapopens_pdf <- file.path(test_output_dir, "calibration_alignments_gapopens.pdf")
alignment_mismatches_csv <- file.path(test_output_dir, "calibration_alignments_mismatches.csv")
alignment_mismatches_pdf <- file.path(test_output_dir, "calibration_alignments_mismatches.pdf")
alignment_evalue_csv <- file.path(test_output_dir, "calibration_alignments_evalue.csv")
alignment_evalue_pdf <- file.path(test_output_dir, "calibration_alignments_evalue.pdf")
alignment_bitscore_csv <- file.path(test_output_dir, "calibration_alignments_bitscore.csv")
alignment_bitscore_pdf <- file.path(test_output_dir, "calibration_alignments_bitscore.pdf")
alignment_composite_pdf <- file.path(test_output_dir, "calibration_alignments_thresholds_composite.pdf")
alignment_thresholds_tsv <- file.path(manual_input_dir, "calibration_alignments.tsv")

# --- Performance summaries and shared plot construction ---
# summarise_threshold(): Compute top-identification correctness and sample retrievability after applying a candidate threshold.
summarise_threshold <- function(data, threshold) {
  data %>%
    group_by(query, target_sp, id_correct_close, query_samples) %>%
    reframe(count = n()) %>%
    group_by(query) %>%
    arrange(query, desc(count)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
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

# alignment_threshold_table(): Evaluate each candidate alignment cutoff and combine its performance summary rows.
alignment_threshold_table <- function(limits, filter_fn) {
  # Callback: Evaluate the filter and performance summary for each candidate cutoff.
  bind_rows(lapply(limits, function(i) summarise_threshold(filter_fn(ids, i), i)))
}

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]])
source(file.path(dirname(script_file), "calibration_plots.R"))
source(file.path(dirname(script_file), "calibration_parameters.R"))
alignment_plot <- calibration_plot

# --- Step 1/8: Assessing minimum alignment similarity thresholds. ---
message("Step 1/8: Assessing minimum alignment similarity thresholds.")
# Callback: Keep alignments meeting each candidate minimum percentage similarity.
alignment_similarity_df <- alignment_threshold_table(seq(50, 100, 1), function(data, i) filter(data, pident >= i))
alignment_similarity_plot <- alignment_plot(alignment_similarity_df, "Min. alignment similarity (%)")
write.csv(alignment_similarity_df, alignment_similarity_csv, row.names = FALSE)
ggsave(alignment_similarity_pdf, alignment_similarity_plot, width = 8, height = 4)

# --- Step 2/8: Assessing minimum alignment length thresholds. ---
message("Step 2/8: Assessing minimum alignment length thresholds.")
# Callback: Keep alignments meeting each candidate minimum aligned length.
alignment_length_df <- alignment_threshold_table(seq(0, 5000, 100), function(data, i) filter(data, length >= i))
alignment_length_plot <- alignment_plot(alignment_length_df, "Min. alignment length")
write.csv(alignment_length_df, alignment_length_csv, row.names = FALSE)
ggsave(alignment_length_pdf, alignment_length_plot, width = 8, height = 4)

# --- Step 3/8: Assessing maximum alignment gap opening thresholds. ---
message("Step 3/8: Assessing maximum alignment gap opening thresholds.")
# Callback: Keep alignments at or below each candidate gap-opening limit.
alignment_gapopens_df <- alignment_threshold_table(seq(0, 100, 1), function(data, i) filter(data, gapopen <= i))
alignment_gapopens_plot <- alignment_plot(alignment_gapopens_df, "Max. alignment gap openings")
write.csv(alignment_gapopens_df, alignment_gapopens_csv, row.names = FALSE)
ggsave(alignment_gapopens_pdf, alignment_gapopens_plot, width = 8, height = 4)

# --- Step 4/8: Assessing maximum alignment mismatch thresholds. ---
message("Step 4/8: Assessing maximum alignment mismatch thresholds.")
# Callback: Keep alignments at or below each candidate mismatch limit.
alignment_mismatch_df <- alignment_threshold_table(seq(0, 100, 1), function(data, i) filter(data, mismatch <= i))
alignment_mismatches_plot <- alignment_plot(alignment_mismatch_df, "Max. alignment mismatches")
write.csv(alignment_mismatch_df, alignment_mismatches_csv, row.names = FALSE)
ggsave(alignment_mismatches_pdf, alignment_mismatches_plot, width = 8, height = 4)

# --- Step 5/8: Assessing maximum E-value thresholds. ---
message("Step 5/8: Assessing maximum E-value thresholds.")
# Callback: Keep alignments at or below each candidate E-value limit.
alignment_evalue_df <- alignment_threshold_table(10^(-seq(0, 200, 10)), function(data, i) filter(data, evalue <= i))
alignment_evalue_plot <- alignment_plot(alignment_evalue_df, "Max. E-value", trans = "log10")
write.csv(alignment_evalue_df, alignment_evalue_csv, row.names = FALSE)
ggsave(alignment_evalue_pdf, alignment_evalue_plot, width = 8, height = 4)

# --- Step 6/8: Assessing minimum Bit-score thresholds. ---
message("Step 6/8: Assessing minimum Bit-score thresholds.")
# Callback: Keep alignments meeting each candidate minimum Bit-score.
alignment_bitscore_df <- alignment_threshold_table(seq(0, 10000, 100), function(data, i) filter(data, bitscore >= i))
alignment_bitscore_plot <- alignment_plot(alignment_bitscore_df, "Min. Bit-score")
write.csv(alignment_bitscore_df, alignment_bitscore_csv, row.names = FALSE)
ggsave(alignment_bitscore_pdf, alignment_bitscore_plot, width = 8, height = 4)

# --- Step 7/8: Combining alignment threshold plots. ---
message("Step 7/8: Combining alignment threshold plots.")
composite_plot <- ggarrange(
  alignment_similarity_plot, alignment_length_plot,
  alignment_gapopens_plot, alignment_mismatches_plot,
  alignment_evalue_plot, alignment_bitscore_plot,
  nrow = 3, ncol = 2, common.legend = TRUE, legend = "bottom"
)
ggsave(alignment_composite_pdf, composite_plot, width = calibration_plot_width, height = 3 * calibration_plot_row_height)

# --- Step 8/8: Writing editable alignment threshold TSV. ---
message("Step 8/8: Writing editable alignment threshold TSV.")
alignment_parameters <- c("min_similarity", "min_length", "max_gapopens", "max_mismatches", "max_evalue", "min_bitscore")
alignment_template <- calibration_threshold_rows(alignment_parameters)
write.table(alignment_template, alignment_thresholds_tsv, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

message("Output files written:")
message(alignment_similarity_csv)
message(alignment_similarity_pdf)
message(alignment_length_csv)
message(alignment_length_pdf)
message(alignment_gapopens_csv)
message(alignment_gapopens_pdf)
message(alignment_mismatches_csv)
message(alignment_mismatches_pdf)
message(alignment_evalue_csv)
message(alignment_evalue_pdf)
message(alignment_bitscore_csv)
message(alignment_bitscore_pdf)
message(alignment_composite_pdf)
message(alignment_thresholds_tsv)
