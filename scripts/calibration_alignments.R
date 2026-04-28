#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Usage: calibration_alignments.R <prepared_input.rds> <output_dir> <test_output_dir> <threshold_template.csv>", call. = FALSE)
}

prepared_input <- args[[1]]
output_dir <- args[[2]]
test_output_dir <- args[[3]]
template_file <- args[[4]]

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

alignment_threshold_table <- function(limits, filter_fn) {
  bind_rows(lapply(limits, function(i) summarise_threshold(filter_fn(ids, i), i)))
}

alignment_plot <- function(df, x_label, breaks = waiver(), trans = "identity") {
  ggplot(df) +
    geom_line(aes(x = threshold, y = proportion_correct * 100), linewidth = 0.5) +
    geom_line(aes(x = threshold, y = retrievability_all * 100), linewidth = 0.5, linetype = "dashed") +
    scale_x_continuous(breaks = breaks, trans = trans) +
    scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100), name = "Accuracy (%)", sec.axis = dup_axis(name = "Retrievability (%)")) +
    theme_bw(base_size = 14) +
    theme(panel.grid = element_blank(), legend.key.width = unit(3, "line")) +
    labs(x = x_label)
}

message("Step 1/8: Assessing minimum alignment similarity thresholds.")
alignment_similarity_df <- alignment_threshold_table(seq(50, 100, 1), function(data, i) filter(data, pident >= i))
alignment_similarity_plot <- alignment_plot(alignment_similarity_df, "Min. alignment similarity (%)", seq(0, 100, 10))
write.csv(alignment_similarity_df, alignment_similarity_csv, row.names = FALSE)
ggsave(alignment_similarity_pdf, alignment_similarity_plot, width = 8, height = 4)

message("Step 2/8: Assessing minimum alignment length thresholds.")
alignment_length_df <- alignment_threshold_table(seq(0, 5000, 100), function(data, i) filter(data, length >= i))
alignment_length_plot <- alignment_plot(alignment_length_df, "Min. alignment length", seq(0, 10000, 1000))
write.csv(alignment_length_df, alignment_length_csv, row.names = FALSE)
ggsave(alignment_length_pdf, alignment_length_plot, width = 8, height = 4)

message("Step 3/8: Assessing maximum alignment gap opening thresholds.")
alignment_gapopens_df <- alignment_threshold_table(seq(0, 100, 1), function(data, i) filter(data, gapopen <= i))
alignment_gapopens_plot <- alignment_plot(alignment_gapopens_df, "Max. alignment gap openings", seq(0, 100, 20))
write.csv(alignment_gapopens_df, alignment_gapopens_csv, row.names = FALSE)
ggsave(alignment_gapopens_pdf, alignment_gapopens_plot, width = 8, height = 4)

message("Step 4/8: Assessing maximum alignment mismatch thresholds.")
alignment_mismatch_df <- alignment_threshold_table(seq(0, 100, 1), function(data, i) filter(data, mismatch <= i))
alignment_mismatches_plot <- alignment_plot(alignment_mismatch_df, "Max. alignment mismatches", seq(0, 100, 20))
write.csv(alignment_mismatch_df, alignment_mismatches_csv, row.names = FALSE)
ggsave(alignment_mismatches_pdf, alignment_mismatches_plot, width = 8, height = 4)

message("Step 5/8: Assessing maximum E-value thresholds.")
alignment_evalue_df <- alignment_threshold_table(10^(-seq(0, 200, 10)), function(data, i) filter(data, evalue <= i))
alignment_evalue_plot <- alignment_plot(alignment_evalue_df, "Max. E-value", 10^(-seq(0, 200, 50)), "log10")
write.csv(alignment_evalue_df, alignment_evalue_csv, row.names = FALSE)
ggsave(alignment_evalue_pdf, alignment_evalue_plot, width = 8, height = 4)

message("Step 6/8: Assessing minimum Bit-score thresholds.")
alignment_bitscore_df <- alignment_threshold_table(seq(0, 10000, 100), function(data, i) filter(data, bitscore >= i))
alignment_bitscore_plot <- alignment_plot(alignment_bitscore_df, "Min. Bit-score", seq(0, 10000, 2000))
write.csv(alignment_bitscore_df, alignment_bitscore_csv, row.names = FALSE)
ggsave(alignment_bitscore_pdf, alignment_bitscore_plot, width = 8, height = 4)

message("Step 7/8: Combining alignment threshold plots.")
composite_plot <- ggarrange(
  alignment_similarity_plot, alignment_length_plot,
  alignment_gapopens_plot, alignment_mismatches_plot,
  alignment_evalue_plot, alignment_bitscore_plot,
  nrow = 3, ncol = 2, common.legend = TRUE, legend = "bottom"
)
ggsave(alignment_composite_pdf, composite_plot, width = 9, height = 10)

message("Step 8/8: Writing editable alignment threshold TSV.")
template <- read.csv(template_file, stringsAsFactors = FALSE, check.names = FALSE)
alignment_parameters <- c("min_similarity", "min_length", "max_gapopens", "max_mismatches", "max_evalue", "min_bitscore")
alignment_template <- data.frame(
  parameter = alignment_parameters,
  value = NA,
  stringsAsFactors = FALSE
)
alignment_template <- alignment_template[alignment_template$parameter %in% names(template), , drop = FALSE]
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
