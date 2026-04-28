#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop("Usage: calibration_genes.R <prepared_input.rds> <alignment_thresholds.tsv> <output_dir> <test_output_dir> <threshold_template.csv>", call. = FALSE)
}

prepared_input <- args[[1]]
alignment_file <- args[[2]]
output_dir <- args[[3]]
test_output_dir <- args[[4]]
template_file <- args[[5]]

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(test_output_dir, recursive = TRUE, showWarnings = FALSE)
manual_input_dir <- file.path(output_dir, "manual_input_needed")
dir.create(manual_input_dir, recursive = TRUE, showWarnings = FALSE)

sumstats_per_gene_csv <- file.path(test_output_dir, "calibration_sumstats_per_gene.csv")
calibration_gene_performance_csv <- file.path(output_dir, "calibration_gene_performance.csv")
gene_performance_csv <- file.path(test_output_dir, "calibration_gene_performance_thresholds.csv")
gene_performance_pdf <- file.path(test_output_dir, "calibration_gene_performance_thresholds.pdf")
gene_thresholds_tsv <- file.path(manual_input_dir, "calibration_genes.tsv")

validate_alignment_thresholds <- function(file) {
  message("Step 1/6: Loading and checking alignment filtering thresholds.")
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
    stop("NAs in alignment thresholds file detected. All alignment thresholds need to be specified before this file can be used.", call. = FALSE)
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

thresholds <- validate_alignment_thresholds(alignment_file)
message("Step 2/6: Loading prepared calibration data.")
ids <- readRDS(prepared_input)

message("Step 3/6: Applying alignment filters.")
filtered_ids <- ids %>%
  filter(
    pident >= thresholds$min_similarity,
    length >= thresholds$min_length,
    gapopen <= thresholds$max_gapopens,
    mismatch <= thresholds$max_mismatches,
    evalue <= thresholds$max_evalue,
    bitscore >= thresholds$min_bitscore
  )

genelist_no_filtering <- ids %>%
  select(gene) %>%
  unique() %>%
  arrange(gene)

message("Step 4/6: Calculating per-gene identification performance.")
sumstats_genes_filtered <- filtered_ids %>%
  group_by(gene) %>%
  reframe(
    samples_retrieved_count = n(),
    id_correct_count = sum(id_correct_close == "correct"),
    id_close_count = sum(id_correct_close == "close"),
    id_wrong_count = sum(id_correct_close == "wrong"),
    id_correct_pct = id_correct_count / samples_retrieved_count * 100,
    id_close_pct = id_close_count / samples_retrieved_count * 100,
    id_wrong_pct = id_wrong_count / samples_retrieved_count * 100,
    retrievability_pct = samples_retrieved_count / query_samples * 100
  ) %>%
  select(gene, id_correct_pct, id_close_pct, id_wrong_pct, retrievability_pct) %>%
  unique() %>%
  arrange(gene)

sumstats_genes_filtered <- left_join(genelist_no_filtering, sumstats_genes_filtered, by = "gene")
write.table(sumstats_genes_filtered, sumstats_per_gene_csv, sep = ",", row.names = FALSE, quote = FALSE)

gene_performance <- sumstats_genes_filtered %>%
  select(gene, id_correct_pct) %>%
  rename(performance = id_correct_pct)

write.table(gene_performance, calibration_gene_performance_csv, sep = ",", row.names = FALSE, quote = FALSE)

message("Step 5/6: Assessing minimum gene performance thresholds.")
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

gene_performance_df <- bind_rows(lapply(seq(0, 100, 1), function(i) {
  filtered_ids %>%
    group_by(gene) %>%
    mutate(
      gene_count_all_postfiltering = n(),
      gene_count_correct_postfiltering = sum(id_correct_close == "correct"),
      gene_proportion_correct_postfiltering = gene_count_correct_postfiltering / gene_count_all_postfiltering * 100
    ) %>%
    ungroup() %>%
    filter(gene_proportion_correct_postfiltering >= i) %>%
    summarise_threshold(i)
}))

gene_performance_plot <- ggplot(gene_performance_df) +
  geom_line(aes(x = threshold, y = proportion_correct * 100), linewidth = 0.5) +
  geom_line(aes(x = threshold, y = retrievability_all * 100), linewidth = 0.5, linetype = "dashed") +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100), name = "Accuracy (%)", sec.axis = dup_axis(name = "Retrievability (%)")) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(), legend.key.width = unit(3, "line")) +
  labs(x = "Min. gene performance (%)")

write.csv(gene_performance_df, gene_performance_csv, row.names = FALSE)
ggsave(gene_performance_pdf, gene_performance_plot, width = 8, height = 4)

message("Step 6/6: Writing editable gene performance threshold TSV.")
template <- read.csv(template_file, stringsAsFactors = FALSE, check.names = FALSE)
gene_template <- data.frame(parameter = "min_gene_performance", value = NA, stringsAsFactors = FALSE)
gene_template <- gene_template[gene_template$parameter %in% names(template), , drop = FALSE]
write.table(gene_template, gene_thresholds_tsv, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

message("Output files written:")
message(sumstats_per_gene_csv)
message(calibration_gene_performance_csv)
message(gene_performance_csv)
message(gene_performance_pdf)
message(gene_thresholds_tsv)
