#!/usr/bin/env Rscript
# --- Purpose ---
# Filter known-sample matches, select top identifications, and compare confidence bin counts.

# --- Read CLI paths and load analysis dependencies ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Usage: confidence_estimate.R <prepared_input.rds> <gene_performance.csv> <thresholds_filtering.csv> <output_dir>", call. = FALSE)
}

prepared_input <- args[[1]]
gene_performance_file <- args[[2]]
filtering_thresholds_file <- args[[3]]
output_dir <- args[[4]]

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(withr)
})

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Input checks used before statistical calculations ---
# abort(): Stop processing with a formatted, user-readable error message.
abort <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

# ensure_columns(): Reject input tables missing required column names.
ensure_columns <- function(data, expected, label) {
  missing <- setdiff(expected, names(data))
  if (length(missing) > 0) {
    abort("%s is missing required column(s): %s", label, paste(missing, collapse = ", "))
  }
}

# ensure_numeric_range(): Convert required columns to numeric values and enforce their permitted ranges.
ensure_numeric_range <- function(data, columns, label, ranges) {
  for (column in columns) {
    if (any(is.na(data[[column]]))) {
      abort("%s contains missing values in column '%s'.", label, column)
    }

    value <- suppressWarnings(as.numeric(data[[column]]))
    if (any(is.na(value))) {
      abort("%s contains non-numeric values in column '%s'.", label, column)
    }

    data[[column]] <- value
    if (!is.null(ranges[[column]]) && any(value < ranges[[column]][[1]] | value > ranges[[column]][[2]])) {
      abort(
        "%s contains values outside the allowed range for '%s' (%s to %s).",
        label,
        column,
        ranges[[column]][[1]],
        ranges[[column]][[2]]
      )
    }
  }

  data
}

# read_csv_checked(): Read an input table and include the source filename in any parsing error.
read_csv_checked <- function(file, label) {
  data <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE)
  if (nrow(data) == 0) {
    abort("%s does not contain any data rows: %s", label, file)
  }
  data
}

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]])
source(file.path(dirname(script_file), "gene_performance.R"))
source(file.path(dirname(script_file), "confidence_helpers.R"))
# --- Read gene performance with shared NA handling ---
gene_performance <- read_gene_performance(gene_performance_file)

# --- Read and reshape the parameter/value threshold table ---
filtering_thresholds <- read_csv_checked(filtering_thresholds_file, "Filtering thresholds file")
expected_thresholds <- c(
  "min_similarity",
  "min_length",
  "max_gapopens",
  "max_mismatches",
  "max_evalue",
  "min_bitscore",
  "min_gene_performance",
  "min_parliament_size"
)
# Normalize the legacy eight-column, single-value-row format before shared checks.
if (ncol(filtering_thresholds) == 8 && !anyDuplicated(names(filtering_thresholds)) &&
    setequal(names(filtering_thresholds), expected_thresholds)) {
  if (nrow(filtering_thresholds) != 1) {
    abort("Legacy filtering thresholds file must contain exactly one value row.")
  }
  filtering_thresholds <- data.frame(parameter = names(filtering_thresholds),
                                    value = unlist(filtering_thresholds[1, ], use.names = FALSE),
                                    stringsAsFactors = FALSE)
}
if (!identical(names(filtering_thresholds), c("parameter", "value"))) {
  abort("Filtering thresholds file must use parameter,value columns or the eight unique parameter names followed by one value row.")
}
filtering_thresholds$parameter <- trimws(filtering_thresholds$parameter)
if (anyNA(filtering_thresholds$parameter) || anyDuplicated(filtering_thresholds$parameter) ||
    !setequal(filtering_thresholds$parameter, expected_thresholds)) {
  abort("Filtering thresholds file must contain each of the eight required parameters exactly once.")
}
filtering_thresholds <- as.data.frame(as.list(setNames(
  filtering_thresholds$value[match(expected_thresholds, filtering_thresholds$parameter)],
  expected_thresholds
)), check.names = FALSE)
filtering_thresholds <- ensure_numeric_range(
  filtering_thresholds,
  expected_thresholds,
  "Filtering thresholds file",
  list(
    min_similarity = c(0, 100),
    min_length = c(0, 99999),
    max_gapopens = c(0, 99999),
    max_mismatches = c(0, 99999),
    max_evalue = c(0, 100),
    min_bitscore = c(0, 99999),
    min_gene_performance = c(0, 100),
    min_parliament_size = c(0, 99999)
  )
)

# --- Load prepared matches and apply calibration filters ---
ids <- readRDS(prepared_input)
ensure_columns(
  ids,
  c("gene", "query", "target_sp", "target_group", "id_correct_close", "pident", "length", "gapopen", "mismatch", "evalue", "bitscore"),
  "Prepared confidence data"
)

ids <- left_join(ids, gene_performance, by = "gene")

missing_gene_performance <- sort(unique(as.character(ids$gene[is.na(ids$performance)])))
if (length(missing_gene_performance) > 0) {
  abort("Gene performance values are missing for gene(s): %s", paste(missing_gene_performance, collapse = ", "))
}

filtered_ids <- ids %>%
  group_by(query) %>%
  mutate(genes_per_sample_count = n_distinct(gene)) %>%
  ungroup() %>%
  filter(
    pident >= filtering_thresholds$min_similarity[[1]],
    length >= filtering_thresholds$min_length[[1]],
    gapopen <= filtering_thresholds$max_gapopens[[1]],
    mismatch <= filtering_thresholds$max_mismatches[[1]],
    evalue <= filtering_thresholds$max_evalue[[1]],
    bitscore >= filtering_thresholds$min_bitscore[[1]],
    performance >= filtering_thresholds$min_gene_performance[[1]]
  )

if (nrow(filtered_ids) == 0) {
  abort("No confidence match passed the calibration filtering thresholds.")
}

# --- Select a top identification per sample with reproducible tie handling ---
top_id <- filtered_ids %>%
  mutate(probe_kit_genes_postfiltering = n_distinct(gene)) %>%
  group_by(query) %>%
  mutate(genecount_query_postfiltering = n_distinct(gene)) %>%
  filter(genecount_query_postfiltering >= filtering_thresholds$min_parliament_size[[1]]) %>%
  ungroup() %>%
  group_by(query, target_sp, target_group, id_correct_close, genecount_query_postfiltering, query_samples) %>%
  reframe(
    support_identification_count = n(),
    support_identification_pct = support_identification_count / genecount_query_postfiltering * 100,
    genes_retrieved_postfiltering_pct = genecount_query_postfiltering / probe_kit_genes_postfiltering * 100
  ) %>%
  distinct() %>%
  ungroup() %>%
  group_by(query) %>%
  mutate(random_order = with_seed(42, sample(seq_len(n())))) %>%
  arrange(query, desc(support_identification_pct), random_order) %>%
  slice_head(n = 1) %>%
  ungroup()

if (nrow(top_id) == 0) {
  abort("No confidence sample passed the minimum parliament size threshold.")
}

# --- Summarize and plot exclusive confidence categories ---
# confidence(): Summarize exclusive correct, close, and wrong classes into a stacked confidence-bin figure.
confidence <- function(top_id, bins) {
  bins_df <- data.frame(
    range_support_main_id = cut(
      seq(0, 100, 100 / bins)[-1],
      breaks = seq(0, 100, 100 / bins),
      include.lowest = TRUE
    )
  )

  confidence_table <- top_id %>%
    mutate(range_support_main_id = cut(support_identification_pct, breaks = seq(0, 100, 100 / bins), include.lowest = TRUE)) %>%
    group_by(range_support_main_id) %>%
    reframe(
      count_all = n_distinct(query),
      count_correct = sum(id_correct_close == "correct"),
      count_close = sum(id_correct_close == "close"),
      count_wrong = sum(id_correct_close == "wrong"),
      percentage_correct = count_correct / count_all * 100,
      percentage_close = count_close / count_all * 100,
      percentage_wrong = count_wrong / count_all * 100
    )

  confidence_table <- left_join(bins_df, confidence_table, by = "range_support_main_id")

  confidence_table_long <- pivot_longer(
    confidence_table,
    cols = c("percentage_correct", "percentage_close", "percentage_wrong"),
    names_to = "identification",
    values_to = "percentage"
  )
  confidence_table_long$identification <- factor(
    confidence_table_long$identification,
    c("percentage_wrong", "percentage_close", "percentage_correct"),
    labels = c("wrong", "close", "correct")
  )
  confidence_table_long$percentage[confidence_table_long$percentage == 0] <- NA

  ggplot(confidence_table_long, aes(x = range_support_main_id, y = percentage, fill = identification)) +
    geom_bar(position = "stack", stat = "identity", na.rm = TRUE) +
    geom_text(aes(range_support_main_id, 100, label = paste0("n=", count_all)), vjust = -0.5, size = 3, na.rm = TRUE) +
    geom_text(aes(label = paste0(round(percentage, 1), "%")), position = position_stack(vjust = 0.5), size = 3, na.rm = TRUE) +
    scale_x_discrete(labels = confidence_bin_labels) +
    scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 110)) +
    scale_fill_manual(values = c("correct" = "#67B891", "close" = "#D9A55A", "wrong" = "#C85A5A")) +
    labs(x = "Genes supporting the top identification", y = "Samples (%)", fill = "Identification", title = paste(bins, "bins")) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13), axis.text.x = element_text(size = 12, angle = 40, hjust = 1))
}

# --- Save top identifications and the candidate-bin figures ---
top_ids_file <- file.path(output_dir, "confidence_top_ids.rds")
confidence_pdf <- file.path(output_dir, "confidence_estimate.pdf")

saveRDS(top_id, top_ids_file)

pdf(file = confidence_pdf)
for (bins in seq(1, 10, 1)) {
  print(confidence(top_id, bins))
}
dev.off()

message("Top confidence identifications written:")
message(top_ids_file)
message("Confidence plot written:")
message(confidence_pdf)
