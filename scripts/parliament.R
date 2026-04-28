#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

abort <- function(...) {
  message_text <- sprintf(...)
  cat(sprintf("Error: %s\n", message_text), file = stderr())
  quit(save = "no", status = 1)
}

warn_user <- function(...) {
  message_text <- sprintf(...)
  cat(sprintf("Warning: %s\n", message_text), file = stderr())
}

required_packages <- c("dplyr", "ggplot2", "ggtext", "stringr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  abort(
    "Missing required R package(s): %s",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr, warn.conflicts = FALSE)
  library(ggplot2)
  library(ggtext)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (!(length(args) %in% c(5L, 6L, 7L, 8L, 9L))) {
  abort(
    "Usage: parliament.R <sample directory> <blast file> <gene performance file> <filtering thresholds file> <confidence support file> [species groups file] [output directory] [output formats] [overwrite output]"
  )
}

sample_dir <- args[[1]]
blast_file <- args[[2]]
gene_performance_file <- args[[3]]
filtering_thresholds_file <- args[[4]]
confidence_support_file <- args[[5]]
species_groups_file <- ""
output_dir <- "."
output_formats_raw <- "csv,pdf"
overwrite_output_raw <- "0"

if (length(args) >= 6L) {
  species_groups_file <- args[[6]]
}

if (length(args) >= 7L) {
  output_dir <- args[[7]]
}

if (length(args) >= 8L) {
  output_formats_raw <- args[[8]]
}

if (length(args) >= 9L) {
  overwrite_output_raw <- args[[9]]
}

sample_name <- basename(normalizePath(sample_dir, winslash = "/", mustWork = FALSE))

required_paths <- c(sample_dir, blast_file, gene_performance_file, filtering_thresholds_file, confidence_support_file)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  abort("Input file or directory not found: %s", paste(missing_paths, collapse = ", "))
}

if (nzchar(species_groups_file) && !file.exists(species_groups_file)) {
  abort("Species groups file not found: %s", species_groups_file)
}

normalised_output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)

if (!dir.exists(normalised_output_dir)) {
  created_output_dir <- tryCatch(
    dir.create(normalised_output_dir, recursive = TRUE, showWarnings = FALSE),
    warning = function(warning) FALSE,
    error = function(error) FALSE
  )

  if (!isTRUE(created_output_dir) && !dir.exists(normalised_output_dir)) {
    abort("Could not create output directory: %s", output_dir)
  }
}

parse_output_formats <- function(raw_formats) {
  tokens <- unlist(strsplit(raw_formats, ",", fixed = TRUE), use.names = FALSE)
  tokens <- tolower(trimws(tokens))
  tokens <- tokens[nzchar(tokens)]

  if (length(tokens) == 0) {
    abort("At least one output format must be supplied.")
  }

  allowed_formats <- c("csv", "jpg", "svg", "pdf", "all")
  invalid_formats <- setdiff(tokens, allowed_formats)

  if (length(invalid_formats) > 0) {
    abort(
      "Unsupported output format(s): %s. Allowed values are csv, jpg, svg, pdf and all.",
      paste(unique(invalid_formats), collapse = ", ")
    )
  }

  if ("all" %in% tokens) {
    return(c("csv", "jpg", "svg", "pdf"))
  }

  unique(tokens)
}

requested_output_formats <- parse_output_formats(output_formats_raw)

parse_overwrite_flag <- function(raw_value) {
  cleaned_value <- tolower(trimws(raw_value))

  if (cleaned_value %in% c("", "0", "false", "no")) {
    return(FALSE)
  }

  if (cleaned_value %in% c("1", "true", "yes")) {
    return(TRUE)
  }

  abort("Invalid overwrite-output flag value: %s", raw_value)
}

overwrite_output <- parse_overwrite_flag(overwrite_output_raw)

output_requested <- function(format_name) {
  format_name %in% requested_output_formats
}

output_path <- function(filename) {
  file.path(normalised_output_dir, filename)
}

planned_output_files <- function() {
  files <- character()

  if (output_requested("csv")) {
    files <- c(files, output_path(paste0(sample_name, "_gpid.csv")))
  }

  if (output_requested("jpg")) {
    files <- c(files, output_path(paste0(sample_name, "_gpid.jpg")))
  }

  if (output_requested("pdf")) {
    files <- c(files, output_path(paste0(sample_name, "_gpid.pdf")))
  }

  if (output_requested("svg")) {
    files <- c(files, output_path(paste0(sample_name, "_gpid.svg")))
  }

  files
}

check_existing_output_files <- function() {
  existing_files <- planned_output_files()
  existing_files <- existing_files[file.exists(existing_files)]

  if (length(existing_files) == 0) {
    return(invisible(NULL))
  }

  if (isTRUE(overwrite_output)) {
    warn_user(
      "Existing parliament output files will be overwritten in %s: %s",
      normalised_output_dir,
      paste(existing_files, collapse = ", ")
    )
    return(invisible(NULL))
  }

  abort(
    "Output files already exist in %s: %s. Re-run with output-directory overwrite enabled to replace them.",
    normalised_output_dir,
    paste(existing_files, collapse = ", ")
  )
}

check_existing_output_files()

read_csv_checked <- function(path, sep = ",", ...) {
  tryCatch(
    read.csv(
      path,
      sep = sep,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      strip.white = TRUE,
      ...
    ),
    error = function(error) {
      abort("Could not read %s: %s", path, conditionMessage(error))
    }
  )
}

ensure_columns <- function(data, required_columns, file_label) {
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    abort(
      "%s is missing required column(s): %s",
      file_label,
      paste(missing_columns, collapse = ", ")
    )
  }
}

ensure_numeric_columns <- function(data, columns, file_label) {
  for (column in columns) {
    converted <- suppressWarnings(as.numeric(data[[column]]))
    invalid <- which(is.na(converted) & !is.na(data[[column]]) & trimws(as.character(data[[column]])) != "")

    if (length(invalid) > 0) {
      abort(
        "%s contains non-numeric values in column '%s'.",
        file_label,
        column
      )
    }

    data[[column]] <- converted
  }

  data
}

extract_species_name <- function(values) {
  match <- str_match(values, "^([^_]+_[^_]+)")
  match[, 2]
}

extract_genus_name <- function(values) {
  match <- str_match(values, "^([^_]+)")
  match[, 2]
}

parse_range_support <- function(labels) {
  parsed <- lapply(labels, function(label) {
    match <- regexec("^([\\[(])\\s*([0-9]+(?:\\.[0-9]+)?)\\s*,\\s*([0-9]+(?:\\.[0-9]+)?)\\s*([\\]\\)])$", label, perl = TRUE)
    groups <- regmatches(label, match)[[1]]

    if (length(groups) == 0) {
      abort("Confidence support file contains an invalid range label: %s", label)
    }

    lower <- as.numeric(groups[[3]])
    upper <- as.numeric(groups[[4]])

    if (is.na(lower) || is.na(upper) || lower > upper) {
      abort("Confidence support file contains an invalid support range: %s", label)
    }

    data.frame(
      range_support = label,
      lower_inclusive = identical(groups[[2]], "["),
      lower = lower,
      upper = upper,
      upper_inclusive = identical(groups[[5]], "]"),
      stringsAsFactors = FALSE
    )
  })

  parsed <- bind_rows(parsed)

  if (nrow(parsed) == 0) {
    abort("Confidence support file does not contain any support ranges.")
  }

  ordered <- order(parsed$lower, parsed$upper)
  parsed <- parsed[ordered, , drop = FALSE]

  if (nrow(parsed) > 1) {
    for (index in seq_len(nrow(parsed) - 1)) {
      left <- parsed[index, , drop = FALSE]
      right <- parsed[index + 1, , drop = FALSE]

      overlap <- left$upper > right$lower ||
        (left$upper == right$lower && left$upper_inclusive && right$lower_inclusive)

      if (overlap) {
        abort(
          "Confidence support ranges overlap between %s and %s.",
          left$range_support,
          right$range_support
        )
      }
    }
  }

  parsed
}

find_matching_range <- function(value, parsed_ranges) {
  matches <- which(
    (value > parsed_ranges$lower | (parsed_ranges$lower_inclusive & value >= parsed_ranges$lower)) &
      (value < parsed_ranges$upper | (parsed_ranges$upper_inclusive & value <= parsed_ranges$upper))
  )

  if (length(matches) == 0) {
    return(NA_integer_)
  }

  matches[[1]]
}

write_minimal_output <- function(sample_name, status_code) {
  if (!output_requested("csv")) {
    return(invisible(NULL))
  }

  output <- data.frame(
    Sample = sample_name,
    Rank = NA,
    Identification = NA,
    Species_group = NA,
    Support_pct = NA,
    Support_count = NA,
    Parliament_size = 0,
    Data_checks = status_code,
    ID_correct_pct = NA,
    ID_close_pct = NA,
    ID_wrong_pct = NA
  )

  write.table(
    output,
    file = output_path(paste0(sample_name, "_gpid.csv")),
    sep = ",",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
}

safe_ggsave <- function(filename, plot_object, width = 8, height = 4) {
  tryCatch(
    ggsave(filename = filename, plot = plot_object, width = width, height = height),
    error = function(error) {
      warn_user("Could not save %s: %s", filename, conditionMessage(error))
    }
  )
}

blast_required_columns <- c("gene", "query", "target", "pident", "length", "mismatch", "gapopen", "evalue", "bitscore")
ids <- read_csv_checked(blast_file, sep = "\t")
ensure_columns(ids, blast_required_columns, "BLAST results file")

if (nrow(ids) == 0) {
  warn_user("BLAST results file is empty. No gene parliament can be calculated.")
  write_minimal_output(sample_name, "FAILED_2")
  quit(save = "no", status = 1)
}

ids <- ensure_numeric_columns(ids, c("pident", "length", "mismatch", "gapopen", "evalue", "bitscore"), "BLAST results file")

ids$gene <- as.character(ids$gene)
ids$query <- as.character(ids$query)
ids$target <- as.character(ids$target)

if (any(ids$gene == "" | is.na(ids$gene))) {
  abort("BLAST results file contains missing gene names.")
}

if (any(ids$query == "" | is.na(ids$query))) {
  abort("BLAST results file contains missing query names.")
}

if (any(ids$target == "" | is.na(ids$target))) {
  abort("BLAST results file contains missing target names.")
}

ids$identification <- extract_species_name(ids$target)

if (any(is.na(ids$identification))) {
  abort("BLAST target names must begin with <Genus>_<species> so the identification can be extracted.")
}

gene_performance <- read_csv_checked(gene_performance_file)
ensure_columns(gene_performance, c("gene", "performance"), "Gene performance file")
gene_performance <- gene_performance[, c("gene", "performance"), drop = FALSE]
gene_performance <- ensure_numeric_columns(gene_performance, "performance", "Gene performance file")
gene_performance$gene <- as.character(gene_performance$gene)

if (anyDuplicated(gene_performance$gene) > 0) {
  duplicated_genes <- unique(gene_performance$gene[duplicated(gene_performance$gene)])
  abort(
    "Gene performance file contains duplicated gene names: %s",
    paste(duplicated_genes, collapse = ", ")
  )
}

gene_performance <- rename(gene_performance, gene_performance_pct = performance)

filtering_thresholds <- read_csv_checked(filtering_thresholds_file)
threshold_columns <- c(
  "min_similarity",
  "min_length",
  "max_gapopens",
  "max_mismatches",
  "max_evalue",
  "min_bitscore",
  "min_gene_performance",
  "min_parliament_size"
)
ensure_columns(filtering_thresholds, threshold_columns, "Filtering thresholds file")

if (nrow(filtering_thresholds) < 1) {
  abort("Filtering thresholds file does not contain any threshold row.")
}

filtering_thresholds <- filtering_thresholds[1, threshold_columns, drop = FALSE]
filtering_thresholds <- ensure_numeric_columns(filtering_thresholds, threshold_columns, "Filtering thresholds file")

confidence_support <- read_csv_checked(confidence_support_file)
confidence_columns <- c("range_support", "probability_correct", "probability_close", "probability_wrong")
ensure_columns(confidence_support, confidence_columns, "Confidence support file")

if (nrow(confidence_support) < 1) {
  abort("Confidence support file does not contain any confidence bins.")
}

confidence_support <- confidence_support[, confidence_columns, drop = FALSE]
confidence_support$range_support <- as.character(confidence_support$range_support)

if (anyDuplicated(confidence_support$range_support) > 0) {
  duplicated_ranges <- unique(confidence_support$range_support[duplicated(confidence_support$range_support)])
  abort(
    "Confidence support file contains duplicated support ranges: %s",
    paste(duplicated_ranges, collapse = ", ")
  )
}

confidence_support <- ensure_numeric_columns(
  confidence_support,
  c("probability_correct", "probability_close", "probability_wrong"),
  "Confidence support file"
)

parsed_ranges <- parse_range_support(confidence_support$range_support)
confidence_support <- left_join(confidence_support, parsed_ranges, by = "range_support")

if (nzchar(species_groups_file)) {
  species_groups <- read_csv_checked(species_groups_file, colClasses = "character")
  ensure_columns(species_groups, c("genus_species", "species_group"), "Species groups file")
  species_groups <- species_groups[, c("genus_species", "species_group"), drop = FALSE]
} else {
  species_groups <- data.frame(
    genus_species = sort(unique(ids$identification)),
    species_group = extract_genus_name(sort(unique(ids$identification))),
    stringsAsFactors = FALSE
  )
}

species_groups$genus_species <- as.character(species_groups$genus_species)
species_groups$species_group <- as.character(species_groups$species_group)

if (anyDuplicated(species_groups$genus_species) > 0) {
  duplicated_species <- unique(species_groups$genus_species[duplicated(species_groups$genus_species)])
  abort(
    "Species groups file contains duplicated species names: %s",
    paste(duplicated_species, collapse = ", ")
  )
}

ids <- left_join(ids, gene_performance[, c("gene", "gene_performance_pct")], by = "gene")

missing_gene_performance <- sort(unique(ids$gene[is.na(ids$gene_performance_pct)]))
if (length(missing_gene_performance) > 0) {
  abort(
    "Gene performance values are missing for gene(s): %s",
    paste(missing_gene_performance, collapse = ", ")
  )
}

ids <- left_join(ids, species_groups, by = c("identification" = "genus_species"))

if (any(is.na(ids$species_group))) {
  missing_species <- sort(unique(ids$identification[is.na(ids$species_group)]))
  warn_user(
    "Species group assignments are missing for: %s. Falling back to 'Unassigned'.",
    paste(missing_species, collapse = ", ")
  )
  ids$species_group <- ifelse(is.na(ids$species_group), "Unassigned", ids$species_group)
}

groups <- unique(ids$species_group)
group_colors <- setNames(grDevices::rainbow(length(groups), s = 0.6, v = 0.8), groups)
ids$group_color <- unname(group_colors[ids$species_group])

query_names <- unique(ids$query)
if (length(query_names) > 1) {
  warn_user(
    "BLAST results contain more than one query name (%s). Output files will use the sample directory name '%s'.",
    paste(query_names, collapse = ", "),
    sample_name
  )
}

cat(sprintf("Genes before alignment filtering: %d\n", dplyr::n_distinct(ids$gene)))

filtered_ids <- ids %>%
  filter(
    pident >= filtering_thresholds$min_similarity[[1]],
    length >= filtering_thresholds$min_length[[1]],
    gapopen <= filtering_thresholds$max_gapopens[[1]],
    mismatch <= filtering_thresholds$max_mismatches[[1]],
    evalue <= filtering_thresholds$max_evalue[[1]],
    bitscore >= filtering_thresholds$min_bitscore[[1]]
  )

cat(sprintf("Genes after alignment filtering: %d\n", dplyr::n_distinct(filtered_ids$gene)))

if (nrow(filtered_ids) == 0) {
  warn_user("No gene passed alignment filtering thresholds.")
  write_minimal_output(sample_name, "FAILED_2")
  cat(
    paste(
      "Minimal gene parliament table written:",
      if (output_requested("csv")) output_path(paste0(sample_name, "_gpid.csv")) else "CSV output not requested.",
      "No gene parliament figure produced.",
      "Terminating.",
      sep = "\n"
    ),
    "\n"
  )
  quit(save = "no", status = 1)
}

parliament_size <- dplyr::n_distinct(filtered_ids$gene)
parliament_size_check <- ifelse(
  parliament_size >= filtering_thresholds$min_parliament_size[[1]],
  "PASSED",
  "FAILED_3"
)
parliament_size_value <- parliament_size
parliament_size_check_value <- parliament_size_check

gene_parliament <- filtered_ids %>%
  count(query, identification, species_group, group_color, name = "support_identification_count") %>%
  mutate(
    support_identification_pct = support_identification_count / parliament_size * 100,
    parliament_size = parliament_size,
    parliament_size_check = parliament_size_check
  )

set.seed(42)
gene_parliament <- gene_parliament %>%
  mutate(
    rank = rank(-support_identification_count, ties.method = "min"),
    random_order = sample.int(n())
  ) %>%
  arrange(desc(support_identification_pct), random_order)

if (identical(parliament_size_check, "FAILED_3")) {
  warn_user("Parliament size below minimum threshold. Treat results with caution.")
}

top_id <- gene_parliament %>% filter(rank == 1)

matched_confidence_index <- vapply(
  top_id$support_identification_pct,
  find_matching_range,
  integer(1),
  parsed_ranges = confidence_support[, c("lower_inclusive", "lower", "upper", "upper_inclusive")]
)

if (any(is.na(matched_confidence_index))) {
  unmatched_support <- top_id$support_identification_pct[is.na(matched_confidence_index)]
  warn_user(
    "No confidence support bin matched support percentage(s): %s",
    paste(round(unmatched_support, 3), collapse = ", ")
  )
}

top_id_confidence_support <- top_id %>%
  mutate(confidence_row = matched_confidence_index) %>%
  mutate(
    probability_correct = ifelse(
      is.na(confidence_row),
      NA_real_,
      confidence_support$probability_correct[confidence_row]
    ),
    probability_close = ifelse(
      is.na(confidence_row),
      NA_real_,
      confidence_support$probability_close[confidence_row]
    ),
    probability_wrong = ifelse(
      is.na(confidence_row),
      NA_real_,
      confidence_support$probability_wrong[confidence_row]
    )
  )

if (identical(parliament_size_check, "FAILED_3")) {
  top_id_confidence_support <- top_id_confidence_support %>%
    mutate(
      probability_correct = NA_real_,
      probability_close = NA_real_,
      probability_wrong = NA_real_
    )
}

top_id_confidence_support <- top_id_confidence_support %>%
  select(
    identification,
    parliament_size,
    parliament_size_check,
    probability_correct,
    probability_close,
    probability_wrong
  ) %>%
  distinct()

gene_parliament_support <- left_join(
  gene_parliament %>%
    select(
      query,
      rank,
      identification,
      species_group,
      support_identification_pct,
      support_identification_count
    ),
  top_id_confidence_support,
  by = "identification"
) %>%
  mutate(
    parliament_size = coalesce(parliament_size, parliament_size_value),
    parliament_size_check = coalesce(parliament_size_check, parliament_size_check_value)
  ) %>%
  relocate(
    Sample = query,
    Rank = rank,
    Identification = identification,
    Species_group = species_group,
    Support_pct = support_identification_pct,
    Support_count = support_identification_count,
    Parliament_size = parliament_size,
    Data_checks = parliament_size_check,
    ID_correct_pct = probability_correct,
    ID_close_pct = probability_close,
    ID_wrong_pct = probability_wrong
  )

if (output_requested("csv")) {
  write.table(
    gene_parliament_support,
    file = output_path(paste0(sample_name, "_gpid.csv")),
    sep = ",",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
}

gene_parliament_top10 <- gene_parliament %>%
  slice_head(n = 10) %>%
  arrange(support_identification_pct, desc(identification)) %>%
  mutate(
    identification_label = str_replace_all(identification, "_", " "),
    species_group_label = str_replace_all(species_group, "_", " "),
    label_html = paste0(
      "<i>",
      identification_label,
      "</i> (",
      species_group_label,
      "), ",
      support_identification_count,
      " genes (",
      round(support_identification_pct, 1),
      "%)"
    ),
    identification_label = factor(identification_label, levels = identification_label)
  )

max_chars <- max(nchar(gene_parliament_top10$label_html))

parliament_plot <- ggplot(gene_parliament_top10, aes(x = identification_label, y = support_identification_pct)) +
  geom_segment(aes(y = 0, yend = support_identification_pct)) +
  geom_richtext(
    aes(label = label_html),
    hjust = -0.05,
    size = 3,
    fill = NA,
    label.color = NA
  ) +
  geom_point(
    aes(fill = species_group),
    pch = 21,
    size = 3,
    colour = "black"
  ) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 100)) +
  coord_flip(clip = "off") +
  theme_classic() +
  labs(
    x = NULL,
    y = "Percentage of genes supporting identification",
    title = paste0("Identification of sample: ", sample_name),
    subtitle = paste0("Parliament size: ", parliament_size, " genes")
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(colour = "black", size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.margin = margin(5, max_chars * 4, 5, 5)
  )

if (output_requested("pdf")) {
  safe_ggsave(output_path(paste0(sample_name, "_gpid.pdf")), parliament_plot)
}

if (output_requested("jpg")) {
  safe_ggsave(output_path(paste0(sample_name, "_gpid.jpg")), parliament_plot)
}

if (output_requested("svg")) {
  safe_ggsave(output_path(paste0(sample_name, "_gpid.svg")), parliament_plot)
}
