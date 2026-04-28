#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: validation_bins.R <validate_top_ids.rds> <bins> <output_dir>", call. = FALSE)
}

top_ids_file <- args[[1]]
bins <- suppressWarnings(as.integer(args[[2]]))
output_dir <- args[[3]]

if (is.na(bins) || bins < 1 || bins > 100) {
  stop("Bins must be an integer between 1 and 100.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

top_id <- readRDS(top_ids_file)

required_columns <- c("query", "support_identification_pct", "id_correct_close")
missing_columns <- setdiff(required_columns, names(top_id))
if (length(missing_columns) > 0) {
  stop(paste("Top IDs file is missing required columns:", paste(missing_columns, collapse = ", ")), call. = FALSE)
}

bins_df <- data.frame(
  range_support_main_id = cut(
    seq(0, 100, 100 / bins)[-1],
    breaks = seq(0, 100, 100 / bins),
    include.lowest = TRUE
  )
)

confidence_support <- top_id %>%
  mutate(range_support_main_id = cut(support_identification_pct, breaks = seq(0, 100, 100 / bins), include.lowest = TRUE)) %>%
  group_by(range_support_main_id) %>%
  reframe(
    count_all = n_distinct(query),
    count_correct = sum(id_correct_close == "correct"),
    count_close = sum(id_correct_close == "close"),
    count_wrong = sum(id_correct_close == "wrong"),
    percentage_correct = round(count_correct / count_all * 100, 2),
    percentage_close = round(count_close / count_all * 100, 2),
    percentage_wrong = round(count_wrong / count_all * 100, 2)
  )

confidence_support <- left_join(bins_df, confidence_support, by = "range_support_main_id") %>%
  rename(
    range_support = range_support_main_id,
    probability_correct = percentage_correct,
    probability_close = percentage_close,
    probability_wrong = percentage_wrong
  )

confidence_support_long <- confidence_support %>%
  pivot_longer(
    cols = c("probability_correct", "probability_close", "probability_wrong"),
    names_to = "identification",
    values_to = "percentage"
  )
confidence_support_long$identification <- factor(
  confidence_support_long$identification,
  c("probability_wrong", "probability_close", "probability_correct"),
  labels = c("wrong", "close", "correct")
)
confidence_support_long$percentage[confidence_support_long$percentage == 0] <- NA

confidence_support_plot <- ggplot(confidence_support_long, aes(x = range_support, y = percentage, fill = identification)) +
  geom_bar(position = "stack", stat = "identity", na.rm = TRUE) +
  geom_text(aes(range_support, 100, label = paste0("n=", count_all)), vjust = -0.5, size = 3, na.rm = TRUE) +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), position = position_stack(vjust = 0.5), size = 3, na.rm = TRUE) +
  scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 110)) +
  scale_fill_manual(values = c("correct" = "#67B891", "close" = "#D9A55A", "wrong" = "#C85A5A")) +
  labs(x = "Support for top identification (%)", y = "Samples (%)", fill = "Identification", title = paste(bins, "bins")) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13), axis.text.x = element_text(size = 12))

output_file <- file.path(output_dir, "validation_confidence_support.csv")
output_pdf <- file.path(output_dir, "validation_confidence_support.pdf")
output_connection <- file(output_file, "wb")
write.table(
  confidence_support[, c("range_support", "probability_correct", "probability_close", "probability_wrong")],
  file = output_connection,
  sep = ",",
  col.names = TRUE,
  row.names = FALSE
)
close(output_connection)
ggsave(output_pdf, confidence_support_plot, width = 8, height = 4)

message("Validation confidence support file written:")
message(output_file)
message("Validation confidence support plot written:")
message(output_pdf)
