# --- Purpose ---
# Shared calibration plotting style and panel specifications used by individual and combined figures.
# Shared plotting definitions for the calibration steps and selected-threshold summary.
calibration_plot_width <- 9
calibration_plot_row_height <- 10 / 3

# calibration_plot(): Draw accuracy and retrievability on one percentage axis with automatic x ticks and a shared legend.
calibration_plot <- function(df, x_label, trans = "identity") {
  ggplot(df) +
    geom_line(aes(x = threshold, y = proportion_correct * 100, linetype = "Accuracy"), linewidth = 0.5) +
    geom_line(aes(x = threshold, y = retrievability_all * 100, linetype = "Retrievability"), linewidth = 0.5) +
    scale_x_continuous(trans = trans) +
    scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100), name = "Performance (%)") +
    scale_linetype_manual(name = NULL, values = c(Accuracy = "solid", Retrievability = "dashed")) +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom", legend.key.width = grid::unit(3, "line")) +
    labs(x = x_label)
}

# Row-major order matches the six-panel alignment figure, then genes/parliament.
calibration_plot_specs <- data.frame(
  parameter = c("min_similarity", "min_length", "max_gapopens", "max_mismatches",
                "max_evalue", "min_bitscore", "min_gene_performance", "min_parliament_size"),
  file = c("alignments/calibration_alignments_similarity.csv",
           "alignments/calibration_alignments_length.csv",
           "alignments/calibration_alignments_gapopens.csv",
           "alignments/calibration_alignments_mismatches.csv",
           "alignments/calibration_alignments_evalue.csv",
           "alignments/calibration_alignments_bitscore.csv",
           "genes/calibration_gene_performance_thresholds.csv",
           "parliament/calibration_parliament_size.csv"),
  label = c("Min. alignment similarity (%)", "Min. alignment length",
            "Max. alignment gap openings", "Max. alignment mismatches",
            "Max. E-value", "Min. Bit-score", "Min. gene performance (%)",
            "Min. parliament size (n genes)"),
  trans = c(rep("identity", 4), "log10", rep("identity", 3)),
  stringsAsFactors = FALSE
)
