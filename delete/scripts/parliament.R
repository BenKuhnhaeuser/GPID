#!/usr/bin/env Rscript

##---------##
## 0 Setup ##
##---------##

# 0.1 Load Libraries
#--------------------
# Data manipulation
library(dplyr, warn.conflicts = FALSE) # Data manipulation
library(tidyr) # Tidying data
library(stringr) # String manipulation
library(strex) # String manipulation
library(withr) # Set random seeds

# Plotting
library(ggplot2) # Plotting
library(ggtext) # Improved text rendering in ggplot


# 0.2 Specify input file names
#------------------------------

args <- commandArgs(trailingOnly=TRUE)

# Sample directory
sample_dir = args[1]

# BLAST results of matching samples against reference data set
blast_file = args[2]

# Gene performance calibration file
gene_performance_file <- args[3]

# Filtering thresholds calibration file
filtering_thresholds_file <- args[4]

# Confidence support calibration file
confidence_support_file <- args[5]

# Species group file that specifies for each species a species group (optional)
species_groups_file = args[6]


# 0.3 Read input files
#----------------------
# BLAST results
ids <- read.csv(blast_file, sep = "\t")

# Species groups
if (file.exists(species_groups_file)) { # Check whether species groups file exists
  species_groups <- read.csv(species_groups_file, colClasses = "character") # If it does exist, read the file
} else { # If species groups file doesn't exist, construct species groups file using genus names of target species
  species_groups <- ids %>% 
    mutate(genus_species = str_before_nth(target, "_", 2),
           species_group = str_before_nth(target, "_", 1)) %>%
    select(genus_species, species_group) %>%
    arrange(genus_species) %>%
    unique()
  } 

# Gene performance
gene_performance <- read.csv(gene_performance_file) %>% 
  select("gene", "performance") %>% # Select variables with gene name and percentage of correct identifications per gene
  rename("gene_performance_pct" = "performance") # Rename variable name of percentage of correct identifications to "gene_performance_pct"

# Filtering thresholds
filtering_thresholds <- read.csv(filtering_thresholds_file)

# Confidence support
confidence_support <- read.csv(confidence_support_file)


# 0.4 Data preparation
#-----------------------
# Get colours for species groups
groups <- unique(species_groups$species_group) # Get unique species group names
# Get palette of 280 different colours. These colours are all RColorBrewer palettes merged, starting with qualitative, then sequential, then diverging, and only keeping unique colours.
colours <- c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#FFED6F", "#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B", "#F7FCFD", "#E5F5F9", "#CCECE6", "#99D8C9", "#66C2A4", "#41AE76", "#238B45", "#006D2C", "#00441B", "#E0ECF4", "#BFD3E6", "#9EBCDA", "#8C96C6", "#8C6BB1", "#88419D", "#810F7C", "#4D004B", "#F7FCF0", "#E0F3DB", "#A8DDB5", "#7BCCC4", "#4EB3D3", "#2B8CBE", "#0868AC", "#084081", "#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#FFFFFF", "#F0F0F0", "#BDBDBD", "#969696", "#737373", "#525252", "#252525", "#000000", "#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603", "#7F2704", "#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000", "#FFF7FB", "#ECE7F2", "#D0D1E6", "#A6BDDB", "#74A9CF", "#3690C0", "#0570B0", "#045A8D", "#023858", "#ECE2F0", "#67A9CF", "#02818A", "#016C59", "#014636", "#F7F4F9", "#E7E1EF", "#D4B9DA", "#C994C7", "#DF65B0", "#CE1256", "#980043", "#67001F", "#FCFBFD", "#EFEDF5", "#DADAEB", "#BCBDDC", "#9E9AC8", "#807DBA", "#6A51A3", "#54278F", "#3F007D", "#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E", "#7A0177", "#49006A", "#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D", "#FFFFE5", "#F7FCB9", "#D9F0A3", "#ADDD8E", "#78C679", "#238443", "#006837", "#004529", "#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8", "#253494", "#081D58", "#FFF7BC", "#FEE391", "#FEC44F", "#FE9929", "#EC7014", "#CC4C02", "#993404", "#662506", "#FFEDA0", "#FED976", "#FEB24C", "#FC4E2A", "#BD0026", "#800026", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5", "#C7EAE5", "#80CDC1", "#35978F", "#01665E", "#C51B7D", "#DE77AE", "#F1B6DA", "#FDE0EF", "#F7F7F7", "#E6F5D0", "#B8E186", "#7FBC41", "#4D9221", "#762A83", "#9970AB", "#C2A5CF", "#E7D4E8", "#D9F0D3", "#A6DBA0", "#5AAE61", "#1B7837", "#B35806", "#E08214", "#FDB863", "#FEE0B6", "#D8DAEB", "#B2ABD2", "#8073AC", "#542788", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#E0E0E0", "#BABABA", "#878787", "#4D4D4D", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#FEE08B", "#D9EF8B", "#A6D96A", "#66BD63", "#1A9850", "#D53E4F", "#E6F598", "#ABDDA4", "#3288BD")
group_colors <- setNames(colours[1:length(groups)], groups) # Assign different colour to each species group
species_groups$group_color <- group_colors[species_groups$species_group] # Add group colour as a new column 

# Create new column with gene performance (percentage of correct identifications)
ids <- left_join(ids, gene_performance[,c("gene", "gene_performance_pct")], by = "gene")

# Create new columns with the species names of the reference samples
ids <- ids %>%
  mutate(identification = str_before_nth(target, "_", 2)) # Genus and species name, which are the first two words

# Create new columns containing the species group of query and reference species by matching against the species groups file (if provided)
ids <- left_join(ids, species_groups, by = c("identification" = "genus_species")) # Target

# Tidy up
factors <- c("gene", "query", "target", "identification", "species_group") # Define which variables should be factors
ids[factors] <- lapply(ids[factors], factor) # Save factor columns as factor


#-------------------------------
# 1.1 Calculate Gene Parliament
#-------------------------------
gene_parliament <- ids %>%

  ## Remove low confidence matches (Apply optimal alignment filtering thresholds)
  filter(pident >= filtering_thresholds$min_similarity &
           length >= filtering_thresholds$min_length &
           gapopen <= filtering_thresholds$max_gapopens &
           mismatch <= filtering_thresholds$max_mismatches &
           evalue <= filtering_thresholds$max_evalue &
           bitscore >= filtering_thresholds$min_bitscore) %>%

  ## Quantify matches per reference taxon (Calculate parliament size, i.e. number of genes retrieved per query sample following alignment filtering)
  mutate(parliament_size = n_distinct(gene)) %>%

  ## Test whether sample is data-deficient (i.e. whether sample has minimum parliament size)
  mutate(parliament_size_check = if_else(parliament_size >= filtering_thresholds$min_parliament_size, "PASSED", "FAILED_3")) %>%

  ## Calculate support for each species identification
  group_by(query, identification, species_group, group_color, parliament_size, parliament_size_check) %>%
  reframe(support_identification_count = n(), # Number of genes supporting ID
          support_identification_pct = support_identification_count / parliament_size * 100) %>% # Percentage of genes supporting ID
  distinct() %>%
  ungroup() %>%

  ## Assess ranking of results
  mutate(rank = rank(-support_identification_count, ties.method = "min")) %>%

  ## Arrange results by support for identification in descending order, then by a random value (if there should be multiple identifications with identical support)
  mutate(random_order = with_seed(42, sample(seq_len(n())))) %>% # Create a new column with random order to be able to break ties in arrange. The seed of the random number is arbitrarily set to 42 to allow reproducibility
  arrange(desc(support_identification_pct), random_order)


# Terminate if no genes passed filtering thresholds
#---------------------------------------------------

if (nrow(gene_parliament) == 0) {

  cat("WARNING: No gene passed alignment filtering thresholds.\n")

  gene_parliament_support <- data.frame(
    Sample = ids$query[1], Rank = NA, Identification = NA, Species_group = NA, Support_pct = NA, Support_count = NA, Parliament_size = 0, Data_checks = "FAILED_2", ID_correct_pct = "", ID_close_pct = "", ID_wrong_pct = "", stringsAsFactors = FALSE)

  write.table(gene_parliament_support, file = paste0(ids$query[1], "_gpid2.csv"), sep = ",", row.names = FALSE, quote = FALSE)

  cat("Minimal gene parliament table written:",
      paste0(ids$query[1], "_gpid.csv"),
      "No gene parliament figure produced.",
      "Terminating.",
      "\n",
      sep = "\n")

  quit(status = 1)
}

# Print warning if minimum parliament size check failed
#-------------------------------------------------------
if (gene_parliament$parliament_size_check[1] == "FAILED_3") {

  cat("WARNING: Parliament size below minimum threshold. Treat results with caution.\n")

}



# Retrieve top identification and assess confidence
#---------------------------------------------------

# Create confidence support table with NAs only if parliament size check failed
confidence_support_datacheck_failed <- confidence_support
confidence_support_datacheck_failed[,c(2:4)] <- NA

# Retrieve top identification
top_id <- gene_parliament %>%
  filter(rank == 1) # Select top identification (most genes)

# Get number of bins from confidence_support file
bins = nrow(confidence_support)

# Assign top identification to a bin
top_id_confidence <- top_id %>%
  mutate(range_support = cut(support_identification_pct, breaks = seq(0,100,100/bins), include.lowest = TRUE))

# Combine top identification with confidence information depending on whether parliament size check was passed
if (gene_parliament$parliament_size_check[1] == "PASSED") { # Check whether parliament size check was passed
  top_id_confidence_support <- left_join(top_id_confidence, confidence_support, by = "range_support") # If passed, link with confidence support table
} else {top_id_confidence_support <- left_join(top_id_confidence, confidence_support_datacheck_failed, by = "range_support")} # Otherwise, link with table that has NAs only

# Select only columns of interest
top_id_confidence_support <- top_id_confidence_support %>% select(identification, parliament_size, parliament_size_check, probability_correct, probability_close, probability_wrong)

# Add confidence estimate to gene parliament table
gene_parliament_support <- left_join(gene_parliament %>% select(!c(random_order, group_color, parliament_size, parliament_size_check)),
                                           top_id_confidence_support,
                                           by = "identification")

# Reorder table and rename variables
gene_parliament_support <- gene_parliament_support %>%
  relocate(Sample = query, Rank = rank, Identification = identification, Species_group = species_group, Support_pct = support_identification_pct, Support_count = support_identification_count, Parliament_size = parliament_size, Data_checks = parliament_size_check, ID_correct_pct = probability_correct, ID_close_pct = probability_close, ID_wrong_pct = probability_wrong
           )


# Save table
#------------
write.table(gene_parliament_support, file = paste0(gene_parliament$query[1], "_gpid.csv"), sep = ",", row.names = FALSE, quote = FALSE, na = "")


#--------------------------
# 1.2 Gene Parliament plot
#--------------------------

# Get top 10 identifications in gene parliament
gene_parliament_top10 <- gene_parliament %>%
  slice_head(n = 10) %>%
  arrange(support_identification_pct, desc(identification)) %>%
  mutate(identification = str_replace_all(identification, "_", " "), # Remove underscores in identification names
         species_group = str_to_sentence(species_group), # Change species group names to sentence case
         label_html = paste0("<i>", identification, "</i>", " (", species_group, "), ", support_identification_count, " genes (", round(support_identification_pct, 1), "%)"),
         identification=factor(identification, levels = identification))

# Estimate label width (to adjust plot margin dynamically)
max_chars <- max(nchar(gene_parliament_top10$label_html))

# Plot gene parliament
parliament_plot <- ggplot(gene_parliament_top10, aes(x = identification, y = support_identification_pct)) +
  geom_segment(aes(y = 0, yend = support_identification_pct)) +
  geom_richtext(aes(label = label_html),
                hjust = -0.05,
                size = 3,
                fill = NA, label.color = NA) +
  geom_point(pch = 21, size = 3, fill = gene_parliament_top10$group_color, colour = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 100)) +
  coord_flip(clip = "off") + # Flip coordinates to show identifications horizontally
  theme_classic() +
  labs(
    x = NULL,
    y = "Percentage of genes supporting identification",
    title = paste0("Identification of sample: ", gene_parliament_top10$query[1]),
    subtitle = paste0("Parliament size: ", gene_parliament_top10$parliament_size[1], " genes")) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(colour = "black", size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.margin = margin(5, max_chars * 4, 5, 5), # Plot margin dynamically depending on maximum label width
    legend.position = "none")

# Save
ggsave(filename = paste0(gene_parliament$query[1], "_gpid.pdf"), parliament_plot, width = 8, height = 4)
ggsave(filename = paste0(gene_parliament$query[1], "_gpid.jpg"), parliament_plot, width = 8, height = 4)
ggsave(filename = paste0(gene_parliament$query[1], "_gpid.svg"), parliament_plot, width = 8, height = 4)


