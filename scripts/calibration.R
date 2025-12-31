######################################################
## GeneParliamentID pipeline                        ##
## Method calibration                               ##
## Benedikt Kuhnhaeuser                             ##
## 3 December 2025                                  ##
######################################################

# This script analyses your test data set to identify the optimal parameters for the GeneParliamentID (GPID) pipeline
# For each calibration step, one or multiple figures are produced to enable you to make an informed decision on the parameters optimal for your dataset
# IMPORTANT: Your input is required at steps 1.7, 2.3, 3.2 and 4.3 to specify the optimal parameters. These steps are marked with (INPUT REQUIRED) 

# Calibration steps
#-------------------
# 1. Identify optimal thresholds for alignment filtering 
# 2. Identify optimal threshold for gene performance
# 3. Identify optimal threshold for parliament size
# 4. Estimate confidence in top identification depending on the relative support for the identification

# Outputs
#---------
# Based on the specified parameters, the following output files are produced:
## calibration_gene_performance.csv: Gene performance under optimal alignment filtering thresholds
## calibration_filtering_thresholds.csv: Optimal filtering thresholds
## calibration_confidence_support.csv: Confidence in top identification depending on percentage of genes supporting identification
# These files can be used as input for the main GPID pipeline script.


##----------##
## 0. Setup ##
##----------##

# 0.1 Load Libraries
#--------------------
## NOTE: If libraries are not installed yet, install them using the following syntax:
## install.packages("dplyr")

# Data manipulation
library(dplyr) # Data manipulation
library(tidyr) # Tidying data
library(strex) # String manipulation
library(withr) # Set random seeds

# Plotting
library(ggplot2) # Plotting
library(ggpubr) # Combing multiple plots


# 0.2 Set up working environment (INPUT REQUIRED)
#--------------------------------
# Specify working directory
getwd() # Get current working directory
list.files() # Check whether working directory contains the script Method_validation.R and input files
# setwd("Path/to/working/directory") # If not, set working directory to folder that contains these files, then double-check contents again using list.files()

# Create folder for outputs
dir.create("./analyses/") # Figures and tables created during method calibration
dir.create("./calibration_files/") # Calibration files with optimal parameters


# 0.3 Specify input file names (INPUT REQUIRED)
#-----------------------------------------------
# BLAST results of matching test samples against reference data set
blast_file = "blast.tsv" 

# Species groups file that specifies for each species a species group (optional)
species_groups_file = "species_groups.csv" 


# 0.4 Read input files
#----------------------
# BLAST results
ids <- read.csv(blast_file, sep = "\t") 

# Species groups (if file exists)
ifelse (file.exists(species_groups_file), # Check whether species groups file exists
        species_groups <- read.csv(species_groups_file, colClasses = "character"), # If it does exist, read the file
        species_groups <- data.frame(genus_species = character(), species_group = character())) # Otherwise, create an empty data frame

# Create empty data frame
species_groups <- data.frame(genus_species = character(), species_group = character())

# TO DO (optional): Specify species groups file that specifies for each species a species group. This overwrites the empty species groups data frame, if provided.
species_groups <- read.csv(species_groups_file, colClasses = "character") 


# 0.5 Data preparation
#-----------------------
# Create new columns with the species names of query and reference, and the number of query samples
ids <- ids %>%
  mutate(query_sp = str_before_nth(query, "_", 2), # Genus and species name for the test samples ("query"), which are the first two words
         target_sp = str_before_nth(target, "_", 2), # Genus and species name for reference ("target")
         query_samples = n_distinct(query)) # Number of query samples

# Create new columns containing the species group of query and reference species by matching against the species groups file (if provided)
ids <- left_join(ids, species_groups, by = c("query_sp" = "genus_species")) %>% rename(query_group = species_group) # Query
ids <- left_join(ids, species_groups, by = c("target_sp" = "genus_species")) %>% rename(target_group = species_group) # Target

# Check whether identification is correct for each query-reference pair
ids <- ids %>% mutate(id_correct = (query_sp == target_sp)) # Species
ids <- ids %>% mutate(id_correct_group = (query_group == target_group)) # Species group

# Combine "correct", "close" and "wrong" ids in a single column
ids$id_correct_close <- ifelse(ids$id_correct == TRUE, "correct",
                               ifelse(is.na(ids$id_correct_group), "wrong",
                                      ifelse(ids$id_correct_group == TRUE, "close", "wrong")))

# Tidy up
factors <- c("gene", "query", "target", "query_sp", "target_sp", "query_group", "target_group", "id_correct", "id_correct_group", "id_correct_close") # Define which variables should be factors
ids[factors] <- lapply(ids[factors], factor) # Save factor columns as factor
ids$id_correct_close <- factor(ids$id_correct_close, c("correct", "close", "wrong")) # Reorder factor levels of "id_correct_close"



##----------------------------------------------------##
## 1. Identify optimal alignment filtering thresholds ##
##----------------------------------------------------##

## This step explores the impacts of different filtering thresholds for the following alignment variables:

## The tested alignment variables are: 
### 1.1 Minimum alignment similarity
### 1.2 Minimum alignment length
### 1.3 Maximum alignment gap openings
### 1.4 Maximum alignment mismatches
### 1.5 Maximum E-value
### 1.6 Minimum Bit-score

## Each variable is assessed independently.

## NOTE - Adjusting thresholds
## Thresholds are always defined as a sequence (seq) of values between a minimum value and maximum value, with the specified step size
## E.g., "seq(50,100,1)" creates thresholds from 50 to 100 in steps of 1
## You can adjust the filtering thresholds for each variable if desired by changing the minimum, maximum or step size


# 1.1 Minimum alignment similarity
#----------------------------------

# Thresholds
pident_lims <- seq(50,100,1) # Alignment similarity thresholds from 50 to 100 in steps of 1

# Table
alignment_similarity_df <- data.frame() # Start from empty data frame

for (i in pident_lims) { # Loop through thresholds
  tmp <- ids %>% # Save results for individual threshold to temporary variable
    filter(pident >= i) %>% # Set filtering threshold
    group_by(query, target_sp, id_correct_close, query_samples) %>%
    reframe(count = n()) %>% # Count number of genes with respective species identification per specimen
    group_by(query) %>%
    arrange(query, desc(count)) %>%
    slice_head() %>% # Retrieve best match per sample
    ungroup() %>%
    reframe(count_all = n(), count_correct = sum(id_correct_close == "correct"), count_close = sum(id_correct_close == "close"), count_wrong = sum(id_correct_close == "wrong"), proportion_correct = count_correct / count_all, proportion_close = count_close / count_all, proportion_wrong = count_wrong / count_all, retrievability_all = count_all / query_samples) %>% # Calculate accuracy and retrievability
    distinct() %>% # Only keep unique rows, as reframe() produces duplicate rows
    mutate(threshold = i) # Add the filtering threshold to each row
  alignment_similarity_df <- rbind(alignment_similarity_df, tmp)} # Add results to data frame

# Figure
alignment_similarity_plot <- ggplot(alignment_similarity_df) +
  geom_line(aes(x = threshold, y = proportion_correct * 100), linewidth = 0.5) + # Accuracy (solid line)
  geom_line(aes(x = threshold, y = retrievability_all * 100), linewidth = 0.5, linetype = "dashed") + # Retrievability (dashed line)
  scale_x_continuous(breaks = seq(0,100,10)) + # Set x axis scale
  scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100), name = "Accuracy (%)", sec.axis = dup_axis(name = "Retrievability")) + 
  theme_bw(base_size = 14) + # Define main theme
  theme(panel.grid = element_blank(), # Remove grid lines
        legend.key.width = unit(3, "line")) + # Increase width of legend lines
  labs(x = "Min. alignment similarity (%)") # Define x axis label and legend header
alignment_similarity_plot

# Save
write.csv(alignment_similarity_df, file = "./analyses/1_alignment_similarity.csv") # Table
ggsave(alignment_similarity_plot, file = "./analyses/1_alignment_similarity.pdf", width = 8, height = 4) # Figure


# 1.2 Minimum alignment length
#------------------------------

# Thresholds
length_lims <- seq(0,5000,100) # Alignment length thresholds from 0 to 5000 in steps of 100

# Table
alignment_length_df <- data.frame() # Start from empty data frame

for (i in length_lims) { # Loop through thresholds
  tmp <- ids %>% # Save results for individual threshold to temporary variable
    filter(length >= i) %>% # Set filtering threshold
    group_by(query, target_sp, id_correct_close, query_samples) %>%
    reframe(count = n()) %>% # Count number of genes with respective species identification per specimen
    group_by(query) %>%
    arrange(query, desc(count)) %>%
    slice_head() %>% # Retrieve best match per sample
    ungroup() %>%
    reframe(count_all = n(), count_correct = sum(id_correct_close == "correct"), count_close = sum(id_correct_close == "close"), count_wrong = sum(id_correct_close == "wrong"), proportion_correct = count_correct / count_all, proportion_close = count_close / count_all, proportion_wrong = count_wrong / count_all, retrievability_all = count_all / query_samples) %>% # Calculate accuracy and retrievability
    distinct() %>% # Only keep unique rows, as reframe() produces duplicate rows
    mutate(threshold = i) # Add the filtering threshold to each row
  alignment_length_df <- rbind(alignment_length_df, tmp)} # Add results to data frame

# Figure
alignment_length_plot <- ggplot(alignment_length_df) +
  geom_line(aes(x = threshold, y = proportion_correct * 100), linewidth = 0.5) + # Accuracy (solid line)
  geom_line(aes(x = threshold, y = retrievability_all * 100), linewidth = 0.5, linetype = "dashed") + # Retrievability (dashed line)
  scale_x_continuous(breaks = seq(0,10000,1000)) + # Set x axis scale
  scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100), name = "Accuracy (%)", sec.axis = dup_axis(name = "Retrievability (%)")) +
  theme_bw(base_size = 14) + # Define main theme
  theme(panel.grid = element_blank(), # Remove grid lines
        legend.key.width = unit(3, "line")) + # Increase width of legend lines
  labs(x = "Min. alignment length") # Define x axis label and legend header
alignment_length_plot

# Save
write.csv(alignment_length_df, file = "./analyses/1_alignment_length.csv") # Table
ggsave(alignment_length_plot, file = "./analyses/1_alignment_length.pdf", width = 8, height = 4) # Figure


# 1.3 Maximum alignment gap openings
#------------------------------------

# Thresholds
gapopen_lims <- seq(0,100,1) # Alignmnent gap openings thresholds from 0 to 100 in steps of 1

# Table
alignment_gapopens_df <- data.frame() # Start from empty data frame

for (i in gapopen_lims) { # Loop through thresholds
  tmp <- ids %>% # Save results for individual threshold to temporary variable
    filter(gapopen <= i) %>% # Set filtering threshold
    group_by(query, target_sp, id_correct_close, query_samples) %>%
    reframe(count = n()) %>% # Count number of genes with respective species identification per specimen
    group_by(query) %>%
    arrange(query, desc(count)) %>%
    slice_head() %>% # Retrieve best match per sample
    ungroup() %>%
    reframe(count_all = n(), count_correct = sum(id_correct_close == "correct"), count_close = sum(id_correct_close == "close"), count_wrong = sum(id_correct_close == "wrong"), proportion_correct = count_correct / count_all, proportion_close = count_close / count_all, proportion_wrong = count_wrong / count_all, retrievability_all = count_all / query_samples) %>% # Calculate accuracy and retrievability
    distinct() %>% # Only keep unique rows, as reframe() produces duplicate rows
    mutate(threshold = i) # Add the filtering threshold to each row
  alignment_gapopens_df <- rbind(alignment_gapopens_df, tmp)} # Add results to data frame

# Figure
alignment_gapopens_plot <- ggplot(alignment_gapopens_df) +
  geom_line(aes(x = threshold, y = proportion_correct * 100), linewidth = 0.5) + # Accuracy (solid line)
  geom_line(aes(x = threshold, y = retrievability_all * 100), linewidth = 0.5, linetype = "dashed") + # Retrievability (dashed line)
  scale_x_continuous(breaks = seq(0,100,20)) + # Set x axis scale
  scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100), name = "Accuracy (%)", sec.axis = dup_axis(name = "Retrievability (%)")) +
  theme_bw(base_size = 14) + # Define main theme
  theme(panel.grid = element_blank(), # Remove grid lines
        legend.key.width = unit(3, "line")) + # Increase width of legend lines
  labs(x = "Max. alignment gap openings") # Define x axis label and legend header
alignment_gapopens_plot

# Save
write.csv(alignment_gapopens_df, file = "./analyses/1_alignment_gapopens.csv") # Table
ggsave(alignment_gapopens_plot, file = "./analyses/1_alignment_gapopens.pdf", width = 8, height = 4) # Figure


# 1.4 Maximum alignment mismatches
#----------------------------------

# Thresholds
mismatch_lims <- seq(0,100,1) # Alignment mismatches thresholds from 0 to 100 in steps of 1

# Table
alignment_mismatch_df <- data.frame() # Start from empty data frame

for (i in mismatch_lims) { # Loop through thresholds
  tmp <- ids %>% # Save results for individual threshold to temporary variable
    filter(mismatch <= i) %>% # Set filtering threshold
    group_by(query, target_sp, id_correct_close, query_samples) %>%
    reframe(count = n()) %>% # Count number of genes with respective species identification per specimen
    group_by(query) %>%
    arrange(query, desc(count)) %>%
    slice_head() %>% # Retrieve best match per sample
    ungroup() %>%
    reframe(count_all = n(), count_correct = sum(id_correct_close == "correct"), count_close = sum(id_correct_close == "close"), count_wrong = sum(id_correct_close == "wrong"), proportion_correct = count_correct / count_all, proportion_close = count_close / count_all, proportion_wrong = count_wrong / count_all, retrievability_all = count_all / query_samples) %>% # Calculate accuracy and retrievability
    distinct() %>% # Only keep unique rows, as reframe() produces duplicate rows
    mutate(threshold = i) # Add the filtering threshold to each row
  alignment_mismatch_df <- rbind(alignment_mismatch_df, tmp)} # Add results to data frame

# Figure
alignment_mismatches_plot <- ggplot(alignment_mismatch_df) +
  geom_line(aes(x = threshold, y = proportion_correct * 100), linewidth = 0.5) + # Accuracy (solid line)
  geom_line(aes(x = threshold, y = retrievability_all * 100), linewidth = 0.5, linetype = "dashed") + # Retrievability (dashed line)
  scale_x_continuous(breaks = seq(0,100,20)) + # Set x axis scale
  scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100), name = "Accuracy (%)", sec.axis = dup_axis(name = "Retrievability (%)")) +
  theme_bw(base_size = 14) + # Define main theme
  theme(panel.grid = element_blank(), # Remove grid lines
        legend.key.width = unit(3, "line")) + # Increase width of legend lines
  labs(x = "Max. alignment mismatches") # Define x axis label and legend header
alignment_mismatches_plot

# Save
write.csv(alignment_mismatch_df, file = "./analyses/1_alignment_mismatches.csv") # Table
ggsave(alignment_mismatches_plot, file = "./analyses/1_alignment_mismatches.pdf", width = 8, height = 4) # Figure


# 1.5 Maximum alignment E-value
#-------------------------------

# Thresholds
evalue_lims <- c(10^(-seq(0,200,10))) # Alignment e-value thresholds from 10^-200 to 10^0 in steps of 10^-10

# Table
alignment_evalue_df <- data.frame() # Start from empty data frame

for (i in evalue_lims) { # Loop through thresholds
  tmp <- ids %>% # Save results for individual threshold to temporary variable
    filter(evalue <= i) %>% # Set filtering threshold
    group_by(query, target_sp, id_correct_close, query_samples) %>%
    reframe(count = n()) %>% # Count number of genes with respective species identification per specimen
    group_by(query) %>%
    arrange(query, desc(count)) %>%
    slice_head() %>% # Retrieve best match per sample
    ungroup() %>%
    reframe(count_all = n(), count_correct = sum(id_correct_close == "correct"), count_close = sum(id_correct_close == "close"), count_wrong = sum(id_correct_close == "wrong"), proportion_correct = count_correct / count_all, proportion_close = count_close / count_all, proportion_wrong = count_wrong / count_all, retrievability_all = count_all / query_samples) %>% # Calculate accuracy and retrievability
    distinct() %>% # Only keep unique rows, as reframe() produces duplicate rows
    mutate(threshold = i) # Add the filtering threshold to each row
  alignment_evalue_df <- rbind(alignment_evalue_df, tmp)} # Add results to data frame

# Figure
alignment_evalue_plot <- ggplot(alignment_evalue_df) +
  geom_line(aes(x = threshold, y = proportion_correct * 100), linewidth = 0.5) + # Accuracy (solid line)
  geom_line(aes(x = threshold, y = retrievability_all * 100), linewidth = 0.5, linetype = "dashed") + # Retrievability (dashed line)
  scale_x_continuous(trans = "log10", breaks = 10^(-seq(0,200,50))) + # Set x axis scale
  scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100), name = "Accuracy (%)", sec.axis = dup_axis(name = "Retrievability (%)")) +
  theme_bw(base_size = 14) + # Define main theme
  theme(panel.grid = element_blank(), # Remove grid lines
        legend.key.width = unit(3, "line")) + # Increase width of legend lines
  labs(x = "Max. E-value") # Define x axis label and legend header
alignment_evalue_plot

# Save
write.csv(alignment_evalue_df, file = "./analyses/1_alignment_evalue.csv") # Table
ggsave(alignment_evalue_plot, file = "./analyses/1_alignment_evalue.pdf", width = 8, height = 4) # Figure


# 1.6 Minimum alignment Bit-score
#---------------------------------

# Thresholds
bitscore_lims <- seq(0, 10000, 100) # Alignment Bit-score thresholds from 0 to 10,000 in steps of 100

# Table
alignment_bitscore_df <- data.frame() # Start from empty data frame

for (i in bitscore_lims) { # Loop through thresholds
  tmp <- ids %>% # Save results for individual threshold to temporary variable
    filter(bitscore >= i) %>% # Set filtering threshold
    group_by(query, target_sp, id_correct_close, query_samples) %>%
    reframe(count = n()) %>% # Count number of genes with respective species identification per specimen
    group_by(query) %>%
    arrange(query, desc(count)) %>%
    slice_head() %>% # Retrieve best match per sample
    ungroup() %>%
    reframe(count_all = n(), count_correct = sum(id_correct_close == "correct"), count_close = sum(id_correct_close == "close"), count_wrong = sum(id_correct_close == "wrong"), proportion_correct = count_correct / count_all, proportion_close = count_close / count_all, proportion_wrong = count_wrong / count_all, retrievability_all = count_all / query_samples) %>% # Calculate accuracy and retrievability
    distinct() %>% # Only keep unique rows, as reframe() produces duplicate rows
    mutate(threshold = i) # Add the filtering threshold to each row
  alignment_bitscore_df <- rbind(alignment_bitscore_df, tmp)} # Add results to data frame

# Figure
alignment_bitscore_plot <- ggplot(alignment_bitscore_df) +
  geom_line(aes(x = threshold, y = proportion_correct * 100), linewidth = 0.5) + # Accuracy (solid line)
  geom_line(aes(x = threshold, y = retrievability_all * 100), linewidth = 0.5, linetype = "dashed") + # Retrievability (dashed line)
  scale_x_continuous(breaks = seq(0,10000,2000)) + # Set x axis scale
  scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100), name = "Accuracy (%)", sec.axis = dup_axis(name = "Retrievability (%)")) +
  theme_bw(base_size = 14) + # Define main theme
  theme(panel.grid = element_blank(), # Remove grid lines
        legend.key.width = unit(3, "line")) + # Increase width of legend lines
  labs(x = "Min. Bit-score") # Define x axis label and legend header
alignment_bitscore_plot

# Save
write.csv(alignment_bitscore_df, file = "./analyses/1_alignment_bitscore.csv") # Table
ggsave(alignment_bitscore_plot, file = "./analyses/1_alignment_bitscore.pdf", width = 8, height = 4) # Figure


# 1.7 Combine all alignment filtering plots into composite figure
#-----------------------------------------------------------------
ggarrange(alignment_similarity_plot, alignment_length_plot, 
          alignment_gapopens_plot, alignment_mismatches_plot, 
          alignment_evalue_plot, alignment_bitscore_plot,
          nrow = 3, ncol = 2, common.legend = TRUE, legend = "bottom")
ggsave("./analyses/1_alignment_thresholds_composite.pdf", width = 9, height = 10)


## 1.7 Specify optimal alignment filtering thresholds (INPUT REQUIRED)
#-----------------------------------------------------

## Based on the filtering tests in Steps 1.1 to 1.6, you need to specify the optimal alignment filtering thresholds that you have identified
## To decide on which alignment filtering thresholds are optimal, we recommend balancing accuracy (solid line) with retrievability (dashed line)
## Keep in mind that stricter filtering thresholds often lead to higher accuracy but lower retrievability
## The values specified below are default values established in the main manuscript across rattans and mahoganies

## TODO: Change below default values to the optimal thresholds for your group (default ranges explored are given in parentheses)

# Minimum alignment similarity (0 to 100)
similarity_threshold = 98

# Minimum alignment length (0 to 10000)
length_threshold = 100

# Maximum alignment gap openings (0 to 100)
gap_threshold = 1

# Maximum alignment mismatches (0 to 100)
mismatch_threshold = 5

# Maximum E-value (1e-200 to 1e+0)
evalue_threshold = 1e-60

# Minimum Bit-score (0 to 10000)
bitscore_threshold = 200



##-----------------------------------------------------------------------------##
## 2. Estimate gene performance and identify optimal gene filtering thresholds ##
##-----------------------------------------------------------------------------##

# 2.1 Estimate gene performance
#-------------------------------

# This step estimates gene performance (percentage of samples correctly identified to species) for each gene following application of the alignment filters specified in Step 1.7
# The output table also includes the percentage of close (if specified) and wrong identifications for each gene, and the percentage of samples retrieved for each gene

# Get list of all genes prior to filtering
genelist_no_filtering <- ids %>%
  select(gene) %>%
  unique() %>%
  arrange(gene)

# Performance with optimal alignment filtering thresholds 
sumstats_genes_filtered <- ids %>%
  ## Apply optimal alignment filtering thresholds specified above (under 1.7)
  filter(pident >= similarity_threshold & 
         length >= length_threshold & 
         gapopen <= gap_threshold & 
         mismatch <= mismatch_threshold & 
         evalue <= evalue_threshold & 
         bitscore >= bitscore_threshold) %>%
  
  ## Calculate stats
  group_by(gene) %>%
  reframe(samples_retrieved_count = n(), # Count number of samples for which the gene was retrieved
          id_correct_count = sum(id_correct_close == "correct"), # Count number of samples that were correctly identified for the gene
          id_close_count = sum(id_correct_close == "close"),
          id_wrong_count = sum(id_correct_close == "wrong"),
          id_correct_pct = id_correct_count / samples_retrieved_count * 100, # Calculate percentage of samples that were correctly identified for the gene
          id_close_pct = id_close_count / samples_retrieved_count * 100,
          id_wrong_pct = id_wrong_count / samples_retrieved_count * 100,
          retrievability_pct = samples_retrieved_count / query_samples * 100) %>% # Calculate percentage of samples that were retrieved for the gene
  select(gene, id_correct_pct, id_close_pct, id_wrong_pct, retrievability_pct) %>% # Select relevant columns
  unique() %>%
  arrange(gene)

# Left join all genes with filtered genes (NA means that gene didn't pass filtering threshold)
sumstats_genes_filtered <- left_join(genelist_no_filtering, sumstats_genes_filtered, by = "gene")

# Save table
write.table(sumstats_genes_filtered, "./analyses/2_sumstats_per_gene.csv", sep = ",", row.names = FALSE, quote = FALSE)


# 2.2 Identify optimal threshold for minimum gene performance
#-------------------------------------------------------------

# Thresholds
gene_correct_lims <- seq(0,100,1) # Gene performance (percentage of correct identifications) thresholds from 0 to 100 in steps of 1 (NOTE: can be adjusted as needed)

# Table
gene_performance_df <- data.frame() # Start from empty data frame

for (i in gene_correct_lims) { # Loop through thresholds
  tmp <- ids %>% # Save results for individual threshold to temporary variable
    ## Apply optimal alignment filtering thresholds specified above (under 1.7)
    filter(pident >= similarity_threshold & 
             length >= length_threshold & 
             gapopen <= gap_threshold & 
             mismatch <= mismatch_threshold & 
             evalue <= evalue_threshold & 
             bitscore >= bitscore_threshold) %>%
    group_by(gene) %>%
    mutate(gene_count_all_postfiltering = n(),
           gene_count_correct_postfiltering = sum(id_correct_close == "correct"),
           gene_proportion_correct_postfiltering = gene_count_correct_postfiltering / gene_count_all_postfiltering * 100) %>%
    ungroup() %>%
    filter(gene_proportion_correct_postfiltering >= i) %>% # Filter genes by performance (proportion of correct identification)
    group_by(query, target_sp, id_correct_close, query_samples) %>%
    reframe(count = n()) %>% # Count number of genes with respective species identification per specimen
    group_by(query) %>%
    arrange(query, desc(count)) %>%
    slice_head() %>% # Retrieve best match per sample
    ungroup() %>%
    reframe(count_all = n(), count_correct = sum(id_correct_close == "correct"), count_close = sum(id_correct_close == "close"), count_wrong = sum(id_correct_close == "wrong"), proportion_correct = count_correct / count_all, proportion_close = count_close / count_all, proportion_wrong = count_wrong / count_all, retrievability_all = count_all / query_samples) %>% # Calculate accuracy and retrievability
    distinct() %>% # Only keep unique rows, as reframe() produces duplicate rows
    mutate(threshold = i) # Add the filtering threshold to each row
  gene_performance_df <- rbind(gene_performance_df, tmp)} # Add results to data frame

# Figure
gene_performance_plot <- ggplot(gene_performance_df) +
  geom_line(aes(x = threshold, y = proportion_correct * 100), linewidth = 0.5) + # Accuracy (solid line)
  geom_line(aes(x = threshold, y = retrievability_all * 100), linewidth = 0.5, linetype = "dashed") + # Retrievability (dashed line)
  scale_x_continuous(breaks = seq(0,100,10)) + # Set x axis scale
  scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100), name = "Accuracy (%)", sec.axis = dup_axis(name = "Retrievability (%)")) +
  theme_bw(base_size = 14) + # Define main theme
  theme(panel.grid = element_blank(), # Remove grid lines
        legend.key.width = unit(3, "line")) + # Increase width of legend lines
  labs(x = "Min. gene performance (%)") # Define x axis label and legend header
gene_performance_plot

# Save
write.csv(gene_performance_df, file = "./analyses/2_gene_performance.csv") # Table
ggsave(gene_performance_plot, file = "./analyses/2_gene_performance.pdf", width = 8, height = 4) # Figure


# 2.3 Specify optimal gene performance threshold (INPUT REQUIRED)
#------------------------------------------------

## Based on the filtering tests in Step 2.2, you need to specify the optimal filtering thresholds that you have identified
## To decide on which alignment filtering thresholds are optimal, we recommend balancing accuracy (solid line) with retrievability (dashed line)
## Keep in mind that stricter filtering thresholds often lead to higher accuracy but lower retrievability
## The value specified below is the default value established in the main manuscript across rattans and mahoganies

## TODO: Change below default value to the optimal thresholds for your group (default range explored is given in parentheses)

# Minimum gene performance (0 to 100)
gene_performance_threshold = 30



##-------------------------------------------------------##
## 3. Identify optimal minimum parliament size threshold ##
##-------------------------------------------------------##

# 3.1 Identify optimal threshold for minimum parliament size
#------------------------------------------------------------

# Thresholds
parliament_size_lims <- seq(0,1000,1) # Parliament size thresholds from 0 to 1000 in steps of 1 (NOTE: can be adjusted as needed)

# Table
parliament_size_df <- data.frame() # Start from empty data frame

for (i in parliament_size_lims) { # Loop through thresholds
  tmp <- ids %>% # Save results for individual threshold to temporary variable
    ## Apply optimal alignment filtering thresholds specified above (under 1.7)
    filter(pident >= similarity_threshold & 
             length >= length_threshold & 
             gapopen <= gap_threshold & 
             mismatch <= mismatch_threshold & 
             evalue <= evalue_threshold & 
             bitscore >= bitscore_threshold) %>%
    group_by(gene) %>%
    mutate(gene_count_all_postfiltering = n(),
           gene_count_correct_postfiltering = sum(id_correct_close == "correct"),
           gene_proportion_correct_postfiltering = gene_count_correct_postfiltering / gene_count_all_postfiltering * 100) %>%
    ungroup() %>%
    ## Apply optimal gene performance threshold specified above (under 2.3)
    filter(gene_proportion_correct_postfiltering >= gene_performance_threshold) %>%
    group_by(query, target_sp, id_correct_close, query_samples) %>%
    reframe(count = n()) %>% # Count number of genes with respective species identification per specimen
    group_by(query) %>%
    mutate(parliament_size = sum(count)) %>% # Calculate parliament size (number of gene identifications per replicate following the above filtering)
    arrange(query, desc(count)) %>%
    slice_head() %>% # Retrieve best match per sample
    ungroup() %>%
    filter(parliament_size >= i) %>% # parliament_size_lims: Filter by number of genes in parliament following above filtering
    reframe(count_all = n(), count_correct = sum(id_correct_close == "correct"), count_close = sum(id_correct_close == "close"), count_wrong = sum(id_correct_close == "wrong"), proportion_correct = count_correct / count_all, proportion_close = count_close / count_all, proportion_wrong = count_wrong / count_all, retrievability_all = count_all / query_samples) %>% # Calculate accuracy and retrievability
    distinct() %>% # Only keep unique rows, as reframe() produces duplicate rows
    mutate(threshold = i) # Add the filtering threshold to each row
  parliament_size_df <- rbind(parliament_size_df, tmp)} # Add results to data frame

# Figure
parliament_size_plot <- ggplot(parliament_size_df) +
  geom_line(aes(x = threshold, y = proportion_correct * 100), linewidth = 0.5) + # Accuracy (solid line)
  geom_line(aes(x = threshold, y = retrievability_all * 100), linewidth = 0.5, linetype = "dashed") + # Retrievability (dashed line)
  scale_x_continuous(breaks = seq(0,1000,50)) + # Set x axis scale
  scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100), name = "Accuracy (%)", sec.axis = dup_axis(name = "Retrievability (%)")) +
  theme_bw(base_size = 14) + # Define main theme
  theme(panel.grid = element_blank(), # Remove grid lines
        legend.key.width = unit(3, "line")) + # Increase width of legend lines
  labs(x = "Min. parliament size (n genes)") # Define x axis label and legend header
parliament_size_plot

# Save
write.csv(parliament_size_df, file = "./analyses/3_parliament_size.csv") # Table
ggsave(parliament_size_plot, file = "./analyses/3_parliament_size.pdf", width = 8, height = 4) # Figure


# 3.2 Specify optimal minimum parliament size threshold (INPUT REQUIRED)
#-------------------------------------------------------

## Based on the filtering tests in Step 3.1, you need to specify the optimal filtering thresholds that you have identified
## To decide on which alignment filtering thresholds are optimal, we recommend balancing accuracy (solid line) with retrievability (dashed line)
## Keep in mind that stricter filtering thresholds often lead to higher accuracy but lower retrievability
## The value specified below is the default value established in the main manuscript across rattans and mahoganies

## TODO: Change below default value to the optimal thresholds for your group (default range explored is given in parentheses)

# Minimum parliament size (0 to 1000)
parliament_size_threshold = 10



##-------------------------------------------------------------------------------------------------------##
## 4. Estimate confidence in top identification depending on the relative support for the identification ##
##-------------------------------------------------------------------------------------------------------##

# 4.1 Retrieve top identification 
#---------------------------------

## This step retrieves the top identification for each test sample using the filtering thresholds specified above

top_id <- ids %>%

  ## Calculate percentage of genes retrieved per sample
  group_by(query) %>%
  mutate(genes_per_sample_count = n_distinct(gene)) %>% # Number of genes retrieved per sample
  ungroup() %>%
  
  ## Apply optimal alignment filtering thresholds specified above (under 1.7)
  filter(pident >= similarity_threshold & 
           length >= length_threshold & 
           gapopen <= gap_threshold & 
           mismatch <= mismatch_threshold & 
           evalue <= evalue_threshold & 
           bitscore >= bitscore_threshold) %>%
  
  ## Calculate gene accuracy following alignment filtering ##
  group_by(gene) %>%
  mutate(samples_per_gene_count_all_postfiltering = n_distinct(query), # Number of all query samples per gene remaining per gene following alignment filtering
         samples_per_gene_count_correct_postfiltering = sum(id_correct_close == "correct"), # Number of correctly identified query samples per gene following alignment filtering
         gene_proportion_correct_postfiltering = samples_per_gene_count_correct_postfiltering / samples_per_gene_count_all_postfiltering * 100) %>% # Gene performance (Percentage of correctly identified query samples following alignment filtering)
  ungroup() %>%
  
  ## Apply optimal gene performance threshold specified above (under 2.3)
  filter(gene_proportion_correct_postfiltering >= gene_performance_threshold) %>%
  
  ## Calculate number of genes per probe kit remaining following alignment filtering
  mutate(probe_kit_genes_postfiltering = n_distinct(gene)) %>% # Number of genes retrieved per probe kit following the above filtering
  # mutate(genecount = n_distinct(gene)) %>% # "genecount" is number of genes per probekit following the above filtering
  ungroup() %>%
  
  ## Calculate number of genes retrieved per query sample following alignment filtering
  group_by(query) %>%
  mutate(genecount_query_postfiltering = n_distinct(gene)) %>%
  
  ## Apply optimal parliament size filter specified above (under 3.2)
  filter(genecount_query_postfiltering >= parliament_size_threshold) %>% # Filter by number of genes in parliament following above filtering
  ungroup() %>%
  
  ## Calculate support for each species identification ##
  group_by(query, target_sp, target_group, id_correct_close, genecount_query_postfiltering, query_samples) %>%
  reframe(support_identification_count = n(), # Number of genes supporting ID
          support_identification_pct = support_identification_count / genecount_query_postfiltering * 100, # Percentage of genes supporting ID
          genes_retrieved_postfiltering_pct = genecount_query_postfiltering/ probe_kit_genes_postfiltering * 100) %>% # Percentage of genes retrieved for sample following filtering relative to number of genes in probe kit remaining following filtering
  distinct() %>%
  ungroup() %>%
  
  ## Arrange results by query, then by support for identification in descending order, then by a random value (if there should be multiple identifications with identical support)
  group_by(query) %>%
  mutate(random_order = with_seed(42, sample(seq_len(n())))) %>% # Create a new column with random order to be able to break ties in arrange. The seed of the random number is arbitrarily set to 42 to allow reproducibility
  arrange(query, desc(support_identification_pct), random_order) %>%
  slice_head() # Select top identification


# 4.2 Estimate confidence depending on support for different number of bins
#---------------------------------------------------------------------------

## The function "confidence(bins)" calculates the accuracy of identification depending on the percentage of genes supporting the top identification
## The range of support is first divided into n equal-sized bins
## By default, support is divided into between 1 and 10 bins (can be adjusted freely)
## The top identifications are assigned to bins depending on the percentage of genes supporting the identification
## For each bin, the percentage of correct, close and wrong identifications is estimated, and the number of samples assigned to the bin is calculated
## For example, if the number of bins is set to 5, five categories of 0-20%, >20-40%, >40-60%, >60-80% and >80-100% are created. A top identification supported by 37% of genes would be assigned the category ">20-40%".
## OUTPUTS: For each number of bins, a plot showing support for each bin is produced

# Define number of bins to explore (Adjust if needed)
bins_minimum <- 1 # Minimum number of bins
bins_maximum <- 10 # Maximum number of bins
bins_step_size <- 1 # Steps size between minimum and maximum

bins_to_plot <- seq(from = bins_minimum, to = bins_maximum, by = bins_step_size) # Combine values into a sequence

# Function confidence(bins)
confidence <- function(bins){

# Create data frame only holding bins
bins_df <- data_frame(range_support_main_id = 
                        cut(seq(0,100,100/bins)[-1], # Get the upper ends of the bins (exclude the first value, which is 0)
                            breaks = seq(0,100,100/bins), # Assign each of the values to its bin
                            include.lowest = TRUE)) # Include 0 in the first bin

# Assign top identifications to bins
confidence_table <- top_id %>% 
  mutate(range_support_main_id = cut(support_identification_pct, breaks = seq(0,100,100/bins), include.lowest = TRUE)) %>% # Add column that defines in which bin data are
  group_by(range_support_main_id) %>%
  reframe(count_all = n_distinct(query),
          count_correct = sum(id_correct_close == "correct"),
          count_close = sum(id_correct_close == "close"),
          count_wrong = sum(id_correct_close == "wrong"),
          percentage_correct = count_correct / count_all * 100,
          percentage_close = count_close / count_all  * 100,
          percentage_wrong = count_wrong / count_all * 100
          )

# Combine bins data frame with binned top identifications
confidence_table <- left_join(bins_df, confidence_table, by = "range_support_main_id")

# Prepare data for plotting
confidence_table_long <- pivot_longer(confidence_table, # Transform data frame to long format
                         cols = c("percentage_correct", "percentage_close", "percentage_wrong"),
                         names_to = "identification",
                         values_to = "percentage")
confidence_table_long$identification <- as.factor(confidence_table_long$identification) # Save as factor
confidence_table_long$identification <- factor(confidence_table_long$identification, c("percentage_correct", "percentage_close", "percentage_wrong")) # Reorder factors
levels(confidence_table_long$identification) <- c("correct", "close", "wrong") # Rename factors
confidence_table_long$percentage[confidence_table_long$percentage == 0] <- NA # Replace 0s with NA to avoid labeling of bars with 0%

# Plot
confidence_plot <- confidence_table_long %>%
  ggplot(aes(x = range_support_main_id, y = percentage, fill = forcats::fct_rev(identification))) +
  geom_bar(position = "stack", stat = "identity", na.rm = TRUE) +
  geom_text(aes(range_support_main_id, 100, label = paste0("n=", count_all)), vjust = -0.5, size = 3, na.rm = TRUE) +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), position = position_stack(vjust = 0.5), size = 3, na.rm = TRUE) +
  # scale_x_discrete(breaks = seq(0,100,100/bins)) +
  # scale_x_discrete(labels = c("0-20%", ">20-40%", ">40-60%", ">60-80%", ">80-100%")) +
  scale_y_continuous(breaks = seq(0,100,20), limits = c(0,110)) +
  scale_fill_manual(values = c("correct" = "#67B891", "close" = "#D9A55A", "wrong" = "#C85A5A")) +  
  labs(x = "Support for top identification (%)", y = "Samples (%)", fill = "Identification", title = str_c(bins, " bins")) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face="bold", size = 13),
        axis.text.x=element_text(size=12))
print(confidence_plot)

}


# Save plots with different bin numbers into single pdf file.

pdf(file = "./analyses/4_confidence_bins.pdf")

# Plot different bins
for (bins in bins_to_plot)(
  confidence(bins))

dev.off()


# 4.3 Identify appropriate number of bins for data set (INPUT REQUIRED)
#------------------------------------------------------

## Based on the results of step 4.2, you need to decide on the number of bins into which support should be divided to calculate support
## To make this decision, compare the plots produced for different numbers of bins
## Note that there is a trade-off between too few bins, where you lose resolution, and too many bins, which may be unrepresentative due to low number of samples

# TO DO: specify optimal number of bins for your data set (default = 5, can be adjusted freeely)
bins = 5 



##---------------------------------------------------------------------##
## 5. Save optimal parameters to calibration files (NO CHANGES NEEDED) ##
##---------------------------------------------------------------------##

## The below code saves all the optimal parameters that you have specified above into calibration files in the folder "calibration_files". 
## You can use these files directly as input for the main pipeline.

# 5.1 Gene list
#---------------

# Get list of genes above gene performance threshold
gene_performance <- sumstats_genes_filtered %>%
  filter(id_correct_pct >= gene_performance_threshold) %>%
  select(gene, id_correct_pct) %>%
  rename(performance = id_correct_pct)

# Save
output.file <- file("./calibration_files/calibration_gene_performance.csv", "wb") # Use this construct to save in Unix file format
write.table(gene_performance, file = output.file, sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
close(output.file)


# 5.2 Filtering thresholds
#--------------------------

# Compile optimal filtering thresholds
filtering_thresholds <- data.frame(
  "min_similarity" = similarity_threshold,
  "min_length" = length_threshold,
  "max_gapopens" = gap_threshold,
  "max_mismatches" = mismatch_threshold,
  "max_evalue" = evalue_threshold,
  "min_bitscore" = bitscore_threshold,
  "min_gene_performance" = gene_performance_threshold,
  "min_parliament_size" = parliament_size_threshold)

# Save
output.file <- file("./calibration_files/calibration_filtering_thresholds.csv", "wb") # Use this construct to save in Unix file format
write.table(filtering_thresholds, file = output.file, sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
close(output.file)


# 5.3 Confidence support
#------------------------

# Create data frame only holding bins
bins_df <- data_frame(range_support_main_id = 
                        cut(seq(0,100,100/bins)[-1], # Get the upper ends of the bins (exclude the first value, which is 0)
                            breaks = seq(0,100,100/bins), # Assign each of the values to its bin
                            include.lowest = TRUE)) # Include 0 in the first bin

# Assign top identifications to bins
confidence_support <- top_id %>% 
  mutate(range_support_main_id = cut(support_identification_pct, breaks = seq(0,100,100/bins), include.lowest = TRUE)) %>% # Add column that defines in which bin data are
  group_by(range_support_main_id) %>%
  reframe(count_all = n_distinct(query),
          count_correct = sum(id_correct_close == "correct"),
          count_close = sum(id_correct_close == "close"),
          count_wrong = sum(id_correct_close == "wrong"),
          percentage_correct = round((count_correct / count_all * 100), 2),
          percentage_close = round((count_close / count_all  * 100), 2),
          percentage_wrong = round((count_wrong / count_all * 100), 2)
  ) %>%
  select(range_support_main_id, percentage_correct, percentage_close, percentage_wrong)

# Combine bins data frame with binned top identifications
confidence_support <- left_join(bins_df, confidence_support, by = "range_support_main_id")

# Rename columns
confidence_support <- confidence_support %>%
  rename(range_support = range_support_main_id,
         probability_correct = percentage_correct,
         probability_close = percentage_close,
         probability_wrong = percentage_wrong)

# Save
output.file <- file("./calibration_files/calibration_confidence_support.csv", "wb") # Use this construct to save in Unix file format
write.table(confidence_support, file = output.file, sep = ",", col.names = TRUE, row.names = FALSE)
close(output.file)
