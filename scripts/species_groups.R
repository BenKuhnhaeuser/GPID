# --- Purpose ---
# Shared reference-species lookup and warnings for missing species-group assignments.
# Only reference species without a usable assignment need attention.
# missing_species_groups(): Return reference species with no nonempty, non-NA group assignment.
missing_species_groups <- function(reference_species, species_groups) {
  assigned <- !is.na(species_groups$species_group) &
    nzchar(trimws(species_groups$species_group)) &
    trimws(species_groups$species_group) != "NA"
  sort(setdiff(unique(reference_species[!is.na(reference_species)]),
               species_groups$genus_species[assigned]))
}

# warn_missing_species_groups(): Warn once per unassigned reference species and return the missing-species list.
warn_missing_species_groups <- function(reference_species, species_groups) {
  missing <- missing_species_groups(reference_species, species_groups)
  for (species in missing) {
    warning(sprintf("No species group was provided for reference species: %s", species), call. = FALSE)
  }
  invisible(missing)
}

# reference_fasta_species(): Collect unique reference species across supported FASTA files.
reference_fasta_species <- function(reference_dir) {
  if (!dir.exists(reference_dir)) stop("Reference directory not found: ", reference_dir, call. = FALSE)
  files <- list.files(reference_dir, pattern = "\\.(fna|fasta|fa)$", ignore.case = TRUE, full.names = TRUE)
  if (!length(files)) stop("No reference FASTA files found in: ", reference_dir, call. = FALSE)
  # read_species(): Read one FASTA in bounded chunks and extract species from its header lines.
  read_species <- function(file) {
    connection <- file(file, "r")
    on.exit(close(connection))
    species <- character()
    repeat {
      lines <- readLines(connection, n = 10000L, warn = FALSE)
      if (!length(lines)) break
      headers <- lines[grepl("^>[^_]+_[^_]+", lines)]
      species <- union(species, sub("^>([^_]+_[^_]+).*", "\\1", headers))
    }
    species
  }
  sort(unique(unlist(lapply(files, read_species), use.names = FALSE)))
}
