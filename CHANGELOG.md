# GPID release notes

## 1.2.0

### Command and file changes

- `gpid validate` is now `gpid confidence`; `gpid validate confidence` is now
  `gpid confidence estimate`. Confidence scripts and output paths use the new names.
- Calibration produces `gene_performance.csv` and `thresholds_filtering.csv`
  directly in the working directory. Confidence binning produces
  `confidence_support.csv` there. These outputs have no output-directory override.
- `thresholds_filtering.csv` uses `parameter,value` columns with one threshold per
  row. Identification and confidence estimation also accept legacy CSVs with
  eight parameter columns and one value row. Readers resolve values by name.
- Identification defaults to one shared `identifications` directory with
  sample-prefixed filenames and no saved `genelist_high_performance.txt`.

### Calibration and confidence

- Calibration accepts selected thresholds as CLI flags while retaining the
  editable TSV workflow and default threshold-file paths.
- Calibration parameter definitions are embedded in the scripts, removing the
  external template CSV dependency while keeping threshold selection unchanged.
- Calibration figures use automatic x ticks, visible grids, one Performance (%)
  axis, and a solid/dashed accuracy/retrievability legend.
- `gpid calibrate combine` also writes `thresholds_filtering.pdf`: eight panels in
  two columns and four rows, with selected thresholds marked in red.
- In `confidence_support.csv`, `probability_close` includes correct and close
  identifications. Plots retain separate correct, close, and wrong categories.
- Confidence bin figures have readable percentage intervals and rotated x labels.

### Input handling and reliability

- Reference FASTA headers exceeding 50 characters are rejected before BLAST
  database creation.
- Gene-performance CSVs accept extra columns and require unique `gene` and
  `performance` column names. Explicit NA performance generates a warning and
  is treated as zero for filtering. Calibration checks its generated CSV.
- Species-group checks ignore extra species and warn for unassigned reference
  species. Confidence and identification use genus grouping when no file is given.
- Identification consumes complete BLAST streams, checks gene/database paths,
  reports pipeline failures, and publishes the BLAST table only after success.
- Help text and output announcements are consistent; scripts include function
  descriptions and workflow annotations.

### Packaging and migration

- Conda metadata uses Apache-2.0 to match the bundled LICENSE, and declares the
  Bash and R package versions needed by the current code.
- R >=4.3 records the development baseline. Bash and coreutils are declared as
  build tools, and installed-package checks cover external tools and JPG/SVG
  rendering as well as PDF generation.
- BLAST+ >=2.16.0 records the development baseline. Package checks build a small
  reference database and verify an exact match with `blastn -task megablast`.
- Replace old commands and paths in job scripts. Regenerate calibration and
  confidence outputs using the current workflow, particularly the long-format
  threshold CSV and the cumulative close probabilities. Previously generated
  files are not converted when this release is installed.
