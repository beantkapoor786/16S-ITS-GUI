# DADA2 16S rRNA Analysis Pipeline — Shiny GUI

A modern, dark-themed R Shiny application that walks you through the complete DADA2 pipeline for 16S rRNA amplicon sequence analysis — from raw FASTQ files to phyloseq visualizations.

## Requirements

- **R** ≥ 4.0
- Internet connection (first run only — for package installation and SILVA database download)

All R packages are installed automatically on first launch:

| Package | Source | Purpose |
|---------|--------|---------|
| dada2 | Bioconductor | Core denoising pipeline |
| phyloseq | Bioconductor | Downstream community analysis |
| Biostrings | Bioconductor | DNA sequence handling |
| DECIPHER | Bioconductor | IdTaxa taxonomy assignment (optional method) |
| ggplot2 | CRAN | Plotting |
| tidyverse | CRAN | Data wrangling |
| data.table | CRAN | Fast data manipulation |
| shiny | CRAN | Web application framework |
| shinyjs | CRAN | JavaScript helpers |
| shinyWidgets | CRAN | Enhanced UI widgets |
| DT | CRAN | Interactive tables |
| shinycssloaders | CRAN | Loading spinners |

## How to Launch

```r
# From R console:
shiny::runApp("app.R")

# Or from terminal:
Rscript -e "shiny::runApp('app.R', launch.browser = TRUE)"
```

The app opens in your default browser at `http://127.0.0.1:XXXX`.

## Pipeline Steps (7-step workflow)

### Step 1 — Setup & Files
- Enter the path to your directory of demultiplexed paired-end FASTQ files
- The app auto-detects common forward/reverse naming patterns (e.g., `_R1_001.fastq`) with manual override
- Define sample name extraction: choose a delimiter and element index

### Step 2 — Quality Profiles
- View `plotQualityProfile` for forward and reverse reads
- Use these plots to decide truncation lengths for filtering

### Step 3 — Filter & Trim
- Set `truncLen`, `maxEE`, `truncQ`, `maxN`, `rm.phix` — all editable, tutorial defaults pre-filled
- Filtered files are saved to a `filtered/` subdirectory inside your data path
- View per-sample read retention table

### Step 4 — Learn Errors & Denoise
- Learns error models (`learnErrors`) for forward and reverse reads
- Runs core DADA2 sample inference (`dada`)
- Visualize error rate plots as sanity check
- Multithread auto-detected (enabled on Mac/Linux, disabled on Windows)

### Step 5 — Merge & Chimera Removal
- Merges paired reads (`mergePairs`)
- Constructs ASV sequence table (`makeSequenceTable`)
- Removes chimeras (`removeBimeraDenovo`)
- Shows sequence length distribution, read tracking table, and tracking line plot

### Step 6 — Taxonomy Assignment
- Downloads the SILVA v138.1 database automatically (stored locally for reuse)
- Choose between **assignTaxonomy** (Naive Bayesian) or **IdTaxa** (DECIPHER)
- Optional species-level assignment via `addSpecies`
- For IdTaxa: place `SILVA_SSU_r138_2024.RData` in the database directory

### Step 7 — Phyloseq & Export
- Optionally upload sample metadata (CSV or TSV, first column = sample names as row names)
- Builds a phyloseq object with ASV table, taxonomy, and metadata
- Visualizations: Alpha Diversity (Shannon/Simpson), NMDS Ordination (Bray-Curtis), Taxonomic Bar Plots
- Download all results: ASV table (CSV), taxonomy (CSV), read tracking (CSV), full R workspace (.RData)

## Input Data Requirements

Your FASTQ files must be:
1. **Demultiplexed** — one file per sample per direction
2. **Primers/adapters removed** — no non-biological nucleotides
3. **Matched order** — forward and reverse files correspond to same samples

Supported extensions: `.fastq`, `.fastq.gz`, `.fq`, `.fq.gz`

## Notes

- First run takes longer due to package installation and database downloads
- Subsequent runs load packages instantly — no reinstallation
- The SILVA database (~130 MB) is downloaded once to your specified directory
- All plots use a dark theme matching the application UI
