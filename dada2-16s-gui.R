# ============================================================================
# DADA2 16S rRNA Analysis Pipeline - Shiny GUI
# Full pipeline: QC -> Filter -> Dereplicate -> Merge -> Chimera -> Taxonomy -> Phyloseq
# ============================================================================

# ── Package Management ───────────────────────────────────────────────────────

required_cran <- c("shiny", "ggplot2", "tidyverse", "shinyjs",
                   "DT", "shinycssloaders", "callr", "vegan")
required_bioc <- c("dada2", "phyloseq", "Biostrings", "DECIPHER", "microbiome", "ANCOMBC")

install_if_missing <- function() {
  # CRAN packages
  missing_cran <- required_cran[!sapply(required_cran, requireNamespace, quietly = TRUE)]
  if (length(missing_cran) > 0) {
    install.packages(missing_cran, repos = "https://cloud.r-project.org")
  }
  # Bioconductor packages
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  missing_bioc <- required_bioc[!sapply(required_bioc, requireNamespace, quietly = TRUE)]
  if (length(missing_bioc) > 0) {
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
  }
}

install_if_missing()

library(shiny)
library(shinyjs)
library(DT)
library(shinycssloaders)
library(ggplot2)
library(tidyverse)
library(dada2)
library(phyloseq)
library(Biostrings)
library(DECIPHER)
library(callr)
library(vegan)
library(microbiome)
library(ANCOMBC)

# ── Multithread Strategy ─────────────────────────────────────────────────────
# DADA2's multithread=TRUE uses fork-based parallelism (mclapply) which
# conflicts with Shiny's reactive context on macOS/Linux.
#
# Solution: Heavy compute steps (filterAndTrim, learnErrors, dada,
# removeBimeraDenovo, assignTaxonomy) are offloaded to a separate R
# subprocess via callr::r(). The subprocess has no Shiny reactives,
# so multithread=TRUE works safely. Lightweight steps stay in-process
# with multithread=FALSE.
can_multithread <- .Platform$OS.type != "windows"

# ── SILVA database URLs (v138.1) ──────────────────────────────────────────

SILVA_GENUS_URL <- "https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz"
SILVA_SPECIES_URL <- "https://zenodo.org/records/4587955/files/silva_species_assignment_v138.1.fa.gz"

# ── Common forward/reverse patterns ─────────────────────────────────────────

common_fwd_patterns <- c("_R1_001.fastq", "_R1_001.fastq.gz", "_R1.fastq",
                         "_R1.fastq.gz", "_R1_001.fq", "_R1_001.fq.gz",
                         "_1.fastq", "_1.fastq.gz", "_1.fq", "_1.fq.gz")
common_rev_patterns <- c("_R2_001.fastq", "_R2_001.fastq.gz", "_R2.fastq",
                         "_R2.fastq.gz", "_R2_001.fq", "_R2_001.fq.gz",
                         "_2.fastq", "_2.fastq.gz", "_2.fq", "_2.fq.gz")

# ── Custom CSS ───────────────────────────────────────────────────────────────

app_css <- "
@import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;600&family=Outfit:wght@300;400;500;600;700&display=swap');

:root {
  --bg-primary: #0a0e17;
  --bg-secondary: #111827;
  --bg-card: #1a2332;
  --bg-card-hover: #1f2b3d;
  --bg-input: #0d1320;
  --border-color: #2a3a52;
  --border-accent: #3b82f6;
  --text-primary: #e8ecf4;
  --text-secondary: #8899b0;
  --text-muted: #5a6a80;
  --accent-blue: #3b82f6;
  --accent-cyan: #06b6d4;
  --accent-emerald: #10b981;
  --accent-amber: #f59e0b;
  --accent-rose: #f43f5e;
  --accent-violet: #8b5cf6;
  --gradient-1: linear-gradient(135deg, #3b82f6 0%, #06b6d4 100%);
  --gradient-2: linear-gradient(135deg, #8b5cf6 0%, #ec4899 100%);
  --gradient-3: linear-gradient(135deg, #10b981 0%, #3b82f6 100%);
  --shadow-sm: 0 1px 3px rgba(0,0,0,0.4);
  --shadow-md: 0 4px 20px rgba(0,0,0,0.5);
  --shadow-lg: 0 10px 40px rgba(0,0,0,0.6);
  --radius: 12px;
  --radius-sm: 8px;
  --radius-lg: 16px;
}

* { box-sizing: border-box; margin: 0; padding: 0; }

body {
  background: var(--bg-primary) !important;
  color: var(--text-primary) !important;
  font-family: 'Outfit', sans-serif !important;
  font-weight: 400;
  line-height: 1.6;
  min-height: 100vh;
}

/* ── Scrollbar ── */
::-webkit-scrollbar { width: 6px; height: 6px; }
::-webkit-scrollbar-track { background: var(--bg-primary); }
::-webkit-scrollbar-thumb { background: var(--border-color); border-radius: 3px; }
::-webkit-scrollbar-thumb:hover { background: var(--accent-blue); }

/* ── Header ── */
.app-header {
  background: var(--bg-secondary);
  border-bottom: 1px solid var(--border-color);
  padding: 16px 32px;
  display: flex;
  align-items: center;
  gap: 16px;
  position: sticky;
  top: 0;
  z-index: 1000;
  overflow: hidden;
}

.app-header-bacteria {
  flex: 1;
  height: 50px;
  opacity: 0.07;
  pointer-events: none;
  margin-left: auto;
}

.app-logo {
  width: 42px; height: 42px;
  background: var(--gradient-1);
  border-radius: 10px;
  display: flex; align-items: center; justify-content: center;
  font-family: 'JetBrains Mono', monospace;
  font-weight: 700; font-size: 16px; color: white;
  flex-shrink: 0;
  box-shadow: 0 0 20px rgba(59,130,246,0.3);
}

.app-title {
  font-size: 22px; font-weight: 600;
  background: var(--gradient-1);
  -webkit-background-clip: text; -webkit-text-fill-color: transparent;
  background-clip: text;
}

.app-subtitle {
  font-size: 13px; color: var(--text-muted); font-weight: 300;
  margin-top: -2px;
}

/* ── Stepper ── */
.stepper-container {
  background: var(--bg-secondary);
  border-bottom: 1px solid var(--border-color);
  padding: 12px 32px;
  display: flex;
  gap: 4px;
  overflow-x: auto;
}

.step-pill {
  padding: 8px 18px;
  border-radius: 100px;
  font-size: 13px;
  font-weight: 500;
  white-space: nowrap;
  cursor: pointer;
  transition: all 0.25s ease;
  border: 1px solid transparent;
  color: var(--text-muted);
  background: transparent;
  display: flex; align-items: center; gap: 8px;
}

.step-pill:hover { color: var(--text-secondary); background: var(--bg-card); }
.step-pill.active {
  background: var(--accent-blue);
  color: white;
  box-shadow: 0 0 20px rgba(59,130,246,0.3);
}
.step-pill.completed {
  color: var(--accent-emerald);
  border-color: rgba(16,185,129,0.3);
  background: rgba(16,185,129,0.08);
}

.step-number {
  width: 22px; height: 22px;
  border-radius: 50%;
  display: flex; align-items: center; justify-content: center;
  font-size: 11px; font-weight: 600;
  background: var(--bg-card);
  border: 1px solid var(--border-color);
}

.step-pill.active .step-number { background: rgba(255,255,255,0.2); border-color: transparent; }
.step-pill.completed .step-number { background: var(--accent-emerald); color: white; border-color: transparent; }

/* ── Main container ── */
.main-container {
  max-width: 1400px;
  margin: 0 auto;
  padding: 24px 32px 60px;
}

/* ── Cards ── */
.card {
  background: var(--bg-card);
  border: 1px solid var(--border-color);
  border-radius: var(--radius);
  padding: 24px;
  margin-bottom: 20px;
  box-shadow: var(--shadow-sm);
  transition: border-color 0.2s ease;
}

.card:hover { border-color: rgba(59,130,246,0.3); }

.card-header {
  font-size: 16px; font-weight: 600;
  margin-bottom: 16px;
  display: flex; align-items: center; gap: 10px;
  color: var(--text-primary);
}

.card-header .icon {
  width: 32px; height: 32px;
  border-radius: var(--radius-sm);
  display: flex; align-items: center; justify-content: center;
  font-size: 16px;
  flex-shrink: 0;
}

.card-header .icon.blue { background: rgba(59,130,246,0.15); color: var(--accent-blue); }
.card-header .icon.cyan { background: rgba(6,182,212,0.15); color: var(--accent-cyan); }
.card-header .icon.emerald { background: rgba(16,185,129,0.15); color: var(--accent-emerald); }
.card-header .icon.amber { background: rgba(245,158,11,0.15); color: var(--accent-amber); }
.card-header .icon.rose { background: rgba(244,63,94,0.15); color: var(--accent-rose); }
.card-header .icon.violet { background: rgba(139,92,246,0.15); color: var(--accent-violet); }

.card-description {
  font-size: 13px; color: var(--text-secondary);
  margin-bottom: 16px; line-height: 1.5;
}

/* ── Inputs ── */
.shiny-input-container { width: 100% !important; }

.form-control, .shiny-input-container input[type='text'],
.shiny-input-container input[type='number'],
.selectize-input, .selectize-control.single .selectize-input {
  background: var(--bg-input) !important;
  border: 1px solid var(--border-color) !important;
  border-radius: var(--radius-sm) !important;
  color: var(--text-primary) !important;
  font-family: 'JetBrains Mono', monospace !important;
  font-size: 13px !important;
  padding: 10px 14px !important;
  transition: border-color 0.2s ease !important;
}

.form-control:focus, .shiny-input-container input:focus, .selectize-input.focus {
  border-color: var(--accent-blue) !important;
  box-shadow: 0 0 0 3px rgba(59,130,246,0.15) !important;
  outline: none !important;
}

.selectize-dropdown {
  background: var(--bg-card) !important;
  border: 1px solid var(--border-color) !important;
  border-radius: var(--radius-sm) !important;
}

.selectize-dropdown .option {
  color: var(--text-primary) !important;
  font-size: 13px !important;
}

.selectize-dropdown .option.active {
  background: var(--accent-blue) !important;
  color: white !important;
}

/* Fix selectize inner elements */
.selectize-input > .item {
  color: var(--text-primary) !important;
  background: transparent !important;
  border: none !important;
  box-shadow: none !important;
}
.selectize-input::after {
  border-color: var(--text-muted) transparent transparent transparent !important;
}
.selectize-control.single .selectize-input::after {
  right: 12px !important;
}
.selectize-input .remove-single {
  display: none !important;
}
.selectize-dropdown-content {
  max-height: 250px !important;
  overflow-y: auto !important;
}

.control-label, label {
  color: var(--text-secondary) !important;
  font-size: 12px !important;
  font-weight: 500 !important;
  text-transform: uppercase !important;
  letter-spacing: 0.05em !important;
  margin-bottom: 6px !important;
}

.checkbox label, .radio label {
  text-transform: none !important;
  letter-spacing: normal !important;
  font-size: 13px !important;
  color: var(--text-primary) !important;
}

/* ── File input (Browse button on right) ── */
.shiny-input-container .input-group {
  display: flex !important;
  align-items: stretch !important;
  flex-direction: row !important;
  width: 100% !important;
  box-sizing: border-box !important;
}
/* Shiny puts btn-group first in DOM, text input second. We reorder with flex. */
.shiny-input-container .input-group .input-group-btn {
  order: 2 !important;
  display: flex !important;
  flex: 0 0 auto !important;
  z-index: 2 !important;
}
.shiny-input-container .input-group .form-control[readonly] {
  background: var(--bg-input) !important;
  border: 1px solid var(--border-color) !important;
  border-right: none !important;
  border-radius: var(--radius-sm) 0 0 var(--radius-sm) !important;
  color: var(--text-secondary) !important;
  font-family: 'JetBrains Mono', monospace !important;
  font-size: 13px !important;
  cursor: default !important;
  order: 1 !important;
  flex: 1 1 100px !important;
  min-width: 100px !important;
  overflow: hidden !important;
  text-overflow: ellipsis !important;
  box-sizing: border-box !important;
}
.shiny-input-container .input-group .btn-default,
.shiny-input-container .input-group .btn-file {
  background: var(--gradient-1) !important;
  border: none !important;
  border-radius: 0 var(--radius-sm) var(--radius-sm) 0 !important;
  color: white !important;
  font-family: 'Outfit', sans-serif !important;
  font-weight: 600 !important;
  font-size: 13px !important;
  padding: 10px 24px !important;
  cursor: pointer !important;
  transition: all 0.2s ease !important;
  white-space: nowrap !important;
  overflow: visible !important;
  min-width: fit-content !important;
}
.shiny-input-container .input-group .btn-default:hover,
.shiny-input-container .input-group .btn-file:hover {
  box-shadow: 0 2px 10px rgba(59,130,246,0.3) !important;
}
.shiny-input-container .progress {
  margin-top: 6px !important;
}

/* ── Buttons ── */
.btn-primary, .btn-run {
  background: var(--gradient-1) !important;
  border: none !important;
  border-radius: var(--radius-sm) !important;
  color: white !important;
  font-family: 'Outfit', sans-serif !important;
  font-weight: 600 !important;
  font-size: 14px !important;
  padding: 12px 28px !important;
  cursor: pointer !important;
  transition: all 0.25s ease !important;
  box-shadow: 0 4px 15px rgba(59,130,246,0.3) !important;
  letter-spacing: 0.02em !important;
}

.btn-primary:hover, .btn-run:hover {
  box-shadow: 0 6px 25px rgba(59,130,246,0.5) !important;
  transform: translateY(-1px) !important;
}

.btn-secondary {
  background: var(--bg-card) !important;
  border: 1px solid var(--border-color) !important;
  border-radius: var(--radius-sm) !important;
  color: var(--text-primary) !important;
  font-family: 'Outfit', sans-serif !important;
  font-weight: 500 !important;
  font-size: 13px !important;
  padding: 10px 20px !important;
  cursor: pointer;
  transition: all 0.2s ease !important;
}

.btn-secondary:hover {
  border-color: var(--accent-blue) !important;
  background: var(--bg-card-hover) !important;
}

.btn-download {
  background: rgba(16,185,129,0.12) !important;
  border: 1px solid rgba(16,185,129,0.3) !important;
  border-radius: var(--radius-sm) !important;
  color: var(--accent-emerald) !important;
  font-family: 'Outfit', sans-serif !important;
  font-weight: 500 !important;
  font-size: 13px !important;
  padding: 10px 20px !important;
  cursor: pointer;
  transition: all 0.2s ease !important;
}

.btn-download:hover {
  background: rgba(16,185,129,0.2) !important;
  border-color: var(--accent-emerald) !important;
}

/* ── Tables ── */
.dataTables_wrapper { color: var(--text-primary) !important; }
table.dataTable { border-collapse: collapse !important; }
table.dataTable thead th {
  background: var(--bg-secondary) !important;
  color: var(--text-secondary) !important;
  border-bottom: 1px solid var(--border-color) !important;
  font-size: 12px !important;
  font-weight: 600 !important;
  text-transform: uppercase !important;
  letter-spacing: 0.05em !important;
  padding: 12px 16px !important;
}
table.dataTable tbody td {
  background: var(--bg-card) !important;
  color: var(--text-primary) !important;
  border-bottom: 1px solid rgba(42,58,82,0.5) !important;
  padding: 10px 16px !important;
  font-size: 13px !important;
  font-family: 'JetBrains Mono', monospace !important;
}
table.dataTable tbody tr:hover td { background: var(--bg-card-hover) !important; }

.dataTables_info, .dataTables_length, .dataTables_filter,
.dataTables_length label, .dataTables_filter label,
.dataTables_paginate { color: var(--text-muted) !important; font-size: 12px !important; }

.dataTables_filter input {
  background: var(--bg-input) !important;
  border: 1px solid var(--border-color) !important;
  color: var(--text-primary) !important;
  border-radius: var(--radius-sm) !important;
  padding: 6px 10px !important;
}

.paginate_button { color: var(--text-muted) !important; }
.paginate_button.current { background: var(--accent-blue) !important; color: white !important; border-radius: 6px !important; }

/* ── Status / Log ── */
.log-panel {
  background: var(--bg-input);
  border: 1px solid var(--border-color);
  border-radius: var(--radius-sm);
  padding: 16px;
  max-height: 300px;
  overflow-y: auto;
  font-family: 'JetBrains Mono', monospace;
  font-size: 12px;
  line-height: 1.8;
  color: var(--text-secondary);
}

.log-panel .log-success { color: var(--accent-emerald); }
.log-panel .log-error { color: var(--accent-rose); }
.log-panel .log-info { color: var(--accent-cyan); }
.log-panel .log-warn { color: var(--accent-amber); }

/* ── Status badge ── */
.status-badge {
  display: inline-flex; align-items: center; gap: 6px;
  padding: 6px 14px;
  border-radius: 100px;
  font-size: 12px; font-weight: 600;
  letter-spacing: 0.03em;
}

.status-badge.running { background: rgba(59,130,246,0.15); color: var(--accent-blue); }
.status-badge.success { background: rgba(16,185,129,0.15); color: var(--accent-emerald); }
.status-badge.error { background: rgba(244,63,94,0.15); color: var(--accent-rose); }
.status-badge.waiting { background: rgba(90,106,128,0.15); color: var(--text-muted); }

/* ── Plots ── */
.shiny-plot-output {
  border-radius: var(--radius-sm);
  overflow: hidden;
}

/* ── Spinner override ── */
.load-container .shiny-spinner-output-container .fa-spinner { color: var(--accent-blue) !important; }

/* ── Tab panels ── */
.nav-tabs {
  border-bottom: 1px solid var(--border-color) !important;
  margin-bottom: 16px !important;
}

.nav-tabs > li > a {
  color: var(--text-muted) !important;
  background: transparent !important;
  border: none !important;
  border-bottom: 2px solid transparent !important;
  font-size: 13px !important;
  font-weight: 500 !important;
  padding: 10px 18px !important;
  transition: all 0.2s ease !important;
}

.nav-tabs > li > a:hover {
  color: var(--text-primary) !important;
  border-bottom-color: var(--border-color) !important;
  background: transparent !important;
}

.nav-tabs > li.active > a, .nav-tabs > li.active > a:focus, .nav-tabs > li.active > a:hover {
  color: var(--accent-blue) !important;
  background: transparent !important;
  border: none !important;
  border-bottom: 2px solid var(--accent-blue) !important;
}

.tab-content { background: transparent !important; }

/* ── Two column grid ── */
.grid-2 { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }
.grid-3 { display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 16px; }
.grid-4 { display: grid; grid-template-columns: repeat(4, 1fr); gap: 16px; }

@media (max-width: 900px) {
  .grid-2, .grid-3, .grid-4 { grid-template-columns: 1fr; }
}

/* ── Well panels ── */
.well {
  background: var(--bg-card) !important;
  border: 1px solid var(--border-color) !important;
  border-radius: var(--radius) !important;
  box-shadow: none !important;
}

/* ── Progress bar ── */
.progress {
  background: var(--bg-input) !important;
  border-radius: 100px !important;
  height: 8px !important;
  overflow: hidden;
}
.progress-bar {
  background: var(--gradient-1) !important;
  border-radius: 100px !important;
}

/* ── Hide default shiny error styling ── */
.shiny-output-error { color: var(--accent-rose) !important; }
.shiny-output-error:before { content: '' !important; }

/* ── Panel animations ── */
.step-panel { animation: fadeIn 0.3s ease; }
@keyframes fadeIn { from { opacity: 0; transform: translateY(8px); } to { opacity: 1; transform: translateY(0); } }

/* ── Proceed button ── */
.proceed-bar {
  margin-top: 24px;
  padding-top: 20px;
  border-top: 1px solid var(--border-color);
  display: flex;
  align-items: center;
  justify-content: flex-end;
  gap: 12px;
}
.btn-proceed {
  background: var(--gradient-3) !important;
  border: none !important;
  border-radius: var(--radius-sm) !important;
  color: white !important;
  font-family: 'Outfit', sans-serif !important;
  font-weight: 600 !important;
  font-size: 14px !important;
  padding: 12px 32px !important;
  cursor: pointer !important;
  transition: all 0.25s ease !important;
  box-shadow: 0 4px 15px rgba(16,185,129,0.3) !important;
  letter-spacing: 0.02em !important;
}
.btn-proceed:hover {
  box-shadow: 0 6px 25px rgba(16,185,129,0.5) !important;
  transform: translateY(-1px) !important;
}
.btn-proceed:disabled, .btn-proceed.disabled {
  opacity: 0.35 !important;
  cursor: not-allowed !important;
  transform: none !important;
  box-shadow: none !important;
}

/* ── Step progress indicator ── */
.step-progress-container {
  margin: 16px 0;
  animation: fadeIn 0.3s ease;
}
.step-progress-label {
  font-size: 13px;
  font-weight: 500;
  color: var(--text-secondary);
  margin-bottom: 8px;
  display: flex;
  align-items: center;
  justify-content: space-between;
}
.step-progress-label .progress-status {
  font-family: 'JetBrains Mono', monospace;
  font-size: 12px;
  color: var(--accent-cyan);
}
.step-progress-track {
  width: 100%;
  height: 10px;
  background: var(--bg-input);
  border: 1px solid var(--border-color);
  border-radius: 100px;
  overflow: hidden;
  position: relative;
}
.step-progress-fill {
  height: 100%;
  border-radius: 100px;
  background: var(--gradient-1);
  transition: width 0.4s ease;
  position: relative;
}
.step-progress-fill.indeterminate {
  width: 30% !important;
  animation: indeterminate 1.8s ease-in-out infinite;
}
@keyframes indeterminate {
  0% { transform: translateX(-100%); }
  100% { transform: translateX(400%); }
}
.step-progress-fill.complete {
  background: var(--gradient-3) !important;
  width: 100% !important;
}

/* ── Running indicator (shows instantly) ── */
.running-indicator {
  display: flex;
  align-items: center;
  gap: 12px;
  margin-top: 16px;
  padding: 14px 18px;
  background: rgba(59,130,246,0.08);
  border: 1px solid rgba(59,130,246,0.25);
  border-radius: var(--radius-sm);
  animation: fadeIn 0.2s ease;
}
.running-indicator .pulse-dot {
  width: 10px; height: 10px;
  border-radius: 50%;
  background: var(--accent-blue);
  animation: pulse 1.2s ease-in-out infinite;
  flex-shrink: 0;
}
@keyframes pulse {
  0%, 100% { opacity: 1; transform: scale(1); }
  50% { opacity: 0.4; transform: scale(0.75); }
}
.running-indicator .running-text {
  font-size: 13px;
  font-weight: 500;
  color: var(--accent-blue);
}
.running-indicator .running-sub {
  font-size: 11px;
  color: var(--text-muted);
  font-family: 'JetBrains Mono', monospace;
}

/* ── Stat cards ── */
.stat-card {
  background: var(--bg-card);
  border: 1px solid var(--border-color);
  border-radius: var(--radius);
  padding: 18px 20px;
  text-align: center;
}
.stat-card .stat-value {
  font-size: 28px; font-weight: 700;
  background: var(--gradient-1);
  -webkit-background-clip: text; -webkit-text-fill-color: transparent;
  background-clip: text;
}
.stat-card .stat-label {
  font-size: 11px; color: var(--text-muted);
  text-transform: uppercase; letter-spacing: 0.08em;
  margin-top: 4px; font-weight: 500;
}

/* ── Numeric input spinner fix ── */
input[type=number] { -moz-appearance: textfield; }
input[type=number]::-webkit-inner-spin-button,
input[type=number]::-webkit-outer-spin-button { opacity: 1; }

/* ── Slider override ── */
.irs--shiny .irs-bar { background: var(--accent-blue) !important; border-top: none !important; border-bottom: none !important; }
.irs--shiny .irs-handle { border-color: var(--accent-blue) !important; background: var(--accent-blue) !important; }
.irs--shiny .irs-from, .irs--shiny .irs-to, .irs--shiny .irs-single {
  background: var(--accent-blue) !important; font-size: 11px !important;
}
.irs--shiny .irs-line { background: var(--bg-input) !important; border: 1px solid var(--border-color) !important; }
.irs--shiny .irs-grid-text { color: var(--text-muted) !important; font-size: 10px !important; }

/* ── Modal styling ── */
.modal-content {
  background: var(--bg-card) !important;
  border: 1px solid var(--border-color) !important;
  border-radius: var(--radius-lg) !important;
  box-shadow: var(--shadow-lg) !important;
  color: var(--text-primary) !important;
}
.modal-header {
  border-bottom: 1px solid var(--border-color) !important;
  padding: 20px 24px !important;
}
.modal-header .modal-title {
  font-family: 'Outfit', sans-serif !important;
  font-weight: 600 !important;
  font-size: 18px !important;
  color: var(--text-primary) !important;
}
.modal-header .close {
  color: var(--text-muted) !important;
  opacity: 0.8 !important;
  text-shadow: none !important;
}
.modal-body {
  padding: 24px !important;
  color: var(--text-secondary) !important;
  font-size: 14px !important;
  line-height: 1.6 !important;
}
.modal-footer {
  border-top: 1px solid var(--border-color) !important;
  padding: 16px 24px !important;
  display: flex !important;
  gap: 12px !important;
  justify-content: flex-end !important;
}
.modal-backdrop { background: #000 !important; }
.modal-backdrop.in { opacity: 0.7 !important; }

.session-info-card {
  background: var(--bg-input);
  border: 1px solid var(--border-color);
  border-radius: var(--radius-sm);
  padding: 16px;
  margin: 12px 0;
  font-family: 'JetBrains Mono', monospace;
  font-size: 13px;
  line-height: 1.8;
}
.session-info-card .info-label {
  color: var(--text-muted);
  font-size: 11px;
  text-transform: uppercase;
  letter-spacing: 0.05em;
}
.session-info-card .info-value {
  color: var(--accent-cyan);
}

/* ── Rerun warning banner ── */
.rerun-warning {
  background: rgba(245,158,11,0.1);
  border: 1px solid rgba(245,158,11,0.3);
  border-radius: var(--radius-sm);
  padding: 12px 16px;
  margin-bottom: 16px;
  font-size: 13px;
  color: var(--accent-amber);
  display: flex;
  align-items: center;
  gap: 10px;
}
"

# ══════════════════════════════════════════════════════════════════════════════
# UI
# ══════════════════════════════════════════════════════════════════════════════

ui <- fluidPage(
  useShinyjs(),
  tags$head(
    tags$style(HTML(app_css)),
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$script(HTML("
      // Instant running indicator: inject DOM element immediately on click,
      // before Shiny's reactive flush. This ensures the user sees feedback instantly
      // even if the observeEvent handler has synchronous blocking work.
      $(document).on('click', '#btn_filter, #btn_denoise, #btn_merge, #btn_taxonomy, #btn_ancombc2', function() {
        var btnId = $(this).attr('id');
        var targetMap = {
          'btn_filter': 'progress_step3',
          'btn_denoise': 'progress_step4',
          'btn_merge': 'progress_step5',
          'btn_taxonomy': 'progress_step6',
          'btn_ancombc2': 'progress_step10'
        };
        var target = targetMap[btnId];
        if (target) {
          var el = document.getElementById(target);
          if (el) {
            el.innerHTML = '<div class=\"step-progress-container\">' +
              '<div class=\"step-progress-label\">' +
                '<span>\\u23f3  Initializing...</span>' +
                '<span class=\"progress-status\"></span>' +
              '</div>' +
              '<div class=\"step-progress-track\">' +
                '<div class=\"step-progress-fill indeterminate\" style=\"width:30%;\"></div>' +
              '</div></div>';
          }
        }
      });
    "))
  ),

  # ── Header ──
  div(class = "app-header",
    div(class = "app-logo", "16S"),
    div(
      div(class = "app-title", "DADA2 Pipeline"),
      div(class = "app-subtitle", "16S rRNA Amplicon Sequence Analysis")
    ),
    HTML('<svg class="app-header-bacteria" viewBox="0 0 900 80" preserveAspectRatio="xMaxYMid meet" xmlns="http://www.w3.org/2000/svg">
      <g fill="#e8ecf4" stroke="#e8ecf4">
        <!-- rod with flagella -->
        <ellipse cx="80" cy="40" rx="38" ry="14" transform="rotate(-12 80 40)" fill="#e8ecf4" stroke="none"/>
        <path d="M48 32 Q32 18 20 24 Q10 30 4 20" stroke-width="2.5" fill="none" stroke-linecap="round"/>
        <path d="M45 38 Q28 30 16 38 Q6 44 0 36" stroke-width="2" fill="none" stroke-linecap="round"/>
        <!-- coccus pair -->
        <circle cx="175" cy="32" r="12" fill="#e8ecf4" stroke="none"/>
        <circle cx="195" cy="28" r="9" fill="#e8ecf4" stroke="none"/>
        <!-- spirillum -->
        <path d="M240 42 Q260 22 280 42 Q300 62 320 42 Q340 22 360 42 Q380 62 400 42" stroke-width="7" fill="none" stroke-linecap="round"/>
        <!-- small rod -->
        <ellipse cx="460" cy="38" rx="30" ry="11" transform="rotate(8 460 38)" fill="#e8ecf4" stroke="none"/>
        <!-- cocci cluster -->
        <circle cx="540" cy="35" r="9" fill="#e8ecf4" stroke="none"/>
        <circle cx="555" cy="28" r="7" fill="#e8ecf4" stroke="none"/>
        <circle cx="548" cy="48" r="8" fill="#e8ecf4" stroke="none"/>
        <!-- another rod with flagella -->
        <ellipse cx="630" cy="42" rx="35" ry="12" transform="rotate(-18 630 42)" fill="#e8ecf4" stroke="none"/>
        <path d="M600 34 Q585 20 572 26" stroke-width="2" fill="none" stroke-linecap="round"/>
        <path d="M598 40 Q582 32 570 40" stroke-width="1.8" fill="none" stroke-linecap="round"/>
        <!-- vibrio -->
        <path d="M710 50 Q730 20 750 35 Q770 50 780 30" stroke-width="8" fill="none" stroke-linecap="round"/>
        <!-- small cocci -->
        <circle cx="830" cy="36" r="10" fill="#e8ecf4" stroke="none"/>
        <circle cx="850" cy="42" r="8" fill="#e8ecf4" stroke="none"/>
        <!-- tiny rod -->
        <ellipse cx="890" cy="35" rx="22" ry="8" transform="rotate(-5 890 35)" fill="#e8ecf4" stroke="none"/>
      </g>
    </svg>')
  ),

  # ── Step navigation ──
  div(class = "stepper-container",
    div(id = "step_nav_1", class = "step-pill active", onclick = "Shiny.setInputValue('nav_step', 1, {priority: 'event'})",
      span(class = "step-number", "1"), "Setup & Files"),
    div(id = "step_nav_2", class = "step-pill", onclick = "Shiny.setInputValue('nav_step', 2, {priority: 'event'})",
      span(class = "step-number", "2"), "Quality Profiles"),
    div(id = "step_nav_3", class = "step-pill", onclick = "Shiny.setInputValue('nav_step', 3, {priority: 'event'})",
      span(class = "step-number", "3"), "Filter & Trim"),
    div(id = "step_nav_4", class = "step-pill", onclick = "Shiny.setInputValue('nav_step', 4, {priority: 'event'})",
      span(class = "step-number", "4"), "Dereplication"),
    div(id = "step_nav_5", class = "step-pill", onclick = "Shiny.setInputValue('nav_step', 5, {priority: 'event'})",
      span(class = "step-number", "5"), "Merge & Chimeras"),
    div(id = "step_nav_6", class = "step-pill", onclick = "Shiny.setInputValue('nav_step', 6, {priority: 'event'})",
      span(class = "step-number", "6"), "Taxonomy"),
    div(id = "step_nav_7", class = "step-pill", onclick = "Shiny.setInputValue('nav_step', 7, {priority: 'event'})",
      span(class = "step-number", "7"), "Phyloseq"),
    div(id = "step_nav_8", class = "step-pill", onclick = "Shiny.setInputValue('nav_step', 8, {priority: 'event'})",
      span(class = "step-number", "8"), "Figures"),
    div(id = "step_nav_9", class = "step-pill", onclick = "Shiny.setInputValue('nav_step', 9, {priority: 'event'})",
      span(class = "step-number", "9"), "PERMANOVA"),
    div(id = "step_nav_10", class = "step-pill", onclick = "Shiny.setInputValue('nav_step', 10, {priority: 'event'})",
      span(class = "step-number", "10"), "ANCOM-BC2")
  ),

  # ── Main content ──
  div(class = "main-container",

    # ═══════════════════════════════════════════════════════════════════════
    # STEP 1: Setup & File Loading
    # ═══════════════════════════════════════════════════════════════════════
    div(id = "panel_step1", class = "step-panel",

      div(class = "card",
        div(class = "card-header",
          div(class = "icon blue", icon("folder-open")),
          "Sequence Data Directory"
        ),
        div(class = "card-description",
          "Provide the full path to the directory containing your demultiplexed, paired-end FASTQ files."
        ),
        div(class = "grid-2",
          div(
            textInput("data_path", "Path to FASTQ Directory", placeholder = "/path/to/your/fastq/files"),
            actionButton("btn_scan_files", "Scan Directory", class = "btn-primary",
                         icon = icon("magnifying-glass"))
          ),
          div(
            textInput("fwd_pattern", "Forward Read Pattern", value = ""),
            textInput("rev_pattern", "Reverse Read Pattern", value = "")
          )
        ),
        uiOutput("detected_pattern_ui")
      ),

      div(class = "card",
        div(class = "card-header",
          div(class = "icon cyan", icon("tag")),
          "Sample Name Extraction"
        ),
        div(class = "card-description",
          "Define how sample names are parsed from filenames. Filenames are split on the delimiter, then the selected element indices are joined to form the sample name. For example, with delimiter '_' and indices 1-2 on 'F3D0_S188_L001_R1_001.fastq', you get 'F3D0_S188'."
        ),
        div(class = "grid-4",
          textInput("sample_delim", "Delimiter", value = "_"),
          numericInput("sample_element_from", "From Element Index", value = 1, min = 1, step = 1),
          numericInput("sample_element_to", "To Element Index", value = 1, min = 1, step = 1),
          div(style = "display:flex; align-items:flex-end; height:100%;",
            actionButton("btn_extract_samples", "Extract Sample Names", class = "btn-primary",
                         icon = icon("flask"))
          )
        ),
        div(class = "card-description", style = "margin-top: 8px;",
          uiOutput("sample_name_preview")
        )
      ),

      div(class = "card",
        div(class = "card-header",
          div(class = "icon emerald", icon("list")),
          "Detected Files & Samples"
        ),
        uiOutput("file_summary_ui"),
        DTOutput("sample_table") %>% withSpinner(type = 6, color = "#3b82f6")
      ),

      div(class = "proceed-bar",
        actionButton("proceed_1", "Proceed to Quality Profiles ->", class = "btn-proceed",
                     icon = icon("arrow-right"))
      ),

      div(class = "log-panel", id = "log_step1", htmlOutput("log1"))
    ),

    # ═══════════════════════════════════════════════════════════════════════
    # STEP 2: Quality Profiles
    # ═══════════════════════════════════════════════════════════════════════
    hidden(div(id = "panel_step2", class = "step-panel",

      div(class = "card",
        div(class = "card-header",
          div(class = "icon cyan", icon("chart-line")),
          "Quality Profile Inspection"
        ),
        div(class = "card-description",
          "Visualize quality scores across read positions. Use these plots to decide truncation lengths for the filter step."
        ),
        div(class = "grid-2",
          numericInput("qp_n_samples", "Number of samples to plot", value = 2, min = 1, max = 20),
          div(style = "display:flex; align-items:flex-end;",
            actionButton("btn_plot_quality", "Generate Quality Plots", class = "btn-primary",
                         icon = icon("chart-area"))
          )
        )
      ),

      uiOutput("qplot_cards"),

      div(class = "proceed-bar",
        actionButton("proceed_2", "Proceed to Filter & Trim ->", class = "btn-proceed",
                     icon = icon("arrow-right"))
      ),

      div(class = "log-panel", id = "log_step2", htmlOutput("log2"))
    )),

    # ═══════════════════════════════════════════════════════════════════════
    # STEP 3: Filter & Trim
    # ═══════════════════════════════════════════════════════════════════════
    hidden(div(id = "panel_step3", class = "step-panel",

      div(class = "card",
        div(class = "card-header",
          div(class = "icon amber", icon("filter")),
          "Filter & Trim Parameters"
        ),
        div(class = "card-description",
          "Set truncation lengths based on quality profiles. Adjust additional parameters below."
        ),
        div(class = "grid-4",
          numericInput("truncLen_fwd", "Truncation Length (Fwd)", value = 240, min = 50, max = 500),
          numericInput("truncLen_rev", "Truncation Length (Rev)", value = 160, min = 50, max = 500),
          numericInput("maxEE_fwd", "Max Expected Errors (Fwd)", value = 2, min = 0, step = 0.5),
          numericInput("maxEE_rev", "Max Expected Errors (Rev)", value = 2, min = 0, step = 0.5)
        ),
        div(class = "grid-4", style = "margin-top: 12px;",
          numericInput("truncQ", "truncQ", value = 2, min = 0),
          numericInput("maxN", "maxN", value = 0, min = 0),
          div(style = "display:flex; align-items:flex-end; padding-bottom: 8px;",
            checkboxInput("rm_phix", "Remove PhiX", value = TRUE)
          ),
          div(style = "display:flex; align-items:flex-end; padding-bottom: 8px;",
            checkboxInput("compress_out", "Compress Output", value = TRUE)
          )
        ),
        div(style = "margin-top: 16px;",
          actionButton("btn_filter", "Run Filter & Trim", class = "btn-primary",
                       icon = icon("scissors"))
        ),
        uiOutput("progress_step3")
      ),

      div(class = "card",
        div(class = "card-header",
          div(class = "icon emerald", icon("table")),
          "Filtering Results"
        ),
        DTOutput("filter_table") %>% withSpinner(type = 6, color = "#3b82f6")
      ),

      div(class = "card",
        div(class = "card-header",
          div(class = "icon cyan", icon("chart-line")),
          "Post-Filter Quality Profiles"
        ),
        div(class = "card-description",
          "Inspect quality profiles of the filtered reads to confirm truncation was effective."
        ),
        div(class = "grid-2",
          numericInput("qp_filt_n_samples", "Number of samples to plot", value = 2, min = 1, max = 20),
          div(style = "display:flex; align-items:flex-end;",
            actionButton("btn_plot_filt_quality", "Generate Post-Filter Quality Plots", class = "btn-primary",
                         icon = icon("chart-area"))
          )
        ),
        uiOutput("qplot_filt_cards")
      ),

      div(class = "proceed-bar",
        actionButton("proceed_3", "Proceed to Dereplication ->", class = "btn-proceed",
                     icon = icon("arrow-right"))
      ),

      div(class = "log-panel", id = "log_step3", htmlOutput("log3"))
    )),

    # ═══════════════════════════════════════════════════════════════════════
    # STEP 4: Learn Errors & Dereplication
    # ═══════════════════════════════════════════════════════════════════════
    hidden(div(id = "panel_step4", class = "step-panel",

      div(class = "card",
        div(class = "card-header",
          div(class = "icon violet", icon("brain")),
          "Error Learning & Dereplication"
        ),
        div(class = "card-description",
          "DADA2 learns error rates from the data and then applies the core sample inference algorithm to dereplicate and denoise sequences."
        ),
        actionButton("btn_denoise", "Learn Errors & Dereplicate", class = "btn-primary",
                     icon = icon("microchip")),
        uiOutput("progress_step4")
      ),

      div(class = "card",
        div(class = "card-header",
          div(class = "icon blue", icon("chart-bar")),
          "Error Rate Plots"
        ),
        tabsetPanel(
          tabPanel("Forward Errors", plotOutput("errplot_fwd", height = "500px") %>% withSpinner(type = 6, color = "#3b82f6")),
          tabPanel("Reverse Errors", plotOutput("errplot_rev", height = "500px") %>% withSpinner(type = 6, color = "#3b82f6"))
        )
      ),

      div(class = "card",
        div(class = "card-header",
          div(class = "icon emerald", icon("info-circle")),
          "Dereplication Summary"
        ),
        DTOutput("denoise_summary") %>% withSpinner(type = 6, color = "#3b82f6")
      ),

      div(class = "proceed-bar",
        actionButton("proceed_4", "Proceed to Merge & Chimeras ->", class = "btn-proceed",
                     icon = icon("arrow-right"))
      ),

      div(class = "log-panel", id = "log_step4", htmlOutput("log4"))
    )),

    # ═══════════════════════════════════════════════════════════════════════
    # STEP 5: Merge & Chimera Removal
    # ═══════════════════════════════════════════════════════════════════════
    hidden(div(id = "panel_step5", class = "step-panel",

      div(class = "card",
        div(class = "card-header",
          div(class = "icon cyan", icon("link")),
          "Merge Paired Reads & Remove Chimeras"
        ),
        div(class = "card-description",
          "Merge denoised forward and reverse reads, construct the ASV table, and remove chimeric sequences."
        ),
        actionButton("btn_merge", "Merge & Remove Chimeras", class = "btn-primary",
                     icon = icon("code-merge")),
        uiOutput("progress_step5")
      ),

      div(class = "card",
        div(class = "card-header",
          div(class = "icon amber", icon("chart-pie")),
          "Sequence Table Summary"
        ),
        uiOutput("seqtab_stats_ui"),
        plotOutput("seqlen_plot", height = "300px") %>% withSpinner(type = 6, color = "#3b82f6")
      ),

      div(class = "card",
        div(class = "card-header",
          div(class = "icon emerald", icon("road")),
          "Read Tracking Through Pipeline"
        ),
        DTOutput("track_table") %>% withSpinner(type = 6, color = "#3b82f6"),
        plotOutput("track_plot", height = "400px") %>% withSpinner(type = 6, color = "#3b82f6"),
        div(style = "margin-top: 16px; display: flex; gap: 12px;",
          downloadButton("dl_track_step5", "Download Tracking Table (CSV)", class = "btn-download"),
          downloadButton("dl_track_plot", "Download Tracking Plot (PDF)", class = "btn-download")
        )
      ),

      div(class = "proceed-bar",
        actionButton("proceed_5", "Proceed to Taxonomy ->", class = "btn-proceed",
                     icon = icon("arrow-right"))
      ),

      div(class = "log-panel", id = "log_step5", htmlOutput("log5"))
    )),

    # ═══════════════════════════════════════════════════════════════════════
    # STEP 6: Taxonomy Assignment
    # ═══════════════════════════════════════════════════════════════════════
    hidden(div(id = "panel_step6", class = "step-panel",

      div(class = "card",
        div(class = "card-header",
          div(class = "icon rose", icon("dna")),
          "Taxonomy Assignment"
        ),
        div(class = "card-description",
          "Assign taxonomy using the SILVA reference database. The app will download the database if not already present."
        ),
        div(class = "grid-2",
          div(
            radioButtons("tax_method", "Taxonomy Method",
              choices = c("assignTaxonomy (Naive Bayesian)" = "bayesian",
                          "IdTaxa (DECIPHER)" = "idtaxa"),
              selected = "bayesian"
            )
          ),
          div(
            textInput("silva_dir", "Database Storage Directory",
                      value = file.path(path.expand("~"), "dada2_databases")),
            checkboxInput("add_species", "Add Species-Level Assignment", value = TRUE)
          )
        ),
        actionButton("btn_taxonomy", "Assign Taxonomy", class = "btn-primary",
                     icon = icon("tags")),
        uiOutput("progress_step6")
      ),

      div(class = "card",
        div(class = "card-header",
          div(class = "icon violet", icon("table")),
          "Taxonomy Table"
        ),
        DTOutput("taxa_table") %>% withSpinner(type = 6, color = "#3b82f6"),
        div(style = "margin-top: 16px; display: flex; gap: 12px;",
          downloadButton("dl_asv_step6", "Download ASV Table (CSV)", class = "btn-download"),
          downloadButton("dl_taxa_step6", "Download Taxonomy Table (CSV)", class = "btn-download")
        )
      ),

      div(class = "proceed-bar",
        actionButton("proceed_6", "Proceed to Phyloseq \u2192", class = "btn-proceed",
                     icon = icon("arrow-right"))
      ),

      div(class = "log-panel", id = "log_step6", htmlOutput("log6"))
    )),

    # ═══════════════════════════════════════════════════════════════════════
    # STEP 7: Phyloseq
    # ═══════════════════════════════════════════════════════════════════════
    hidden(div(id = "panel_step7", class = "step-panel",

      div(class = "card",
        div(class = "card-header",
          div(class = "icon emerald", icon("project-diagram")),
          "Phyloseq Construction"
        ),
        div(class = "card-description",
          "Build a phyloseq object for downstream analysis. Optionally upload sample metadata."
        ),
        div(style = "margin-bottom: 12px;",
          fileInput("metadata_file", "Upload Sample Metadata (CSV/TSV)",
                    accept = c(".csv", ".tsv", ".txt"))
        ),
        uiOutput("metadata_vartype_ui"),
        div(style = "margin-top: 16px;",
          actionButton("btn_phyloseq", "Build Phyloseq Object", class = "btn-primary",
                       icon = icon("cubes"))
        )
      ),

      div(class = "card",
        div(class = "card-header",
          div(class = "icon cyan", icon("balance-scale")),
          "Data Transformation"
        ),
        div(class = "card-description",
          "Choose how to normalize your data before visualization. Rarefaction subsamples all samples to an equal depth. Relative abundance converts counts to proportions."
        ),
        uiOutput("transform_ui"),
        div(style = "margin-top: 16px;",
          actionButton("btn_transform", "Apply Transformation", class = "btn-primary",
                       icon = icon("sync"))
        ),
        uiOutput("transform_status")
      ),

      div(class = "card",
        div(class = "card-header",
          div(class = "icon amber", icon("download")),
          "Export Results"
        ),
        div(class = "card-description",
          "Download all pipeline outputs."
        ),
        div(class = "grid-4",
          downloadButton("dl_asv", "ASV Table (CSV)", class = "btn-download"),
          downloadButton("dl_taxa", "Taxonomy (CSV)", class = "btn-download"),
          downloadButton("dl_track", "Tracking Table (CSV)", class = "btn-download"),
          downloadButton("dl_rdata", "Full Workspace (.RData)", class = "btn-download")
        ),
        div(class = "grid-2", style = "margin-top: 12px;",
          downloadButton("dl_phyloseq_rds", "Phyloseq Object (.rds)", class = "btn-download"),
          div()
        )
      ),

      div(class = "proceed-bar",
        actionButton("proceed_7", "Proceed to Figures \u2192", class = "btn-proceed",
                     icon = icon("arrow-right"))
      ),

      div(class = "log-panel", id = "log_step7", htmlOutput("log7"))
    )),

    # ═══════════════════════════════════════════════════════════════════════
    # STEP 8: Figures
    # ═══════════════════════════════════════════════════════════════════════
    hidden(div(id = "panel_step8", class = "step-panel",

      div(class = "card",
        div(class = "card-header",
          div(class = "icon blue", icon("chart-bar")),
          "Visualizations"
        ),
        tabsetPanel(
          tabPanel("Rarefaction Curves",
            div(class = "grid-2", style = "margin-top:12px; margin-bottom:12px;",
              selectInput("rare_color", "Color by", choices = NULL),
              numericInput("rare_step", "Step Size", value = 100, min = 10, max = 1000)
            ),
            plotOutput("rarefaction_plot", height = "500px") %>% withSpinner(type = 6, color = "#3b82f6"),
            div(style = "margin-top: 12px;",
              downloadButton("dl_rarefaction_png", "Download Rarefaction Curves (PNG)", class = "btn-download")
            )
          ),
          tabPanel("Alpha Diversity",
            div(style = "margin-top:12px; margin-bottom:12px; max-width: 300px;",
              selectInput("alpha_x", "Group by", choices = NULL)
            ),
            plotOutput("alpha_plot", height = "500px") %>% withSpinner(type = 6, color = "#3b82f6"),
            div(style = "margin-top: 12px;",
              downloadButton("dl_alpha_png", "Download Alpha Diversity Plot (PNG)", class = "btn-download")
            )
          ),
          tabPanel("Ordination (NMDS)",
            div(class = "grid-2", style = "margin-top:12px; margin-bottom:12px;",
              selectInput("ord_color", "Color by", choices = NULL),
              selectInput("ord_distance", "Distance Method",
                choices = c("bray", "jaccard", "unifrac", "wunifrac"), selected = "bray")
            ),
            plotOutput("ordination_plot", height = "500px") %>% withSpinner(type = 6, color = "#3b82f6"),
            div(style = "margin-top: 12px;",
              downloadButton("dl_ord_png", "Download NMDS Plot (PNG)", class = "btn-download")
            )
          ),
          tabPanel("Ordination (PCoA)",
            div(class = "grid-2", style = "margin-top:12px; margin-bottom:12px;",
              selectInput("pcoa_color", "Color by", choices = NULL),
              selectInput("pcoa_distance", "Distance Method",
                choices = c("bray", "jaccard", "unifrac", "wunifrac"), selected = "bray")
            ),
            plotOutput("pcoa_plot", height = "500px") %>% withSpinner(type = 6, color = "#3b82f6"),
            div(style = "margin-top: 12px;",
              downloadButton("dl_pcoa_png", "Download PCoA Plot (PNG)", class = "btn-download")
            )
          ),
          tabPanel("Abundance Plot",
            div(class = "grid-3", style = "margin-top:12px; margin-bottom:12px;",
              selectInput("bar_x", "X-axis variable", choices = NULL),
              selectInput("bar_fill", "Fill by Taxonomic Rank",
                choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
                selected = "Family"),
              numericInput("bar_top_n", "Top N Taxa", value = 10, min = 5, max = 100)
            ),
            plotOutput("bar_plot", height = "500px") %>% withSpinner(type = 6, color = "#3b82f6"),
            div(style = "margin-top: 12px;",
              downloadButton("dl_bar_png", "Download Abundance Plot (PNG)", class = "btn-download")
            )
          )
        )
      ),

      div(class = "proceed-bar",
        actionButton("proceed_8", "Proceed to PERMANOVA \u2192", class = "btn-proceed",
                     icon = icon("arrow-right"))
      ),

      div(class = "log-panel", id = "log_step8", htmlOutput("log8"))
    )),

    # ═══════════════════════════════════════════════════════════════════════
    # STEP 9: PERMANOVA
    # ═══════════════════════════════════════════════════════════════════════
    hidden(div(id = "panel_step9", class = "step-panel",

      div(class = "card",
        div(class = "card-header",
          div(class = "icon violet", icon("flask")),
          "PERMANOVA Analysis"
        ),
        div(class = "card-description",
          "Permutational Multivariate Analysis of Variance (PERMANOVA) tests whether community composition differs significantly between groups. Uses the vegan package."
        ),
        div(class = "grid-3",
          div(
            textInput("permanova_formula", "Formula",
                      placeholder = "e.g. Treatment + Age"),
            div(style = "font-size: 11px; color: var(--text-muted); margin-top: -8px; margin-bottom: 4px;",
              "Variables from your metadata separated by +")
          ),
          selectInput("permanova_distance", "Distance Method",
            choices = c("bray", "jaccard", "euclidean"), selected = "bray"),
          numericInput("permanova_perm", "Number of Permutations", value = 9999, min = 99, max = 99999, step = 100)
        ),
        div(class = "card-description", style = "font-size: 12px; color: var(--text-muted);",
          "Tip: For CLR-transformed data, Euclidean distance (Aitchison distance) is recommended."
        ),
        uiOutput("permanova_cmd_preview"),
        div(style = "margin-top: 12px;",
          actionButton("btn_permanova", "Run PERMANOVA", class = "btn-primary",
                       icon = icon("calculator"))
        )
      ),

      div(class = "card",
        div(class = "card-header",
          div(class = "icon blue", icon("table")),
          "PERMANOVA Results"
        ),
        div(class = "card-description",
          "Adonis2 test results showing whether community composition differs significantly between groups."
        ),
        uiOutput("permanova_results_ui")
      ),

      div(class = "card",
        div(class = "card-header",
          div(class = "icon amber", icon("balance-scale")),
          "Betadisper - Homogeneity of Dispersion"
        ),
        div(class = "card-description",
          "Tests whether group dispersions are homogeneous. A significant result here means differences detected by PERMANOVA may be due to differences in dispersion rather than location."
        ),
        uiOutput("betadisper_results_ui")
      ),

      uiOutput("pairwise_permanova_card"),

      div(class = "card",
        div(class = "card-header",
          div(class = "icon emerald", icon("download")),
          "Export PERMANOVA Results"
        ),
        div(style = "display: flex; gap: 12px;",
          downloadButton("dl_permanova_csv", "PERMANOVA Results (CSV)", class = "btn-download"),
          downloadButton("dl_betadisper_csv", "Betadisper Results (CSV)", class = "btn-download"),
          downloadButton("dl_pairwise_csv", "Pairwise PERMANOVA (CSV)", class = "btn-download")
        )
      ),

      div(class = "proceed-bar",
        actionButton("proceed_9", "Proceed to ANCOM-BC2 \u2192", class = "btn-proceed",
                     icon = icon("arrow-right"))
      ),

      div(class = "log-panel", id = "log_step9", htmlOutput("log9"))
    )),

    # =====================================================================
    # STEP 10: ANCOM-BC2
    # =====================================================================
    hidden(div(id = "panel_step10", class = "step-panel",

      div(class = "card",
        div(class = "card-header",
          div(class = "icon rose", icon("not-equal")),
          "ANCOM-BC2 - Differential Abundance Analysis"
        ),
        div(class = "card-description",
          "Analysis of Compositions of Microbiomes with Bias Correction 2 (ANCOM-BC2) identifies taxa whose absolute abundances differ significantly between groups. ANCOM-BC2 handles its own bias correction internally."
        ),
        div(class = "rerun-warning",
          icon("info-circle"),
          "ANCOM-BC2 uses the raw (untransformed) phyloseq object. Any transformation applied in Step 7 does not affect this analysis."
        ),
        div(class = "grid-2",
          div(
            textInput("ancom_fix_formula", "Fixed Effects Formula",
                      placeholder = "e.g. Treatment + Age + Sex"),
            div(style = "font-size: 11px; color: var(--text-muted); margin-top: -8px; margin-bottom: 12px;",
              "Variables from your metadata separated by +")
          ),
          div(
            textInput("ancom_rand_formula", "Random Effects Formula (optional)",
                      placeholder = "e.g. (1 | Subject)"),
            div(style = "font-size: 11px; color: var(--text-muted); margin-top: -8px; margin-bottom: 12px;",
              "lme4-style random effects for repeated measures")
          )
        ),
        div(class = "grid-3",
          selectInput("ancom_group", "Group Variable (for global/pairwise tests)", choices = NULL),
          selectInput("ancom_tax_level", "Taxonomic Level",
            choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
            selected = "Genus"),
          selectInput("ancom_p_adj", "P-value Adjustment Method",
            choices = c("holm", "BH", "bonferroni", "fdr", "none"),
            selected = "holm")
        ),
        div(class = "grid-3",
          numericInput("ancom_prv_cut", "Prevalence Filter", value = 0.10, min = 0, max = 1, step = 0.05),
          numericInput("ancom_lib_cut", "Library Size Cutoff", value = 1000, min = 0, step = 100),
          numericInput("ancom_alpha", "Significance Level (alpha)", value = 0.05, min = 0.001, max = 0.2, step = 0.01)
        ),

        # ── Advanced Settings (collapsed by default) ──
        tags$details(style = "margin-top: 16px; margin-bottom: 16px;",
          tags$summary(style = "cursor: pointer; font-weight: 500; color: var(--text-secondary); font-size: 13px;",
            icon("cog"), " Advanced Settings"
          ),
          div(style = "padding: 16px 0 0 0;",
            div(class = "grid-3",
              numericInput("ancom_pseudo", "Pseudo Count", value = 0, min = 0, step = 0.1),
              checkboxInput("ancom_pseudo_sens", "Sensitivity Analysis", value = TRUE),
              numericInput("ancom_s0_perc", "Regularization (s0_perc)", value = 0.05, min = 0, max = 1, step = 0.01)
            ),
            div(class = "grid-3",
              checkboxInput("ancom_struc_zero", "Structural Zero Detection", value = TRUE),
              checkboxInput("ancom_neg_lb", "Negative Lower Bound", value = TRUE),
              numericInput("ancom_n_cl", "Parallel Clusters", value = 1, min = 1, max = 16, step = 1)
            ),
            div(class = "grid-2",
              checkboxInput("ancom_global", "Global Test", value = TRUE),
              checkboxInput("ancom_pairwise", "Pairwise Directional Test", value = TRUE)
            ),
            div(class = "grid-2",
              checkboxInput("ancom_dunnet", "Dunnett's Type Test", value = FALSE),
              checkboxInput("ancom_trend", "Trend Test", value = FALSE)
            ),
            div(class = "card-description", style = "font-weight: 500; color: var(--text-primary); margin-top: 8px; margin-bottom: 8px;",
              "Iteration & EM Controls"),
            div(class = "grid-4",
              numericInput("ancom_iter_tol", "Iter Tolerance", value = 0.01, min = 1e-6, max = 1, step = 0.001),
              numericInput("ancom_iter_max", "Iter Max Steps", value = 20, min = 1, max = 200, step = 1),
              numericInput("ancom_em_tol", "EM Tolerance", value = 1e-5, min = 1e-8, max = 1, step = 1e-5),
              numericInput("ancom_em_max", "EM Max Steps", value = 100, min = 1, max = 500, step = 10)
            ),
            div(class = "card-description", style = "font-weight: 500; color: var(--text-primary); margin-top: 8px; margin-bottom: 8px;",
              "mdFDR Control"),
            div(class = "grid-2",
              selectInput("ancom_mdfdr_method", "FWER Control Method",
                choices = c("holm", "bonferroni", "BH", "hochberg", "hommel", "none"),
                selected = "holm"),
              numericInput("ancom_mdfdr_B", "Bootstrap Samples (B)", value = 100, min = 10, max = 1000, step = 10)
            )
          )
        ),

        uiOutput("ancom_cmd_preview"),
        div(style = "margin-top: 16px;",
          actionButton("btn_ancombc2", "Run ANCOM-BC2", class = "btn-primary",
                       icon = icon("play"))
        ),
        uiOutput("progress_step10")
      ),

      div(class = "card",
        div(class = "card-header",
          div(class = "icon blue", icon("globe")),
          "Global Test"
        ),
        div(class = "card-description",
          "Tests whether each taxon is differentially abundant across any of the groups."
        ),
        DTOutput("ancom_global_dt") %>% withSpinner(type = 6, color = "#3b82f6"),
        div(style = "margin-top: 12px;",
          downloadButton("dl_ancom_global", "Download Global Test (CSV)", class = "btn-download")
        )
      ),

      div(class = "card",
        div(class = "card-header",
          div(class = "icon cyan", icon("exchange-alt")),
          "Pairwise Directional Test"
        ),
        div(class = "card-description",
          "Identifies which taxa are differentially abundant between each pair of groups, with direction of change."
        ),
        DTOutput("ancom_pairwise_dt") %>% withSpinner(type = 6, color = "#3b82f6"),
        div(style = "margin-top: 12px;",
          downloadButton("dl_ancom_pairwise", "Download Pairwise Results (CSV)", class = "btn-download")
        )
      ),

      div(class = "card",
        div(class = "card-header",
          div(class = "icon amber", icon("chart-area")),
          "Volcano Plot"
        ),
        div(class = "card-description",
          "Log fold-change vs -log10(adjusted p-value). Taxa above the significance threshold are colored."
        ),
        uiOutput("ancom_volcano_selector"),
        plotOutput("ancom_volcano_plot", height = "500px") %>% withSpinner(type = 6, color = "#3b82f6"),
        div(style = "margin-top: 12px;",
          downloadButton("dl_ancom_volcano", "Download Volcano Plot (PNG)", class = "btn-download")
        )
      ),

      div(class = "card",
        div(class = "card-header",
          div(class = "icon violet", icon("th")),
          "Differential Abundance Heatmap"
        ),
        div(class = "card-description",
          "Heatmap of log fold-changes for significantly differentially abundant taxa across comparisons."
        ),
        plotOutput("ancom_heatmap", height = "600px") %>% withSpinner(type = 6, color = "#3b82f6"),
        div(style = "margin-top: 12px;",
          downloadButton("dl_ancom_heatmap", "Download Heatmap (PNG)", class = "btn-download")
        )
      ),

      div(class = "log-panel", id = "log_step10", htmlOutput("log10"))
    ))

  ) # end main-container
) # end fluidPage


# ══════════════════════════════════════════════════════════════════════════════
# SERVER
# ══════════════════════════════════════════════════════════════════════════════

server <- function(input, output, session) {

  # ── Reactive values ──────────────────────────────────────────────────────
  rv <- reactiveValues(
    current_step = 1,
    completed_steps = c(),
    startup_done = NULL,
    # File data
    fnFs = NULL, fnRs = NULL, sample_names = NULL,
    all_files = NULL,
    # Quality plots
    qplots_ready = FALSE,
    filt_qplots_ready = FALSE,
    # Background processes
    bg_filter = NULL, bg_filter_start = NULL,
    pending_filtFs = NULL, pending_filtRs = NULL,
    bg_denoise = NULL, bg_denoise_start = NULL, denoise_stage = NULL,
    bg_merge = NULL, bg_merge_start = NULL, merge_stage = NULL,
    bg_tax = NULL, bg_tax_start = NULL,
    # Filter
    filtFs = NULL, filtRs = NULL, filter_out = NULL,
    # Denoise
    errF = NULL, errR = NULL,
    dadaFs = NULL, dadaRs = NULL,
    # Merge
    mergers = NULL, seqtab = NULL, seqtab_nochim = NULL,
    # Taxonomy
    taxa = NULL,
    # Phyloseq
    ps = NULL, ps_transformed = NULL, transform_method = NULL, samdf = NULL,
    # PERMANOVA
    permanova_result = NULL, betadisper_result = NULL, pairwise_result = NULL,
    # ANCOM-BC2
    ancom_result = NULL, bg_ancom = NULL, bg_ancom_start = NULL,
    # Tracking
    track = NULL,
    # Logs
    log1 = "", log2 = "", log3 = "", log4 = "", log5 = "", log6 = "", log7 = "", log8 = "", log9 = "", log10 = "",
    # Step progress: list(pct, label, status)
    prog3 = list(pct = 0, label = "", status = "idle"),
    prog4 = list(pct = 0, label = "", status = "idle"),
    prog5 = list(pct = 0, label = "", status = "idle"),
    prog6 = list(pct = 0, label = "", status = "idle"),
    prog10 = list(pct = 0, label = "", status = "idle")
  )

  # ── Helper: update step progress ──
  set_progress <- function(step, pct, label, status = "running") {
    rv[[paste0("prog", step)]] <- list(pct = pct, label = label, status = status)
  }

  # ── Helper: render progress bar UI ──
  render_progress_bar <- function(prog) {
    if (prog$status == "idle") return(NULL)
    fill_class <- switch(prog$status,
      "done" = "step-progress-fill complete",
      "error" = "step-progress-fill",
      "running" = if (prog$pct == 0) "step-progress-fill indeterminate" else "step-progress-fill",
      "step-progress-fill"
    )
    status_icon <- switch(prog$status,
      "running" = "\u23f3",
      "done" = "\u2713",
      "error" = "\u2717",
      ""
    )
    pct_text <- if (prog$status == "running" && prog$pct > 0) paste0(prog$pct, "%") else
                if (prog$status == "done") "Done" else ""

    div(class = "step-progress-container",
      div(class = "step-progress-label",
        span(paste0(status_icon, "  ", prog$label)),
        span(class = "progress-status", pct_text)
      ),
      div(class = "step-progress-track",
        div(class = fill_class, style = paste0("width: ", prog$pct, "%;"))
      )
    )
  }

  output$progress_step3 <- renderUI({ render_progress_bar(rv$prog3) })
  output$progress_step4 <- renderUI({ render_progress_bar(rv$prog4) })
  output$progress_step5 <- renderUI({ render_progress_bar(rv$prog5) })
  output$progress_step6 <- renderUI({ render_progress_bar(rv$prog6) })
  output$progress_step10 <- renderUI({ render_progress_bar(rv$prog10) })

  # ── Session Save/Load ──────────────────────────────────────────────────

  SESSION_FILENAME <- "dada2_session.RData"

  get_session_path <- function() {
    path <- trimws(isolate(input$data_path))
    if (nzchar(path) && dir.exists(path)) {
      return(file.path(path, SESSION_FILENAME))
    }
    return(NULL)
  }

  auto_save_session <- function() {
    session_path <- get_session_path()
    if (is.null(session_path)) return()

    tryCatch({
      # Collect all saveable state (exclude background process handles)
      session_data <- list(
        data_path = isolate(input$data_path),
        completed_steps = rv$completed_steps,
        current_step = rv$current_step,
        fnFs = rv$fnFs, fnRs = rv$fnRs,
        sample_names = rv$sample_names,
        all_files = rv$all_files,
        filtFs = rv$filtFs, filtRs = rv$filtRs,
        filter_out = rv$filter_out,
        errF = rv$errF, errR = rv$errR,
        dadaFs = rv$dadaFs, dadaRs = rv$dadaRs,
        mergers = rv$mergers,
        seqtab = rv$seqtab, seqtab_nochim = rv$seqtab_nochim,
        taxa = rv$taxa,
        ps = rv$ps, ps_transformed = rv$ps_transformed,
        transform_method = rv$transform_method, samdf = rv$samdf,
        track = rv$track,
        log1 = rv$log1, log2 = rv$log2, log3 = rv$log3,
        log4 = rv$log4, log5 = rv$log5, log6 = rv$log6, log7 = rv$log7, log8 = rv$log8, log9 = rv$log9, log10 = rv$log10,
        # Save input parameters for reproducibility
        fwd_pattern = isolate(input$fwd_pattern),
        rev_pattern = isolate(input$rev_pattern),
        sample_delim = isolate(input$sample_delim),
        sample_element_from = isolate(input$sample_element_from),
        sample_element_to = isolate(input$sample_element_to),
        save_timestamp = Sys.time()
      )
      save(session_data, file = session_path)
      add_log(rv$current_step, paste("Session auto-saved to:", session_path), "success")
    }, error = function(e) {
      add_log(rv$current_step, paste("Session save failed:", e$message), "warn")
    })
  }

  restore_session <- function(session_path) {
    tryCatch({
      load(session_path)  # loads session_data

      # Restore reactive values
      rv$completed_steps <- session_data$completed_steps
      rv$fnFs <- session_data$fnFs
      rv$fnRs <- session_data$fnRs
      rv$sample_names <- session_data$sample_names

      # Update quality profile sample count max
      if (!is.null(session_data$fnFs)) {
        n_samples <- length(session_data$fnFs)
        updateNumericInput(session, "qp_n_samples", max = n_samples,
                           value = min(2, n_samples))
        updateNumericInput(session, "qp_filt_n_samples", max = n_samples,
                           value = min(2, n_samples))
      }
      rv$all_files <- session_data$all_files
      rv$filtFs <- session_data$filtFs
      rv$filtRs <- session_data$filtRs
      rv$filter_out <- session_data$filter_out
      rv$errF <- session_data$errF
      rv$errR <- session_data$errR
      rv$dadaFs <- session_data$dadaFs
      rv$dadaRs <- session_data$dadaRs
      rv$mergers <- session_data$mergers
      rv$seqtab <- session_data$seqtab
      rv$seqtab_nochim <- session_data$seqtab_nochim
      rv$taxa <- session_data$taxa
      rv$ps <- session_data$ps
      rv$ps_transformed <- session_data$ps_transformed
      rv$transform_method <- session_data$transform_method
      rv$samdf <- session_data$samdf
      rv$track <- session_data$track
      rv$log1 <- session_data$log1
      rv$log2 <- session_data$log2
      rv$log3 <- session_data$log3
      rv$log4 <- session_data$log4
      rv$log5 <- session_data$log5
      rv$log6 <- session_data$log6
      rv$log7 <- session_data$log7
      rv$log8 <- if (!is.null(session_data$log8)) session_data$log8 else ""
      rv$log9 <- if (!is.null(session_data$log9)) session_data$log9 else ""
      rv$log10 <- if (!is.null(session_data$log10)) session_data$log10 else ""

      # Restore UI inputs
      updateTextInput(session, "data_path", value = session_data$data_path)
      if (!is.null(session_data$fwd_pattern))
        updateTextInput(session, "fwd_pattern", value = session_data$fwd_pattern)
      if (!is.null(session_data$rev_pattern))
        updateTextInput(session, "rev_pattern", value = session_data$rev_pattern)
      if (!is.null(session_data$sample_delim))
        updateTextInput(session, "sample_delim", value = session_data$sample_delim)
      if (!is.null(session_data$sample_element_from))
        updateNumericInput(session, "sample_element_from", value = session_data$sample_element_from)
      if (!is.null(session_data$sample_element_to))
        updateNumericInput(session, "sample_element_to", value = session_data$sample_element_to)

      # Set quality plot flags if those steps were completed
      if (2 %in% rv$completed_steps) rv$qplots_ready <- TRUE
      if (3 %in% rv$completed_steps) rv$filt_qplots_ready <- TRUE

      # Navigate to first incomplete step
      all_steps <- 1:10
      next_step <- min(setdiff(all_steps, rv$completed_steps), 10)
      rv$current_step <- next_step

      # Validate file paths
      data_path <- session_data$data_path
      paths_valid <- TRUE
      if (!is.null(data_path) && !dir.exists(data_path)) {
        paths_valid <- FALSE
      }
      if (!is.null(rv$fnFs) && length(rv$fnFs) > 0 && !file.exists(rv$fnFs[1])) {
        paths_valid <- FALSE
      }

      if (!paths_valid) {
        showModal(modalDialog(
          title = "File Paths Changed",
          div(class = "rerun-warning",
            icon("exclamation-triangle"),
            "The original FASTQ directory or files could not be found. Please update the path in Step 1."
          ),
          div(style = "color: var(--text-secondary); font-size: 13px;",
            p(paste("Original path:", data_path)),
            p("Your pipeline state has been restored, but you may need to update the directory path and re-scan files if the data has moved.")
          ),
          footer = tagList(
            actionButton("modal_path_ok", "Go to Step 1", class = "btn-primary")
          ),
          easyClose = FALSE
        ))
      }

      add_log(next_step, paste("Session restored. Resuming from Step", next_step, "."), "success")
      return(TRUE)
    }, error = function(e) {
      showNotification(paste("Failed to restore session:", e$message), type = "error")
      return(FALSE)
    })
  }

  # Handle path update modal dismiss
  observeEvent(input$modal_path_ok, {
    removeModal()
    rv$current_step <- 1
  })

  # ── Startup: check for existing session ──
  observe({
    # Run once on startup
    if (!is.null(rv$startup_done)) return()
    rv$startup_done <- TRUE

    # Check common locations for session files
    # We can't check data_path yet (empty on startup), so check home dir
    # The user will see the startup modal
    showModal(modalDialog(
      title = "DADA2 Pipeline",
      div(style = "text-align: center; margin-bottom: 16px;",
        div(class = "app-logo", style = "margin: 0 auto 12px; width: 56px; height: 56px; font-size: 20px;", "16S"),
        div(style = "font-size: 16px; font-weight: 500; color: var(--text-primary);", "Welcome to the DADA2 Analysis Pipeline")
      ),
      div(style = "color: var(--text-secondary); font-size: 14px; margin-bottom: 16px;",
        p("Start a new analysis or load a previous session."),
        p("To resume a previous session, enter the path to your FASTQ directory below. If a saved session exists there, it will be restored.")
      ),
      textInput("resume_path", "FASTQ Directory (to check for saved session)",
                placeholder = "/path/to/your/fastq/files"),
      uiOutput("resume_session_info"),
      footer = tagList(
        actionButton("btn_new_analysis", "New Analysis", class = "btn-secondary"),
        actionButton("btn_resume_session", "Resume Session", class = "btn-primary")
      ),
      easyClose = FALSE
    ))
  })

  # Show session info when path is entered
  output$resume_session_info <- renderUI({
    req(input$resume_path)
    path <- trimws(input$resume_path)
    session_file <- file.path(path, SESSION_FILENAME)

    if (file.exists(session_file)) {
      tryCatch({
        load(session_file)
        step_names <- c("Setup & Files", "Quality Profiles", "Filter & Trim",
                        "Dereplication", "Merge & Chimeras", "Taxonomy", "Phyloseq", "Figures", "PERMANOVA", "ANCOM-BC2")
        completed <- session_data$completed_steps
        last_step <- if (length(completed) > 0) max(completed) else 0
        last_step_name <- if (last_step > 0) step_names[last_step] else "None"
        timestamp <- format(session_data$save_timestamp, "%Y-%m-%d %H:%M:%S")
        n_samples <- length(session_data$sample_names)

        div(class = "session-info-card",
          div(class = "info-label", "Saved Session Found"),
          div(class = "info-value", style = "margin-top: 6px;",
            tags$strong("Last completed step: "), last_step_name, tags$br(),
            tags$strong("Samples: "), n_samples, tags$br(),
            tags$strong("Saved: "), timestamp, tags$br(),
            tags$strong("Steps completed: "), paste(completed, collapse = ", ")
          )
        )
      }, error = function(e) {
        div(style = "color: var(--accent-rose); font-size: 13px; margin-top: 8px;",
          paste("Found session file but could not read it:", e$message))
      })
    } else if (dir.exists(path)) {
      div(style = "color: var(--text-muted); font-size: 13px; margin-top: 8px;",
        "No saved session found in this directory.")
    } else {
      div(style = "color: var(--accent-rose); font-size: 13px; margin-top: 8px;",
        "Directory does not exist.")
    }
  })

  # New analysis button
  observeEvent(input$btn_new_analysis, {
    removeModal()
  })

  # Resume session button
  observeEvent(input$btn_resume_session, {
    req(input$resume_path)
    path <- trimws(input$resume_path)
    session_file <- file.path(path, SESSION_FILENAME)

    if (file.exists(session_file)) {
      removeModal()
      restore_session(session_file)
    } else {
      showNotification("No saved session found at that path.", type = "warning")
    }
  })

  # ── Re-run warning when going back to completed steps ──
  observeEvent(input$nav_step, {
    step <- input$nav_step
    if (step %in% rv$completed_steps && step < max(rv$completed_steps, 0)) {
      # User is going back to a completed step - warn about invalidation
      later_steps <- rv$completed_steps[rv$completed_steps > step]
      if (length(later_steps) > 0) {
        step_names <- c("Setup & Files", "Quality Profiles", "Filter & Trim",
                        "Dereplication", "Merge & Chimeras", "Taxonomy", "Phyloseq", "Figures", "PERMANOVA", "ANCOM-BC2")
        invalidated <- paste(step_names[later_steps], collapse = ", ")
        showNotification(
          paste0("Re-running this step will invalidate subsequent completed steps: ", invalidated, "."),
          type = "warning", duration = 8
        )
      }
    }
    rv$current_step <- step
  })

  # ── Helper: append log ──
  add_log <- function(step, msg, type = "info") {
    icon_map <- c(info = "&#9432;", success = "&#10003;", error = "&#10007;", warn = "&#9888;")
    class_map <- c(info = "log-info", success = "log-success", error = "log-error", warn = "log-warn")
    timestamp <- format(Sys.time(), "%H:%M:%S")
    html <- sprintf('<div class="%s">[%s] %s %s</div>',
                    class_map[type], timestamp, icon_map[type], msg)
    log_name <- paste0("log", step)
    rv[[log_name]] <- paste0(rv[[log_name]], html)
  }

  # ── Render logs ──
  output$log1 <- renderUI(HTML(rv$log1))
  output$log2 <- renderUI(HTML(rv$log2))
  output$log3 <- renderUI(HTML(rv$log3))
  output$log4 <- renderUI(HTML(rv$log4))
  output$log5 <- renderUI(HTML(rv$log5))
  output$log6 <- renderUI(HTML(rv$log6))
  output$log7 <- renderUI(HTML(rv$log7))
  output$log8 <- renderUI(HTML(rv$log8))
  output$log9 <- renderUI(HTML(rv$log9))
  output$log10 <- renderUI(HTML(rv$log10))

  # ── Navigation ──────────────────────────────────────────────────────────
  observe({
    step <- rv$current_step
    completed <- rv$completed_steps
    for (i in 1:10) {
      toggleClass(id = paste0("step_nav_", i), class = "active", condition = (i == step))
      toggleClass(id = paste0("step_nav_", i), class = "completed", condition = (i %in% completed & i != step))
      toggle(id = paste0("panel_step", i), condition = (i == step))
    }
    # Enable/disable proceed buttons based on step completion
    for (i in 1:7) {
      toggleState(id = paste0("proceed_", i), condition = (i %in% completed))
    }
    # Proceed to PERMANOVA and ANCOM-BC2 are always enabled
    shinyjs::enable("proceed_8")
    shinyjs::enable("proceed_9")
  })

  # Proceed button click handlers
  observeEvent(input$proceed_1, { rv$current_step <- 2 })
  observeEvent(input$proceed_2, { rv$current_step <- 3 })
  observeEvent(input$proceed_3, { rv$current_step <- 4 })
  observeEvent(input$proceed_4, { rv$current_step <- 5 })
  observeEvent(input$proceed_5, { rv$current_step <- 6 })
  observeEvent(input$proceed_6, { rv$current_step <- 7 })
  observeEvent(input$proceed_7, { rv$current_step <- 8 })
  observeEvent(input$proceed_8, {
    rv$completed_steps <- union(rv$completed_steps, 8)
    rv$current_step <- 9
  })
  observeEvent(input$proceed_9, { rv$current_step <- 10 })

  # ═══════════════════════════════════════════════════════════════════════
  # STEP 1: Scan files & extract samples
  # ═══════════════════════════════════════════════════════════════════════

  observeEvent(input$btn_scan_files, {
    req(input$data_path)
    path <- trimws(input$data_path)
    add_log(1, paste("Scanning directory:", path))

    if (!dir.exists(path)) {
      add_log(1, "Directory does not exist!", "error")
      return()
    }

    all_files <- list.files(path, full.names = FALSE)
    fastq_files <- all_files[grepl("\\.(fastq|fastq\\.gz|fq|fq\\.gz)$", all_files, ignore.case = TRUE)]

    if (length(fastq_files) == 0) {
      add_log(1, "No FASTQ files found in directory.", "error")
      return()
    }

    rv$all_files <- fastq_files
    add_log(1, paste("Found", length(fastq_files), "FASTQ files."), "success")

    # Auto-detect patterns
    detected_fwd <- NULL
    detected_rev <- NULL
    for (p in common_fwd_patterns) {
      matches <- fastq_files[grepl(gsub("\\.", "\\\\.", p), fastq_files)]
      if (length(matches) > 0) { detected_fwd <- p; break }
    }
    for (p in common_rev_patterns) {
      matches <- fastq_files[grepl(gsub("\\.", "\\\\.", p), fastq_files)]
      if (length(matches) > 0) { detected_rev <- p; break }
    }

    if (!is.null(detected_fwd) && !is.null(detected_rev)) {
      updateTextInput(session, "fwd_pattern", value = detected_fwd)
      updateTextInput(session, "rev_pattern", value = detected_rev)
      add_log(1, paste("Auto-detected patterns - Fwd:", detected_fwd, "  Rev:", detected_rev), "success")

      # Load files with detected patterns
      fnFs <- sort(list.files(path, pattern = gsub("\\.", "\\\\.", detected_fwd), full.names = TRUE))
      fnRs <- sort(list.files(path, pattern = gsub("\\.", "\\\\.", detected_rev), full.names = TRUE))
      rv$fnFs <- fnFs
      rv$fnRs <- fnRs
      add_log(1, paste("Forward files:", length(fnFs), " | Reverse files:", length(fnRs)), "info")
    } else {
      add_log(1, "Could not auto-detect patterns. Please enter them manually and click 'Extract Sample Names'.", "warn")
    }
  })

  output$detected_pattern_ui <- renderUI({
    req(rv$all_files)
    div(class = "card-description", style = "margin-top:8px;",
      paste("Detected file extensions:", paste(unique(tools::file_ext(rv$all_files)), collapse = ", "))
    )
  })

  observeEvent(input$btn_extract_samples, {
    path <- trimws(input$data_path)
    req(path, input$fwd_pattern, input$rev_pattern)

    fwd_pat <- gsub("\\.", "\\\\.", input$fwd_pattern)
    rev_pat <- gsub("\\.", "\\\\.", input$rev_pattern)

    fnFs <- sort(list.files(path, pattern = fwd_pat, full.names = TRUE))
    fnRs <- sort(list.files(path, pattern = rev_pat, full.names = TRUE))

    if (length(fnFs) == 0 || length(fnRs) == 0) {
      add_log(1, "No files matched the specified patterns.", "error")
      return()
    }

    if (length(fnFs) != length(fnRs)) {
      add_log(1, paste("Warning: Unequal file counts - Fwd:", length(fnFs), "Rev:", length(fnRs)), "warn")
    }

    rv$fnFs <- fnFs
    rv$fnRs <- fnRs

    delim <- input$sample_delim
    elem_from <- input$sample_element_from
    elem_to <- input$sample_element_to

    # Extract sample names using range of indices, joined by delimiter
    sample_names <- sapply(strsplit(basename(fnFs), delim, fixed = TRUE), function(parts) {
      idx <- seq(elem_from, min(elem_to, length(parts)))
      paste(parts[idx], collapse = delim)
    })
    rv$sample_names <- sample_names

    # Update quality profile sample count max to total number of samples
    n_samples <- length(fnFs)
    updateNumericInput(session, "qp_n_samples", max = n_samples,
                       value = min(2, n_samples))
    updateNumericInput(session, "qp_filt_n_samples", max = n_samples,
                       value = min(2, n_samples))

    add_log(1, paste("Extracted", length(sample_names), "sample names using elements",
                     elem_from, "to", elem_to, "(delimiter: '", delim, "')."), "success")
    rv$completed_steps <- union(rv$completed_steps, 1)
      auto_save_session()
  })

  # Live preview of sample name extraction
  output$sample_name_preview <- renderUI({
    req(rv$fnFs, input$sample_delim)
    example_file <- basename(rv$fnFs[1])
    delim <- input$sample_delim
    elem_from <- input$sample_element_from
    elem_to <- input$sample_element_to
    parts <- strsplit(example_file, delim, fixed = TRUE)[[1]]
    idx <- seq(elem_from, min(elem_to, length(parts)))
    preview <- paste(parts[idx], collapse = delim)
    tags$span(
      style = "font-family: 'JetBrains Mono', monospace; font-size: 13px;",
      tags$span(style = "color: var(--text-muted);", paste0("Preview: \"", example_file, "\" -> ")),
      tags$span(style = "color: var(--accent-cyan); font-weight: 600;", paste0("\"", preview, "\""))
    )
  })

  output$file_summary_ui <- renderUI({
    req(rv$fnFs)
    n <- length(rv$fnFs)
    div(class = "grid-3", style = "margin-bottom: 16px;",
      div(class = "stat-card",
        div(class = "stat-value", length(rv$fnFs)),
        div(class = "stat-label", "Forward Files")
      ),
      div(class = "stat-card",
        div(class = "stat-value", length(rv$fnRs)),
        div(class = "stat-label", "Reverse Files")
      ),
      div(class = "stat-card",
        div(class = "stat-value", ifelse(is.null(rv$sample_names), "-", length(rv$sample_names))),
        div(class = "stat-label", "Samples")
      )
    )
  })

  output$sample_table <- renderDT({
    req(rv$sample_names, rv$fnFs, rv$fnRs)
    df <- data.frame(
      Sample = rv$sample_names,
      Forward = basename(rv$fnFs),
      Reverse = basename(rv$fnRs),
      stringsAsFactors = FALSE
    )
    datatable(df, options = list(pageLength = 10, scrollX = TRUE, dom = 'frtip'),
              rownames = FALSE, class = 'compact')
  })

  # ═══════════════════════════════════════════════════════════════════════
  # STEP 2: Quality Profiles
  # ═══════════════════════════════════════════════════════════════════════

  observeEvent(input$btn_plot_quality, {
    req(rv$fnFs, rv$fnRs)
    n <- min(input$qp_n_samples, length(rv$fnFs))
    add_log(2, paste("Generating quality plots for", n, "samples..."))

    rv$qplots_ready <- TRUE

    output$qplot_fwd <- renderPlot({
      plotQualityProfile(rv$fnFs[1:n]) + ggtitle("Forward Reads Quality Profile")
    })

    output$qplot_rev <- renderPlot({
      plotQualityProfile(rv$fnRs[1:n]) + ggtitle("Reverse Reads Quality Profile")
    })

    add_log(2, "Quality plots generated. Use these to set truncation lengths in Step 3.", "success")
    rv$completed_steps <- union(rv$completed_steps, 2)
      auto_save_session()
  })

  # Only show quality plot cards after button is pressed
  output$qplot_cards <- renderUI({
    req(rv$qplots_ready)
    tagList(
      div(class = "card",
        div(class = "card-header",
          div(class = "icon blue", icon("arrow-right")),
          "Forward Reads Quality"
        ),
        plotOutput("qplot_fwd", height = "400px") %>% withSpinner(type = 6, color = "#3b82f6")
      ),
      div(class = "card",
        div(class = "card-header",
          div(class = "icon rose", icon("arrow-left")),
          "Reverse Reads Quality"
        ),
        plotOutput("qplot_rev", height = "400px") %>% withSpinner(type = 6, color = "#3b82f6")
      ),
      div(style = "display: flex; gap: 12px; margin-top: 4px;",
        downloadButton("dl_raw_qplot_fwd", "Download Forward Quality Plot (PNG)", class = "btn-download"),
        downloadButton("dl_raw_qplot_rev", "Download Reverse Quality Plot (PNG)", class = "btn-download")
      )
    )
  })

  # ═══════════════════════════════════════════════════════════════════════
  # STEP 3: Filter & Trim (async with r_bg for instant UI feedback)
  # ═══════════════════════════════════════════════════════════════════════

  observeEvent(input$btn_filter, {
    req(rv$fnFs, rv$fnRs, rv$sample_names)
    path <- trimws(input$data_path)
    add_log(3, "Starting filter & trim...")

    filtFs <- file.path(path, "filtered", paste0(rv$sample_names, "_F_filt.fastq.gz"))
    filtRs <- file.path(path, "filtered", paste0(rv$sample_names, "_R_filt.fastq.gz"))
    names(filtFs) <- rv$sample_names
    names(filtRs) <- rv$sample_names

    # Store for use after async completes
    rv$pending_filtFs <- filtFs
    rv$pending_filtRs <- filtRs

    # Show running indicator instantly
    set_progress(3, 0, "Launching filterAndTrim (multithreaded)...", "running")
    shinyjs::disable("btn_filter")

    # Launch background process
    rv$bg_filter <- callr::r_bg(
      function(fnFs, filtFs, fnRs, filtRs, truncLen, maxN, maxEE, truncQ, rm_phix, compress, mt) {
        library(dada2)
        filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                      truncLen = truncLen, maxN = maxN, maxEE = maxEE,
                      truncQ = truncQ, rm.phix = rm_phix,
                      compress = compress, multithread = mt)
      },
      args = list(
        fnFs = rv$fnFs, filtFs = filtFs, fnRs = rv$fnRs, filtRs = filtRs,
        truncLen = c(input$truncLen_fwd, input$truncLen_rev),
        maxN = input$maxN,
        maxEE = c(input$maxEE_fwd, input$maxEE_rev),
        truncQ = input$truncQ,
        rm_phix = input$rm_phix,
        compress = input$compress_out,
        mt = can_multithread
      ),
      supervise = TRUE
    )
    rv$bg_filter_start <- Sys.time()
  })

  # Poll for filter completion
  observe({
    req(rv$bg_filter)
    if (rv$bg_filter$is_alive()) {
      elapsed <- as.numeric(difftime(Sys.time(), rv$bg_filter_start, units = "secs"))
      elapsed_txt <- if (elapsed > 60) paste0(round(elapsed/60, 1), " min") else paste0(round(elapsed), "s")
      set_progress(3, 0, paste0("Running filterAndTrim... (", elapsed_txt, " elapsed)"), "running")
      invalidateLater(1000)
    } else {
      tryCatch({
        out <- rv$bg_filter$get_result()
        filtFs <- rv$pending_filtFs
        filtRs <- rv$pending_filtRs

        rv$filtFs <- filtFs[file.exists(filtFs)]
        rv$filtRs <- filtRs[file.exists(filtRs)]
        rv$filter_out <- out

        existing <- file.exists(filtFs)
        if (sum(existing) < length(rv$sample_names)) {
          dropped <- rv$sample_names[!existing]
          add_log(3, paste("Dropped", length(dropped), "samples with 0 reads after filtering:", paste(dropped, collapse = ", ")), "warn")
          rv$sample_names <- rv$sample_names[existing]
          rv$fnFs <- rv$fnFs[existing]
          rv$fnRs <- rv$fnRs[existing]
          rv$filtFs <- filtFs[existing]
          rv$filtRs <- filtRs[existing]
        }

        total_in <- sum(out[, "reads.in"])
        total_out <- sum(out[, "reads.out"])
        pct <- round(total_out / total_in * 100, 1)
        add_log(3, paste("Filtering complete.", total_out, "/", total_in, "reads passed (", pct, "%)."), "success")
        set_progress(3, 100, "Filter & trim complete", "done")
        rv$completed_steps <- union(rv$completed_steps, 3)
      auto_save_session()
      }, error = function(e) {
        add_log(3, paste("Error:", e$message), "error")
        set_progress(3, 0, paste("Error:", e$message), "error")
      })
      shinyjs::enable("btn_filter")
      rv$bg_filter <- NULL
    }
  })

  output$filter_table <- renderDT({
    req(rv$filter_out)
    df <- as.data.frame(rv$filter_out)
    df$Forward_File <- rownames(df)
    # Match reverse file names by sample
    rev_basenames <- basename(rv$fnRs)
    df$Reverse_File <- rev_basenames[seq_len(nrow(df))]
    df$Sample <- rv$sample_names[seq_len(nrow(df))]
    df$Retention <- paste0(round(df$reads.out / df$reads.in * 100, 1), "%")
    df <- df[, c("Sample", "Forward_File", "Reverse_File", "reads.in", "reads.out", "Retention")]
    datatable(df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        scrollY = "400px",
        scroller = TRUE,
        dom = 'frtip'
      ),
      rownames = FALSE, class = 'compact'
    )
  })

  # ── Post-filter quality plots ──
  observeEvent(input$btn_plot_filt_quality, {
    req(rv$filtFs, rv$filtRs)
    n <- min(input$qp_filt_n_samples, length(rv$filtFs))
    add_log(3, paste("Generating post-filter quality plots for", n, "samples..."))
    rv$filt_qplots_ready <- TRUE

    output$qplot_filt_fwd <- renderPlot({
      plotQualityProfile(rv$filtFs[1:n]) + ggtitle("Forward Reads Quality (Post-Filter)")
    })
    output$qplot_filt_rev <- renderPlot({
      plotQualityProfile(rv$filtRs[1:n]) + ggtitle("Reverse Reads Quality (Post-Filter)")
    })
    add_log(3, "Post-filter quality plots generated.", "success")
  })

  output$qplot_filt_cards <- renderUI({
    req(rv$filt_qplots_ready)
    tagList(
      div(style = "margin-top: 16px;",
        div(class = "card-header", style = "margin-bottom: 8px;",
          div(class = "icon blue", icon("arrow-right")),
          "Forward Reads (Filtered)"
        ),
        plotOutput("qplot_filt_fwd", height = "400px") %>% withSpinner(type = 6, color = "#3b82f6")
      ),
      div(style = "margin-top: 16px;",
        div(class = "card-header", style = "margin-bottom: 8px;",
          div(class = "icon rose", icon("arrow-left")),
          "Reverse Reads (Filtered)"
        ),
        plotOutput("qplot_filt_rev", height = "400px") %>% withSpinner(type = 6, color = "#3b82f6")
      ),
      div(style = "display: flex; gap: 12px; margin-top: 12px;",
        downloadButton("dl_filt_qplot_fwd", "Download Fwd Filtered Quality (PNG)", class = "btn-download"),
        downloadButton("dl_filt_qplot_rev", "Download Rev Filtered Quality (PNG)", class = "btn-download")
      )
    )
  })

  # ── Quality plot download handlers (raw - Step 2) ──
  output$dl_raw_qplot_fwd <- downloadHandler(
    filename = function() paste0("Raw_Fwd_QualityProfile_", Sys.Date(), ".png"),
    content = function(file) {
      req(rv$fnFs)
      n <- min(input$qp_n_samples, length(rv$fnFs))
      p <- plotQualityProfile(rv$fnFs[1:n]) + ggtitle("Forward Reads Quality Profile")
      ggsave(file, plot = p, width = 10, height = 6, dpi = 300, bg = "white")
    }
  )
  output$dl_raw_qplot_rev <- downloadHandler(
    filename = function() paste0("Raw_Rev_QualityProfile_", Sys.Date(), ".png"),
    content = function(file) {
      req(rv$fnRs)
      n <- min(input$qp_n_samples, length(rv$fnRs))
      p <- plotQualityProfile(rv$fnRs[1:n]) + ggtitle("Reverse Reads Quality Profile")
      ggsave(file, plot = p, width = 10, height = 6, dpi = 300, bg = "white")
    }
  )

  # ── Quality plot download handlers (filtered - Step 3) ──
  output$dl_filt_qplot_fwd <- downloadHandler(
    filename = function() paste0("Filtered_Fwd_QualityProfile_", Sys.Date(), ".png"),
    content = function(file) {
      req(rv$filtFs)
      n <- min(input$qp_filt_n_samples, length(rv$filtFs))
      p <- plotQualityProfile(rv$filtFs[1:n]) + ggtitle("Forward Reads Quality (Post-Filter)")
      ggsave(file, plot = p, width = 10, height = 6, dpi = 300, bg = "white")
    }
  )
  output$dl_filt_qplot_rev <- downloadHandler(
    filename = function() paste0("Filtered_Rev_QualityProfile_", Sys.Date(), ".png"),
    content = function(file) {
      req(rv$filtRs)
      n <- min(input$qp_filt_n_samples, length(rv$filtRs))
      p <- plotQualityProfile(rv$filtRs[1:n]) + ggtitle("Reverse Reads Quality (Post-Filter)")
      ggsave(file, plot = p, width = 10, height = 6, dpi = 300, bg = "white")
    }
  )

  # ═══════════════════════════════════════════════════════════════════════
  # STEP 4: Learn Errors & Dereplication (async state machine)
  # ═══════════════════════════════════════════════════════════════════════

  observeEvent(input$btn_denoise, {
    req(rv$filtFs, rv$filtRs)
    add_log(4, "Starting error learning & dereplication pipeline...")
    set_progress(4, 0, "Launching forward error learning...", "running")
    shinyjs::disable("btn_denoise")
    rv$denoise_stage <- 1  # 1=errF, 2=errR, 3=dadaF, 4=dadaR
    rv$bg_denoise_start <- Sys.time()

    # Stage 1: learn forward errors
    rv$bg_denoise <- callr::r_bg(
      function(filtFs, mt) { library(dada2); learnErrors(filtFs, multithread = mt) },
      args = list(filtFs = rv$filtFs, mt = can_multithread),
      supervise = TRUE
    )
  })

  # Poll dereplication stages
  observe({
    req(rv$bg_denoise)
    if (rv$bg_denoise$is_alive()) {
      elapsed <- as.numeric(difftime(Sys.time(), rv$bg_denoise_start, units = "secs"))
      elapsed_txt <- if (elapsed > 60) paste0(round(elapsed/60, 1), " min") else paste0(round(elapsed), "s")
      stage_labels <- c("Learning forward error model...", "Learning reverse error model...",
                        "Dereplicating forward reads...", "Dereplicating reverse reads...")
      stage_pcts <- c(5, 30, 55, 80)
      set_progress(4, stage_pcts[rv$denoise_stage],
                   paste0(stage_labels[rv$denoise_stage], " (", elapsed_txt, " elapsed)"), "running")
      invalidateLater(1000)
    } else {
      tryCatch({
        result <- rv$bg_denoise$get_result()
        stage <- rv$denoise_stage

        if (stage == 1) {
          rv$errF <- result
          add_log(4, "Forward error model learned.", "success")
          rv$denoise_stage <- 2
          rv$bg_denoise_start <- Sys.time()
          rv$bg_denoise <- callr::r_bg(
            function(filtRs, mt) { library(dada2); learnErrors(filtRs, multithread = mt) },
            args = list(filtRs = rv$filtRs, mt = can_multithread), supervise = TRUE
          )
        } else if (stage == 2) {
          rv$errR <- result
          add_log(4, "Reverse error model learned.", "success")
          rv$denoise_stage <- 3
          rv$bg_denoise_start <- Sys.time()
          rv$bg_denoise <- callr::r_bg(
            function(filtFs, errF, mt) { library(dada2); dada(filtFs, err = errF, multithread = mt) },
            args = list(filtFs = rv$filtFs, errF = rv$errF, mt = can_multithread), supervise = TRUE
          )
        } else if (stage == 3) {
          rv$dadaFs <- result
          add_log(4, "Forward reads dereplicated.", "success")
          rv$denoise_stage <- 4
          rv$bg_denoise_start <- Sys.time()
          rv$bg_denoise <- callr::r_bg(
            function(filtRs, errR, mt) { library(dada2); dada(filtRs, err = errR, multithread = mt) },
            args = list(filtRs = rv$filtRs, errR = rv$errR, mt = can_multithread), supervise = TRUE
          )
        } else if (stage == 4) {
          rv$dadaRs <- result
          add_log(4, "Reverse reads dereplicated.", "success")
          set_progress(4, 100, "Error learning & dereplication complete", "done")
          rv$completed_steps <- union(rv$completed_steps, 4)
      auto_save_session()
          shinyjs::enable("btn_denoise")
          rv$bg_denoise <- NULL
        }
      }, error = function(e) {
        add_log(4, paste("Error:", e$message), "error")
        set_progress(4, 0, paste("Error:", e$message), "error")
        shinyjs::enable("btn_denoise")
        rv$bg_denoise <- NULL
      })
    }
  })

  output$errplot_fwd <- renderPlot({
    req(rv$errF)
    plotErrors(rv$errF, nominalQ = TRUE) + ggtitle("Forward Error Rates")
  })

  output$errplot_rev <- renderPlot({
    req(rv$errR)
    plotErrors(rv$errR, nominalQ = TRUE) + ggtitle("Reverse Error Rates")
  })

  output$denoise_summary <- renderDT({
    req(rv$dadaFs, rv$dadaRs, rv$sample_names)
    df <- data.frame(
      Sample = rv$sample_names,
      Fwd_Denoised = sapply(rv$dadaFs, function(x) sum(getUniques(x))),
      Fwd_ASVs = sapply(rv$dadaFs, function(x) length(getUniques(x))),
      Rev_Denoised = sapply(rv$dadaRs, function(x) sum(getUniques(x))),
      Rev_ASVs = sapply(rv$dadaRs, function(x) length(getUniques(x))),
      stringsAsFactors = FALSE
    )
    datatable(df, options = list(pageLength = 15, scrollX = TRUE, dom = 'frtip'),
              rownames = FALSE, class = 'compact')
  })

  # ═══════════════════════════════════════════════════════════════════════
  # STEP 5: Merge & Chimera Removal
  # ═══════════════════════════════════════════════════════════════════════

  observeEvent(input$btn_merge, {
    req(rv$dadaFs, rv$dadaRs, rv$filtFs, rv$filtRs)
    add_log(5, "Merging paired reads...")
    set_progress(5, 0, "Merging paired reads...", "running")
    shinyjs::disable("btn_merge")
    rv$bg_merge_start <- Sys.time()
    rv$merge_stage <- 1  # 1=merge+seqtab, 2=chimera removal

    # Stage 1: merge + build seqtab (fast, in-process)
    tryCatch({
      mergers <- mergePairs(rv$dadaFs, rv$filtFs, rv$dadaRs, rv$filtRs, verbose = FALSE)
      rv$mergers <- mergers
      add_log(5, "Paired reads merged.", "success")

      seqtab <- makeSequenceTable(mergers)
      rv$seqtab <- seqtab
      add_log(5, paste("Sequence table:", nrow(seqtab), "samples,", ncol(seqtab), "ASVs."), "info")

      # Stage 2: chimera removal (slow, run in background)
      set_progress(5, 40, "Removing chimeras (multithreaded)...", "running")
      rv$merge_stage <- 2
      rv$bg_merge_start <- Sys.time()
      rv$bg_merge <- callr::r_bg(
        function(seqtab, mt) {
          library(dada2)
          removeBimeraDenovo(seqtab, method = "consensus", multithread = mt, verbose = FALSE)
        },
        args = list(seqtab = seqtab, mt = can_multithread),
        supervise = TRUE
      )
    }, error = function(e) {
      add_log(5, paste("Error during merge:", e$message), "error")
      set_progress(5, 0, paste("Error:", e$message), "error")
      shinyjs::enable("btn_merge")
    })
  })

  # Poll for chimera removal completion
  observe({
    req(rv$bg_merge)
    if (rv$bg_merge$is_alive()) {
      elapsed <- as.numeric(difftime(Sys.time(), rv$bg_merge_start, units = "secs"))
      elapsed_txt <- if (elapsed > 60) paste0(round(elapsed/60, 1), " min") else paste0(round(elapsed), "s")
      set_progress(5, 40, paste0("Removing chimeras... (", elapsed_txt, " elapsed)"), "running")
      invalidateLater(1000)
    } else {
      tryCatch({
        seqtab_nochim <- rv$bg_merge$get_result()
        rv$seqtab_nochim <- seqtab_nochim
        pct <- round(sum(seqtab_nochim) / sum(rv$seqtab) * 100, 1)
        add_log(5, paste("Chimera removal complete.", ncol(seqtab_nochim), "ASVs retained.",
                         pct, "% of reads retained."), "success")

        # Build tracking table
        getN <- function(x) sum(getUniques(x))
        track <- cbind(
          rv$filter_out,
          sapply(rv$dadaFs, getN),
          sapply(rv$dadaRs, getN),
          sapply(rv$mergers, getN),
          rowSums(seqtab_nochim)
        )
        colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
        rownames(track) <- rv$sample_names
        rv$track <- track

        set_progress(5, 100, "Merge & chimera removal complete", "done")
        rv$completed_steps <- union(rv$completed_steps, 5)
      auto_save_session()
      }, error = function(e) {
        add_log(5, paste("Error:", e$message), "error")
        set_progress(5, 0, paste("Error:", e$message), "error")
      })
      shinyjs::enable("btn_merge")
      rv$bg_merge <- NULL
    }
  })

  output$seqtab_stats_ui <- renderUI({
    req(rv$seqtab, rv$seqtab_nochim)
    div(class = "grid-4", style = "margin-bottom: 16px;",
      div(class = "stat-card",
        div(class = "stat-value", nrow(rv$seqtab_nochim)),
        div(class = "stat-label", "Samples")
      ),
      div(class = "stat-card",
        div(class = "stat-value", ncol(rv$seqtab_nochim)),
        div(class = "stat-label", "ASVs (no chimeras)")
      ),
      div(class = "stat-card",
        div(class = "stat-value", ncol(rv$seqtab)),
        div(class = "stat-label", "ASVs (pre-chimera)")
      ),
      div(class = "stat-card",
        div(class = "stat-value", paste0(round(sum(rv$seqtab_nochim)/sum(rv$seqtab)*100, 1), "%")),
        div(class = "stat-label", "Reads Retained")
      )
    )
  })

  output$seqlen_plot <- renderPlot({
    req(rv$seqtab_nochim)
    seq_lengths <- nchar(getSequences(rv$seqtab_nochim))
    df <- data.frame(Length = seq_lengths)
    ggplot(df, aes(x = Length)) +
      geom_histogram(binwidth = 1, fill = "#3b82f6", color = "#1e40af", alpha = 0.8) +
      labs(title = "Distribution of ASV Sequence Lengths", x = "Length (bp)", y = "Count") +
      theme_minimal(base_size = 14) +
      theme(
        plot.background = element_rect(fill = "#1a2332", color = NA),
        panel.background = element_rect(fill = "#1a2332", color = NA),
        text = element_text(color = "#e8ecf4"),
        axis.text = element_text(color = "#8899b0"),
        panel.grid.major = element_line(color = "#2a3a52"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold")
      )
  })

  output$track_table <- renderDT({
    req(rv$track)
    df <- as.data.frame(rv$track)
    df$Sample <- rownames(df)
    df <- df[, c("Sample", "input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")]
    datatable(df, options = list(pageLength = 15, scrollX = TRUE, dom = 'frtip'),
              rownames = FALSE, class = 'compact')
  })

  output$track_plot <- renderPlot({
    req(rv$track)
    df <- as.data.frame(rv$track)
    df$Sample <- rownames(df)
    df_long <- df %>%
      pivot_longer(cols = -Sample, names_to = "Step", values_to = "Reads") %>%
      mutate(Step = factor(Step, levels = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")))

    ggplot(df_long, aes(x = Step, y = Reads, group = Sample, color = Sample)) +
      geom_line(alpha = 0.7, linewidth = 0.8) +
      geom_point(size = 2, alpha = 0.8) +
      labs(title = "Read Tracking Through Pipeline", x = "Pipeline Step", y = "Number of Reads") +
      theme_minimal(base_size = 14) +
      theme(
        plot.background = element_rect(fill = "#1a2332", color = NA),
        panel.background = element_rect(fill = "#1a2332", color = NA),
        legend.background = element_rect(fill = "#1a2332", color = NA),
        legend.key = element_rect(fill = "#1a2332", color = NA),
        text = element_text(color = "#e8ecf4"),
        axis.text = element_text(color = "#8899b0"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.major = element_line(color = "#2a3a52"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold")
      )
  })

  # ═══════════════════════════════════════════════════════════════════════
  # STEP 6: Taxonomy
  # ═══════════════════════════════════════════════════════════════════════

  observeEvent(input$btn_taxonomy, {
    req(rv$seqtab_nochim)
    add_log(6, "Starting taxonomy assignment...")
    set_progress(6, 0, "Preparing taxonomy assignment...", "running")
    shinyjs::disable("btn_taxonomy")

    tryCatch({
      db_dir <- input$silva_dir
      if (!dir.exists(db_dir)) dir.create(db_dir, recursive = TRUE)

      genus_file <- file.path(db_dir, basename(SILVA_GENUS_URL))
      species_file <- file.path(db_dir, basename(SILVA_SPECIES_URL))

      # Download if needed (synchronous, usually fast if cached)
      if (!file.exists(genus_file)) {
        set_progress(6, 5, "Downloading SILVA genus database...", "running")
        add_log(6, "Downloading SILVA genus training set...")
        download.file(SILVA_GENUS_URL, genus_file, mode = "wb", quiet = TRUE)
        add_log(6, "SILVA genus database downloaded.", "success")
      } else {
        add_log(6, "SILVA genus database found locally.", "info")
      }

      if (input$add_species && !file.exists(species_file)) {
        set_progress(6, 10, "Downloading SILVA species database...", "running")
        add_log(6, "Downloading SILVA species assignment set...")
        download.file(SILVA_SPECIES_URL, species_file, mode = "wb", quiet = TRUE)
        add_log(6, "SILVA species database downloaded.", "success")
      }

      set_progress(6, 15, "Launching taxonomy assignment (multithreaded)...", "running")
      rv$bg_tax_start <- Sys.time()

      # Determine which function to run in background
      if (input$tax_method == "bayesian") {
        rv$bg_tax <- callr::r_bg(
          function(seqtab, genus_file, species_file, add_species, mt) {
            library(dada2)
            tx <- assignTaxonomy(seqtab, genus_file, multithread = mt)
            if (add_species) tx <- addSpecies(tx, species_file)
            tx
          },
          args = list(
            seqtab = rv$seqtab_nochim, genus_file = genus_file,
            species_file = species_file, add_species = input$add_species,
            mt = can_multithread
          ),
          supervise = TRUE
        )
      } else {
        decipher_db <- file.path(db_dir, "SILVA_SSU_r138_2024.RData")
        if (file.exists(decipher_db)) {
          add_log(6, "Running IdTaxa in background subprocess...", "info")
          rv$bg_tax <- callr::r_bg(
            function(seqtab, decipher_db) {
              library(dada2); library(DECIPHER); library(Biostrings)
              dna <- DNAStringSet(getSequences(seqtab))
              load(decipher_db)
              ids <- IdTaxa(dna, trainingSet, strand = "top", processors = NULL, verbose = FALSE)
              ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
              taxid <- t(sapply(ids, function(x) {
                m <- match(ranks, x$rank)
                tx <- x$taxon[m]
                tx[startsWith(tx, "unclassified_")] <- NA
                tx
              }))
              colnames(taxid) <- ranks
              rownames(taxid) <- getSequences(seqtab)
              taxid
            },
            args = list(seqtab = rv$seqtab_nochim, decipher_db = decipher_db),
            supervise = TRUE
          )
        } else {
          add_log(6, paste("DECIPHER training set not found. Falling back to assignTaxonomy."), "warn")
          rv$bg_tax <- callr::r_bg(
            function(seqtab, genus_file, species_file, add_species, mt) {
              library(dada2)
              tx <- assignTaxonomy(seqtab, genus_file, multithread = mt)
              if (add_species) tx <- addSpecies(tx, species_file)
              tx
            },
            args = list(
              seqtab = rv$seqtab_nochim, genus_file = genus_file,
              species_file = species_file, add_species = input$add_species,
              mt = can_multithread
            ),
            supervise = TRUE
          )
        }
      }
    }, error = function(e) {
      add_log(6, paste("Error:", e$message), "error")
      set_progress(6, 0, paste("Error:", e$message), "error")
      shinyjs::enable("btn_taxonomy")
    })
  })

  # Poll for taxonomy completion
  observe({
    req(rv$bg_tax)
    if (rv$bg_tax$is_alive()) {
      elapsed <- as.numeric(difftime(Sys.time(), rv$bg_tax_start, units = "secs"))
      elapsed_txt <- if (elapsed > 60) paste0(round(elapsed/60, 1), " min") else paste0(round(elapsed), "s")
      set_progress(6, 15, paste0("Assigning taxonomy... (", elapsed_txt, " elapsed)"), "running")
      invalidateLater(1000)
    } else {
      tryCatch({
        taxa <- rv$bg_tax$get_result()
        rv$taxa <- taxa
        add_log(6, paste("Taxonomy assigned to", nrow(taxa), "ASVs."), "success")
        set_progress(6, 100, "Taxonomy assignment complete", "done")
        rv$completed_steps <- union(rv$completed_steps, 6)
      auto_save_session()
      }, error = function(e) {
        add_log(6, paste("Error:", e$message), "error")
        set_progress(6, 0, paste("Error:", e$message), "error")
      })
      shinyjs::enable("btn_taxonomy")
      rv$bg_tax <- NULL
    }
  })

  output$taxa_table <- renderDT({
    req(rv$taxa)
    df <- as.data.frame(rv$taxa)
    df$ASV <- paste0("ASV", seq_len(nrow(df)))
    df <- df[, c("ASV", setdiff(names(df), "ASV"))]
    datatable(df, options = list(pageLength = 15, scrollX = TRUE, dom = 'frtip'),
              rownames = FALSE, class = 'compact')
  })

  # ═══════════════════════════════════════════════════════════════════════
  # STEP 7: Phyloseq & Visualization
  # ═══════════════════════════════════════════════════════════════════════

  # ── Metadata variable type selection (factor vs continuous) ──
  output$metadata_vartype_ui <- renderUI({
    req(input$metadata_file)
    meta_path <- input$metadata_file$datapath
    ext <- tools::file_ext(input$metadata_file$name)
    samdf <- tryCatch({
      if (ext %in% c("csv")) read.csv(meta_path, row.names = 1, stringsAsFactors = FALSE)
      else read.delim(meta_path, row.names = 1, stringsAsFactors = FALSE)
    }, error = function(e) NULL)
    if (is.null(samdf) || ncol(samdf) == 0) return(NULL)

    rv$meta_preview <- samdf
    var_names <- colnames(samdf)

    # Auto-detect: if numeric with many unique values -> continuous, else factor
    auto_types <- sapply(var_names, function(v) {
      col <- samdf[[v]]
      if (is.numeric(col) && length(unique(col)) > 5) "continuous" else "factor"
    })

    div(
      div(class = "card-header", style = "margin-bottom: 8px;",
        div(class = "icon amber", icon("sliders-h")),
        "Variable Types"
      ),
      div(class = "card-description",
        "Specify which metadata variables should be treated as categorical (factor) vs continuous."
      ),
      div(style = "display: grid; grid-template-columns: 1fr 1fr; gap: 8px;",
        lapply(var_names, function(v) {
          div(style = "display: flex; align-items: center; gap: 8px; padding: 6px 0;",
            tags$span(style = "font-size: 13px; color: var(--text-primary); min-width: 120px;
                              font-family: 'JetBrains Mono', monospace;", v),
            radioButtons(paste0("vartype_", v), label = NULL, inline = TRUE,
              choices = c("Factor" = "factor", "Continuous" = "continuous"),
              selected = auto_types[v])
          )
        })
      )
    )
  })

  observeEvent(input$btn_phyloseq, {
    req(rv$seqtab_nochim, rv$taxa)
    add_log(7, "Building phyloseq object...")

    tryCatch({
      # Handle metadata
      samdf <- NULL
      if (!is.null(input$metadata_file)) {
        meta_path <- input$metadata_file$datapath
        ext <- tools::file_ext(input$metadata_file$name)
        if (ext %in% c("csv")) {
          samdf <- read.csv(meta_path, row.names = 1, stringsAsFactors = FALSE)
        } else {
          samdf <- read.delim(meta_path, row.names = 1, stringsAsFactors = FALSE)
        }

        # Apply factor/continuous types from user selection
        for (v in colnames(samdf)) {
          vartype_input <- input[[paste0("vartype_", v)]]
          if (!is.null(vartype_input)) {
            if (vartype_input == "factor") {
              samdf[[v]] <- as.factor(samdf[[v]])
            } else {
              samdf[[v]] <- as.numeric(as.character(samdf[[v]]))
            }
          }
        }

        add_log(7, paste("Loaded metadata with", nrow(samdf), "samples and", ncol(samdf), "variables."), "success")
      } else {
        samdf <- data.frame(SampleID = rv$sample_names, row.names = rv$sample_names)
        add_log(7, "No metadata uploaded. Using sample names as metadata.", "info")
      }

      rv$samdf <- samdf

      # Build phyloseq
      ps <- phyloseq(
        otu_table(rv$seqtab_nochim, taxa_are_rows = FALSE),
        sample_data(samdf),
        tax_table(rv$taxa)
      )

      # Store DNA sequences and rename ASVs
      dna <- Biostrings::DNAStringSet(taxa_names(ps))
      names(dna) <- taxa_names(ps)
      ps <- merge_phyloseq(ps, dna)
      taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

      rv$ps <- ps
      add_log(7, paste("Phyloseq object created:", ntaxa(ps), "taxa,", nsamples(ps), "samples."), "success")

      # Update dropdown choices - only factor variables for grouping
      meta_vars <- colnames(samdf)
      factor_vars <- meta_vars[sapply(samdf, is.factor)]
      if (length(factor_vars) == 0) factor_vars <- meta_vars  # fallback

      updateSelectInput(session, "alpha_x", choices = factor_vars, selected = factor_vars[1])
      updateSelectInput(session, "ord_color", choices = c("None", meta_vars), selected = ifelse(length(meta_vars) > 1, meta_vars[1], "None"))
      updateSelectInput(session, "pcoa_color", choices = c("None", meta_vars), selected = ifelse(length(meta_vars) > 1, meta_vars[1], "None"))
      updateSelectInput(session, "rare_color", choices = c("None", meta_vars), selected = ifelse(length(meta_vars) > 1, meta_vars[1], "None"))
      updateSelectInput(session, "bar_x", choices = factor_vars, selected = factor_vars[1])
      updateSelectInput(session, "ancom_group", choices = factor_vars, selected = factor_vars[1])

      rv$completed_steps <- union(rv$completed_steps, 7)
      auto_save_session()
    }, error = function(e) {
      add_log(7, paste("Error:", e$message), "error")
    })
  })

  # ── Transform UI (shown after phyloseq is built) ──
  output$transform_ui <- renderUI({
    req(rv$ps)
    min_depth <- min(sample_sums(rv$ps))
    tagList(
      div(class = "grid-2",
        div(
          radioButtons("transform_method", "Transformation Method",
            choices = c("Rarefaction (subsample to even depth)" = "rarefy",
                        "Relative Abundance (proportions)" = "relative",
                        "CLR (Centered Log-Ratio)" = "clr",
                        "No transformation (raw counts)" = "none"),
            selected = "rarefy")
        ),
        div(
          conditionalPanel(
            condition = "input.transform_method == 'rarefy'",
            numericInput("rarefy_depth", "Rarefaction Depth",
                         value = min_depth, min = 1, step = 100),
            div(style = "font-size: 12px; color: var(--text-muted); margin-top: 4px;",
              paste0("Minimum sample depth: ", format(min_depth, big.mark = ","),
                     " reads. Samples below this depth will be removed."))
          ),
          conditionalPanel(
            condition = "input.transform_method == 'clr'",
            div(style = "font-size: 12px; color: var(--text-muted); margin-top: 8px;",
              "Centered Log-Ratio transformation addresses compositionality. Recommended distance method for ordination: Euclidean (Aitchison distance).")
          )
        )
      )
    )
  })

  output$transform_status <- renderUI({
    req(rv$ps_transformed)
    method <- rv$transform_method
    n_samples <- nsamples(rv$ps_transformed)
    n_taxa <- ntaxa(rv$ps_transformed)
    method_label <- switch(method,
      "rarefy" = paste0("Rarefied to ", format(min(sample_sums(rv$ps_transformed)), big.mark = ","), " reads/sample"),
      "relative" = "Relative abundance (proportions)",
      "clr" = "Centered Log-Ratio (CLR) transformation",
      "none" = "Raw counts (no transformation)",
      method
    )
    div(class = "stat-card", style = "margin-top: 16px; text-align: left; padding: 14px 18px;",
      div(style = "display: flex; align-items: center; gap: 8px; margin-bottom: 4px;",
        icon("check-circle", style = "color: var(--accent-emerald);"),
        span(style = "font-size: 14px; font-weight: 500; color: var(--accent-emerald);",
             "Transformation Applied")
      ),
      div(style = "font-size: 13px; color: var(--text-secondary); font-family: 'JetBrains Mono', monospace;",
        paste0(method_label, " \u2014 ", n_samples, " samples, ", n_taxa, " ASVs"))
    )
  })

  observeEvent(input$btn_transform, {
    req(rv$ps)
    method <- input$transform_method
    add_log(7, paste("Applying transformation:", method, "..."))

    tryCatch({
      if (method == "rarefy") {
        depth <- input$rarefy_depth
        ps_t <- rarefy_even_depth(rv$ps, sample.size = depth,
                                   rngseed = 42, replace = FALSE, trimOTUs = TRUE, verbose = FALSE)
        add_log(7, paste("Rarefied to", depth, "reads/sample.", nsamples(ps_t), "samples,",
                         ntaxa(ps_t), "ASVs retained."), "success")
      } else if (method == "relative") {
        ps_t <- transform_sample_counts(rv$ps, function(x) x / sum(x))
        add_log(7, paste("Converted to relative abundance.", nsamples(ps_t), "samples."), "success")
      } else if (method == "clr") {
        ps_t <- microbiome::transform(rv$ps, "clr")
        add_log(7, paste("CLR transformation applied.", nsamples(ps_t), "samples."), "success")
      } else {
        ps_t <- rv$ps
        add_log(7, "Using raw counts (no transformation).", "info")
      }
      rv$ps_transformed <- ps_t
      rv$transform_method <- method
    }, error = function(e) {
      add_log(7, paste("Transformation error:", e$message), "error")
    })
  })

  # ── Helper: get the active phyloseq object (transformed if available) ──
  get_active_ps <- function() {
    if (!is.null(rv$ps_transformed)) rv$ps_transformed else rv$ps
  }

  # ── Rarefaction Curves ──
  make_rarefaction_plot <- function(ps, color_var, step_size) {
    # Use raw (untransformed) phyloseq for rarefaction curves
    ps_raw <- rv$ps
    otu <- as(otu_table(ps_raw), "matrix")
    if (taxa_are_rows(ps_raw)) otu <- t(otu)

    sdata <- as(sample_data(ps_raw), "data.frame")
    max_depth <- max(rowSums(otu))
    steps <- seq(1, max_depth, by = step_size)

    # Calculate rarefied richness at each depth for each sample
    rare_list <- lapply(seq_len(nrow(otu)), function(i) {
      sample_counts <- otu[i, ]
      total <- sum(sample_counts)
      richness <- sapply(steps[steps <= total], function(d) {
        # Use rarefaction formula: S_rare = S - sum(choose(N-Ni, d) / choose(N, d))
        # Approximate with vegan::rarefy if available, else simple subsampling estimate
        sum(1 - exp(lchoose(total - sample_counts[sample_counts > 0], d) -
                    lchoose(total, d)))
      })
      data.frame(
        Sample = rownames(otu)[i],
        Depth = steps[steps <= total],
        Richness = richness,
        stringsAsFactors = FALSE
      )
    })
    rare_df <- do.call(rbind, rare_list)

    # Add metadata
    rare_df <- merge(rare_df, cbind(Sample = rownames(sdata), sdata), by = "Sample")

    cv <- if (color_var == "None") NULL else color_var

    p <- ggplot(rare_df, aes(x = Depth, y = Richness, group = Sample))
    if (!is.null(cv)) {
      p <- p + geom_line(aes_string(color = cv), alpha = 0.8, linewidth = 0.7)
    } else {
      p <- p + geom_line(alpha = 0.8, linewidth = 0.7, color = "#3b82f6")
    }
    p + labs(title = "Rarefaction Curves", x = "Sequencing Depth", y = "Observed ASVs") +
      theme_minimal(base_size = 14) +
      theme(
        plot.background = element_rect(fill = "#1a2332", color = NA),
        panel.background = element_rect(fill = "#1a2332", color = NA),
        legend.background = element_rect(fill = "#1a2332", color = NA),
        legend.key = element_rect(fill = "#1a2332", color = NA),
        text = element_text(color = "#e8ecf4"),
        axis.text = element_text(color = "#8899b0"),
        panel.grid.major = element_line(color = "#2a3a52"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold")
      )
  }

  output$rarefaction_plot <- renderPlot({
    req(rv$ps, input$rare_color)
    make_rarefaction_plot(get_active_ps(), input$rare_color, input$rare_step)
  })

  # ── Alpha Diversity Plot (manual calculation, one point per sample) ──
  make_alpha_plot <- function(ps, x_var) {
    # Calculate diversity metrics manually - one value per sample
    otu <- as(otu_table(ps), "matrix")
    if (taxa_are_rows(ps)) otu <- t(otu)

    # Observed = number of non-zero ASVs per sample
    observed <- apply(otu, 1, function(x) sum(x > 0))
    # Shannon
    shannon <- apply(otu, 1, function(x) {
      x <- x[x > 0]
      p <- x / sum(x)
      -sum(p * log(p))
    })
    # Simpson
    simpson <- apply(otu, 1, function(x) {
      x <- x[x > 0]
      p <- x / sum(x)
      1 - sum(p^2)
    })

    sdata <- as(sample_data(ps), "data.frame")
    df <- data.frame(
      Sample = sample_names(ps),
      Observed = observed,
      Shannon = shannon,
      Simpson = simpson,
      stringsAsFactors = FALSE
    )
    df <- cbind(df, sdata)

    df_long <- df %>%
      pivot_longer(cols = c("Observed", "Shannon", "Simpson"),
                   names_to = "Metric", values_to = "Value") %>%
      mutate(Metric = factor(Metric, levels = c("Observed", "Shannon", "Simpson")))

    group_var <- df_long[[x_var]]

    ggplot(df_long, aes(x = .data[[x_var]], y = Value, fill = .data[[x_var]])) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, color = "#8899b0") +
      geom_jitter(aes(color = .data[[x_var]]), width = 0.2, size = 3, alpha = 0.85) +
      facet_wrap(~ Metric, scales = "free_y") +
      labs(x = x_var, y = "Value") +
      theme_minimal(base_size = 14) +
      theme(
        plot.background = element_rect(fill = "#1a2332", color = NA),
        panel.background = element_rect(fill = "#1a2332", color = NA),
        legend.position = "none",
        strip.background = element_rect(fill = "#111827", color = NA),
        strip.text = element_text(color = "#e8ecf4", size = 12, face = "bold"),
        text = element_text(color = "#e8ecf4"),
        axis.text = element_text(color = "#8899b0"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.major = element_line(color = "#2a3a52"),
        panel.grid.minor = element_blank()
      )
  }

  output$alpha_plot <- renderPlot({
    req(rv$ps, input$alpha_x)
    make_alpha_plot(get_active_ps(), input$alpha_x)
  })

  # ── Ordination Plot ──
  make_ord_plot <- function(ps, color_var, distance) {
    ps_prop <- transform_sample_counts(ps, function(otu) otu / sum(otu))
    ord <- ordinate(ps_prop, method = "NMDS", distance = distance)
    cv <- if (color_var == "None") NULL else color_var
    p <- plot_ordination(ps_prop, ord, color = cv, title = paste("NMDS -", distance))
    p <- p + geom_point(size = 4, alpha = 0.8)
    # Add dashed ellipses if color variable is set
    if (!is.null(cv)) {
      p <- p + stat_ellipse(aes_string(group = cv), type = "norm",
                            linetype = "dashed", linewidth = 0.7, level = 0.95)
    }
    p + theme_minimal(base_size = 14) +
      theme(
        plot.background = element_rect(fill = "#1a2332", color = NA),
        panel.background = element_rect(fill = "#1a2332", color = NA),
        legend.background = element_rect(fill = "#1a2332", color = NA),
        legend.key = element_rect(fill = "#1a2332", color = NA),
        text = element_text(color = "#e8ecf4"),
        axis.text = element_text(color = "#8899b0"),
        panel.grid.major = element_line(color = "#2a3a52"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold")
      )
  }

  output$ordination_plot <- renderPlot({
    req(rv$ps, input$ord_color)
    tryCatch({
      make_ord_plot(get_active_ps(), input$ord_color, input$ord_distance)
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Ordination error:", e$message), col = "#f43f5e", cex = 1.2)
    })
  })

  # ── Abundance Plot ──
  make_bar_plot <- function(ps, x_var, fill_var, top_n) {
    top_taxa <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:top_n]
    ps_top <- transform_sample_counts(ps, function(OTU) OTU / sum(OTU))
    ps_top <- prune_taxa(top_taxa, ps_top)
    p <- plot_bar(ps_top, x = x_var, fill = fill_var)
    # Remove black bar outlines by overriding geom_bar color
    p$layers[[1]]$aes_params$colour <- NA
    p + theme_minimal(base_size = 14) +
      theme(
        plot.background = element_rect(fill = "#1a2332", color = NA),
        panel.background = element_rect(fill = "#1a2332", color = NA),
        legend.background = element_rect(fill = "#1a2332", color = NA),
        legend.key = element_rect(fill = "#1a2332", color = NA),
        text = element_text(color = "#e8ecf4"),
        axis.text = element_text(color = "#8899b0"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(color = "#2a3a52"),
        panel.grid.minor = element_blank()
      )
  }

  output$bar_plot <- renderPlot({
    req(rv$ps, input$bar_x, input$bar_fill, input$bar_top_n)
    make_bar_plot(get_active_ps(), input$bar_x, input$bar_fill, input$bar_top_n)
  })

  # ── PCoA Ordination Plot ──
  make_pcoa_plot <- function(ps, color_var, distance) {
    ps_prop <- transform_sample_counts(ps, function(otu) otu / sum(otu))
    ord <- ordinate(ps_prop, method = "PCoA", distance = distance)
    cv <- if (color_var == "None") NULL else color_var
    # Get axis labels with variance explained
    evals <- ord$values$Eigenvalues
    var_explained <- round(evals / sum(evals) * 100, 1)
    ax1_lab <- paste0("PCoA1 [", var_explained[1], "%]")
    ax2_lab <- paste0("PCoA2 [", var_explained[2], "%]")

    p <- plot_ordination(ps_prop, ord, color = cv, title = paste("PCoA -", distance))
    p <- p + geom_point(size = 4, alpha = 0.8)
    # Add dashed ellipses if color variable is set
    if (!is.null(cv)) {
      p <- p + stat_ellipse(aes_string(group = cv), type = "norm",
                            linetype = "dashed", linewidth = 0.7, level = 0.95)
    }
    p + labs(x = ax1_lab, y = ax2_lab) +
      theme_minimal(base_size = 14) +
      theme(
        plot.background = element_rect(fill = "#1a2332", color = NA),
        panel.background = element_rect(fill = "#1a2332", color = NA),
        legend.background = element_rect(fill = "#1a2332", color = NA),
        legend.key = element_rect(fill = "#1a2332", color = NA),
        text = element_text(color = "#e8ecf4"),
        axis.text = element_text(color = "#8899b0"),
        panel.grid.major = element_line(color = "#2a3a52"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold")
      )
  }

  output$pcoa_plot <- renderPlot({
    req(rv$ps, input$pcoa_color)
    tryCatch({
      make_pcoa_plot(get_active_ps(), input$pcoa_color, input$pcoa_distance)
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("PCoA error:", e$message), col = "#f43f5e", cex = 1.2)
    })
  })

  # ── Downloads ──────────────────────────────────────────────────────────

  output$dl_asv <- downloadHandler(
    filename = function() paste0("ASV_table_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$seqtab_nochim)
      write.csv(rv$seqtab_nochim, file)
    }
  )

  output$dl_taxa <- downloadHandler(
    filename = function() paste0("Taxonomy_table_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$taxa)
      write.csv(rv$taxa, file)
    }
  )

  output$dl_track <- downloadHandler(
    filename = function() paste0("Read_tracking_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$track)
      write.csv(rv$track, file)
    }
  )

  output$dl_rdata <- downloadHandler(
    filename = function() paste0("dada2_workspace_", Sys.Date(), ".RData"),
    content = function(file) {
      seqtab_nochim <- rv$seqtab_nochim
      taxa <- rv$taxa
      track <- rv$track
      ps <- rv$ps
      sample_names <- rv$sample_names
      filter_out <- rv$filter_out
      errF <- rv$errF
      errR <- rv$errR
      dadaFs <- rv$dadaFs
      dadaRs <- rv$dadaRs
      mergers <- rv$mergers
      seqtab <- rv$seqtab
      samdf <- rv$samdf
      save(seqtab_nochim, taxa, track, ps, sample_names,
           filter_out, errF, errR, dadaFs, dadaRs,
           mergers, seqtab, samdf,
           file = file)
    }
  )

  # ── Step 5 downloads (tracking table + plot) ──
  output$dl_track_step5 <- downloadHandler(
    filename = function() paste0("Read_tracking_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$track)
      write.csv(rv$track, file)
    }
  )

  output$dl_track_plot <- downloadHandler(
    filename = function() paste0("Read_tracking_plot_", Sys.Date(), ".pdf"),
    content = function(file) {
      req(rv$track)
      df <- as.data.frame(rv$track)
      df$Sample <- rownames(df)
      df_long <- df %>%
        pivot_longer(cols = -Sample, names_to = "Step", values_to = "Reads") %>%
        mutate(Step = factor(Step, levels = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")))
      p <- ggplot(df_long, aes(x = Step, y = Reads, group = Sample, color = Sample)) +
        geom_line(alpha = 0.7, linewidth = 0.8) +
        geom_point(size = 2, alpha = 0.8) +
        labs(title = "Read Tracking Through Pipeline", x = "Pipeline Step", y = "Number of Reads") +
        theme_minimal(base_size = 14) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
      ggsave(file, plot = p, width = 10, height = 6)
    }
  )

  # ── Step 6 downloads (ASV + taxonomy) ──
  output$dl_asv_step6 <- downloadHandler(
    filename = function() paste0("ASV_table_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$seqtab_nochim)
      write.csv(rv$seqtab_nochim, file)
    }
  )

  output$dl_taxa_step6 <- downloadHandler(
    filename = function() paste0("Taxonomy_table_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$taxa)
      write.csv(rv$taxa, file)
    }
  )

  # ── Phyloseq plot PNG downloads ──
  output$dl_alpha_png <- downloadHandler(
    filename = function() paste0("Alpha_Diversity_", Sys.Date(), ".png"),
    content = function(file) {
      req(rv$ps, input$alpha_x)
      p <- make_alpha_plot(get_active_ps(), input$alpha_x)
      ggsave(file, plot = p, width = 12, height = 6, dpi = 300, bg = "#1a2332")
    }
  )

  output$dl_ord_png <- downloadHandler(
    filename = function() paste0("NMDS_Ordination_", Sys.Date(), ".png"),
    content = function(file) {
      req(rv$ps, input$ord_color)
      tryCatch({
        p <- make_ord_plot(get_active_ps(), input$ord_color, input$ord_distance)
        ggsave(file, plot = p, width = 10, height = 8, dpi = 300, bg = "#1a2332")
      }, error = function(e) {
        png(file, width = 800, height = 600)
        plot.new(); text(0.5, 0.5, paste("Error:", e$message)); dev.off()
      })
    }
  )

  output$dl_bar_png <- downloadHandler(
    filename = function() paste0("Abundance_Plot_", Sys.Date(), ".png"),
    content = function(file) {
      req(rv$ps, input$bar_x, input$bar_fill, input$bar_top_n)
      p <- make_bar_plot(get_active_ps(), input$bar_x, input$bar_fill, input$bar_top_n)
      ggsave(file, plot = p, width = 12, height = 8, dpi = 300, bg = "#1a2332")
    }
  )

  output$dl_pcoa_png <- downloadHandler(
    filename = function() paste0("PCoA_Ordination_", Sys.Date(), ".png"),
    content = function(file) {
      req(rv$ps, input$pcoa_color)
      tryCatch({
        p <- make_pcoa_plot(get_active_ps(), input$pcoa_color, input$pcoa_distance)
        ggsave(file, plot = p, width = 10, height = 8, dpi = 300, bg = "#1a2332")
      }, error = function(e) {
        png(file, width = 800, height = 600)
        plot.new(); text(0.5, 0.5, paste("Error:", e$message)); dev.off()
      })
    }
  )

  output$dl_rarefaction_png <- downloadHandler(
    filename = function() paste0("Rarefaction_Curves_", Sys.Date(), ".png"),
    content = function(file) {
      req(rv$ps, input$rare_color)
      p <- make_rarefaction_plot(get_active_ps(), input$rare_color, input$rare_step)
      ggsave(file, plot = p, width = 12, height = 7, dpi = 300, bg = "#1a2332")
    }
  )

  # ── Phyloseq object download (.rds) - saves transformed version ──
  output$dl_phyloseq_rds <- downloadHandler(
    filename = function() paste0("phyloseq_object_", Sys.Date(), ".rds"),
    content = function(file) {
      req(rv$ps)
      saveRDS(get_active_ps(), file)
    }
  )

  # ═══════════════════════════════════════════════════════════════════════
  # STEP 9: PERMANOVA
  # ═══════════════════════════════════════════════════════════════════════

  # Command preview
  output$permanova_cmd_preview <- renderUI({
    formula_text <- trimws(input$permanova_formula)
    dist_method <- input$permanova_distance
    n_perm <- input$permanova_perm
    if (!nzchar(formula_text)) return(NULL)

    cmd1 <- paste0("dist_mat <- vegdist(otu_table, method = \"", dist_method, "\")")
    cmd2 <- paste0("adonis2(dist_mat ~ ", formula_text, ", data = sample_data, permutations = ", n_perm, ")")
    # betadisper uses the first variable in the formula
    first_var <- trimws(strsplit(formula_text, "\\+|\\*|:")[[1]][1])
    cmd3 <- paste0("betadisper(dist_mat, groups = sample_data$", first_var, ")")

    div(class = "log-panel", style = "margin-top: 12px; max-height: none; padding: 14px;",
      div(style = "font-size: 11px; color: var(--text-muted); text-transform: uppercase; letter-spacing: 0.05em; margin-bottom: 8px;",
        "Command Preview"),
      div(class = "log-info", style = "margin-bottom: 4px;", cmd1),
      div(class = "log-info", style = "margin-bottom: 4px;", cmd2),
      div(class = "log-info", cmd3)
    )
  })

  observeEvent(input$btn_permanova, {
    req(rv$ps)
    formula_text <- trimws(input$permanova_formula)
    if (!nzchar(formula_text)) {
      add_log(9, "Please enter a formula.", "error")
      return()
    }
    add_log(9, paste("Running PERMANOVA with formula:", formula_text, "..."))

    tryCatch({
      ps_active <- get_active_ps()
      sdata <- as(sample_data(ps_active), "data.frame")
      dist_method <- input$permanova_distance
      n_perm <- input$permanova_perm

      # Get OTU matrix
      otu <- as(otu_table(ps_active), "matrix")
      if (taxa_are_rows(ps_active)) otu <- t(otu)

      # Calculate distance matrix
      dist_mat <- vegdist(otu, method = dist_method)

      # Run PERMANOVA (adonis2) with user formula
      formula_str <- as.formula(paste("dist_mat ~", formula_text))
      perm_result <- adonis2(formula_str, data = sdata, permutations = n_perm)
      rv$permanova_result <- perm_result
      add_log(9, paste("PERMANOVA complete. p-value:", round(perm_result$`Pr(>F)`[1], 3)), "success")

      # Run betadisper using the first variable in the formula
      first_var <- trimws(strsplit(formula_text, "\\+|\\*|:")[[1]][1])
      if (first_var %in% colnames(sdata)) {
        groups <- sdata[[first_var]]
        bd <- betadisper(dist_mat, groups)
        bd_test <- permutest(bd, pairwise = TRUE, permutations = n_perm)
        rv$betadisper_result <- list(betadisper = bd, permutest = bd_test)
        add_log(9, paste("Betadisper complete (on", first_var, "). p-value:", round(bd_test$tab$`Pr(>F)`[1], 3)), "success")

        # Pairwise PERMANOVA if >2 groups (using first variable only)
        group_levels <- levels(as.factor(groups))
        if (length(group_levels) > 2) {
          add_log(9, "Running pairwise PERMANOVA...")
          pairs <- combn(group_levels, 2, simplify = FALSE)
          pair_results <- lapply(pairs, function(pair) {
            idx <- groups %in% pair
            sub_otu <- otu[idx, , drop = FALSE]
            sub_sdata <- sdata[idx, , drop = FALSE]
            sub_dist <- vegdist(sub_otu, method = dist_method)
            sub_formula <- as.formula(paste("sub_dist ~", first_var))
            res <- adonis2(sub_formula, data = sub_sdata, permutations = n_perm)
            data.frame(
              Group1 = pair[1], Group2 = pair[2],
              F_value = round(res$F[1], 3),
              R2 = round(res$R2[1], 3),
              p_value = round(res$`Pr(>F)`[1], 3),
              stringsAsFactors = FALSE
            )
          })
          pairwise_df <- do.call(rbind, pair_results)
          pairwise_df$p_adjusted <- round(p.adjust(pairwise_df$p_value, method = "bonferroni"), 3)
          rv$pairwise_result <- pairwise_df
          add_log(9, paste("Pairwise PERMANOVA complete.", nrow(pairwise_df), "comparisons."), "success")
        } else {
          rv$pairwise_result <- NULL
          add_log(9, "Only 2 groups - pairwise PERMANOVA not needed.", "info")
        }
      } else {
        add_log(9, paste("Variable", first_var, "not found in metadata. Betadisper skipped."), "warn")
        rv$betadisper_result <- NULL
        rv$pairwise_result <- NULL
      }

      rv$completed_steps <- union(rv$completed_steps, 9)
      auto_save_session()
    }, error = function(e) {
      add_log(9, paste("PERMANOVA error:", e$message), "error")
    })
  })

  # ── PERMANOVA results renderers ──
  output$permanova_results_ui <- renderUI({
    req(rv$permanova_result)
    res <- rv$permanova_result
    df <- as.data.frame(res)
    df$Term <- rownames(df)
    df <- df[, c("Term", "Df", "SumOfSqs", "R2", "F", "Pr(>F)")]
    colnames(df) <- c("Term", "Df", "Sum of Squares", "R2", "F statistic", "p-value")

    # Round all numeric columns to 3 decimal places
    num_cols <- c("Sum of Squares", "R2", "F statistic", "p-value")
    for (col in num_cols) df[[col]] <- round(df[[col]], 3)

    sig <- df$`p-value`[1]
    sig_text <- if (!is.na(sig) && sig < 0.001) "***" else if (!is.na(sig) && sig < 0.01) "**" else if (!is.na(sig) && sig < 0.05) "*" else "ns"
    sig_color <- if (!is.na(sig) && sig < 0.05) "var(--accent-emerald)" else "var(--accent-amber)"

    tagList(
      div(class = "stat-card", style = "text-align: left; padding: 14px 18px; margin-bottom: 16px;",
        div(style = paste0("font-size: 14px; font-weight: 600; color: ", sig_color, ";"),
          paste0("p-value = ", format(sig, nsmall = 3), "  ", sig_text)
        ),
        div(style = "font-size: 13px; color: var(--text-secondary); margin-top: 4px;",
          paste0("R2 = ", format(df$R2[1], nsmall = 3),
                 " | F = ", format(df$`F statistic`[1], nsmall = 3),
                 " | Permutations: ", input$permanova_perm)
        )
      ),
      DTOutput("permanova_dt")
    )
  })

  output$permanova_dt <- renderDT({
    req(rv$permanova_result)
    res <- rv$permanova_result
    df <- as.data.frame(res)
    df$Term <- rownames(df)
    df <- df[, c("Term", "Df", "SumOfSqs", "R2", "F", "Pr(>F)")]
    colnames(df) <- c("Term", "Df", "Sum of Squares", "R2", "F statistic", "p-value")
    num_cols <- c("Sum of Squares", "R2", "F statistic", "p-value")
    for (col in num_cols) df[[col]] <- round(df[[col]], 3)
    datatable(df, options = list(dom = "t", pageLength = 10), rownames = FALSE, class = "compact")
  })

  output$betadisper_results_ui <- renderUI({
    req(rv$betadisper_result)
    bd_test <- rv$betadisper_result$permutest
    tab <- as.data.frame(bd_test$tab)
    tab$Term <- rownames(tab)
    tab <- tab[, c("Term", "Df", "Sum Sq", "Mean Sq", "F", "N.Perm", "Pr(>F)")]

    sig <- tab$`Pr(>F)`[1]
    sig_color <- if (!is.na(sig) && sig < 0.05) "var(--accent-rose)" else "var(--accent-emerald)"
    interp <- if (!is.na(sig) && sig < 0.05) {
      "Significant - group dispersions are NOT homogeneous. PERMANOVA results should be interpreted with caution."
    } else {
      "Not significant - group dispersions are homogeneous. PERMANOVA results are reliable."
    }

    tagList(
      div(class = "stat-card", style = "text-align: left; padding: 14px 18px; margin-bottom: 16px;",
        div(style = paste0("font-size: 14px; font-weight: 600; color: ", sig_color, ";"),
          paste0("Betadisper p-value = ", format(round(sig, 3), nsmall = 3))
        ),
        div(style = "font-size: 13px; color: var(--text-secondary); margin-top: 4px;", interp)
      ),
      DTOutput("betadisper_dt")
    )
  })

  output$betadisper_dt <- renderDT({
    req(rv$betadisper_result)
    bd_test <- rv$betadisper_result$permutest
    tab <- as.data.frame(bd_test$tab)
    tab$Term <- rownames(tab)
    tab <- tab[, c("Term", "Df", "Sum Sq", "Mean Sq", "F", "N.Perm", "Pr(>F)")]
    num_cols <- c("Sum Sq", "Mean Sq", "F", "Pr(>F)")
    for (col in num_cols) tab[[col]] <- round(tab[[col]], 3)
    datatable(tab, options = list(dom = "t", pageLength = 10), rownames = FALSE, class = "compact")
  })

  output$pairwise_permanova_card <- renderUI({
    req(rv$pairwise_result)
    div(class = "card",
      div(class = "card-header",
        div(class = "icon cyan", icon("exchange-alt")),
        "Pairwise PERMANOVA"
      ),
      div(class = "card-description",
        "Pairwise comparisons between all group pairs with Bonferroni-adjusted p-values."
      ),
      DTOutput("pairwise_dt")
    )
  })

  output$pairwise_dt <- renderDT({
    req(rv$pairwise_result)
    datatable(rv$pairwise_result, options = list(dom = "t", pageLength = 20, scrollX = TRUE),
              rownames = FALSE, class = "compact")
  })

  # ── PERMANOVA download handlers ──
  output$dl_permanova_csv <- downloadHandler(
    filename = function() paste0("PERMANOVA_results_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$permanova_result)
      df <- as.data.frame(rv$permanova_result)
      df$Term <- rownames(df)
      write.csv(df, file, row.names = FALSE)
    }
  )

  output$dl_betadisper_csv <- downloadHandler(
    filename = function() paste0("Betadisper_results_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$betadisper_result)
      df <- as.data.frame(rv$betadisper_result$permutest$tab)
      df$Term <- rownames(df)
      write.csv(df, file, row.names = FALSE)
    }
  )

  output$dl_pairwise_csv <- downloadHandler(
    filename = function() paste0("Pairwise_PERMANOVA_", Sys.Date(), ".csv"),
    content = function(file) {
      if (!is.null(rv$pairwise_result)) {
        write.csv(rv$pairwise_result, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Note = "Pairwise PERMANOVA not applicable (2 or fewer groups)"), file, row.names = FALSE)
      }
    }
  )

  # =====================================================================
  # STEP 10: ANCOM-BC2
  # =====================================================================

  # Command preview
  output$ancom_cmd_preview <- renderUI({
    fix_f <- trimws(input$ancom_fix_formula)
    rand_f <- trimws(input$ancom_rand_formula)
    grp <- input$ancom_group
    tax <- input$ancom_tax_level
    if (!nzchar(fix_f)) return(NULL)

    rand_line <- if (nzchar(rand_f)) paste0('         rand_formula = "', rand_f, '",\n') else ""

    cmd <- paste0('ancombc2(data = ps, tax_level = "', tax, '",\n',
                  '         fix_formula = "', fix_f, '",\n',
                  rand_line,
                  '         group = "', grp, '",\n',
                  '         p_adj_method = "', input$ancom_p_adj, '",\n',
                  '         pseudo = ', input$ancom_pseudo, ', pseudo_sens = ', input$ancom_pseudo_sens, ',\n',
                  '         prv_cut = ', input$ancom_prv_cut, ', lib_cut = ', input$ancom_lib_cut, ',\n',
                  '         s0_perc = ', input$ancom_s0_perc, ', alpha = ', input$ancom_alpha, ',\n',
                  '         struc_zero = ', input$ancom_struc_zero, ', neg_lb = ', input$ancom_neg_lb, ',\n',
                  '         n_cl = ', input$ancom_n_cl, ',\n',
                  '         global = ', input$ancom_global, ', pairwise = ', input$ancom_pairwise, ',\n',
                  '         dunnet = ', input$ancom_dunnet, ', trend = ', input$ancom_trend, ')')

    div(class = "log-panel", style = "margin-top: 12px; max-height: none; padding: 14px;",
      div(style = "font-size: 11px; color: var(--text-muted); text-transform: uppercase; letter-spacing: 0.05em; margin-bottom: 8px;",
        "Command Preview"),
      div(class = "log-info", tags$pre(style = "margin: 0; white-space: pre-wrap; color: var(--accent-cyan); font-size: 12px;", cmd))
    )
  })

  # Run ANCOM-BC2 asynchronously
  observeEvent(input$btn_ancombc2, {
    req(rv$ps, input$ancom_fix_formula)
    fix_f <- trimws(input$ancom_fix_formula)
    rand_f <- trimws(input$ancom_rand_formula)
    grp <- input$ancom_group
    tax <- input$ancom_tax_level

    if (!nzchar(fix_f)) {
      add_log(10, "Please enter a fixed effects formula.", "error")
      return()
    }

    rand_formula <- if (nzchar(rand_f)) rand_f else NULL

    add_log(10, paste("Running ANCOM-BC2 at", tax, "level with fix_formula:", fix_f,
                      if (!is.null(rand_formula)) paste("| rand_formula:", rand_formula) else ""))
    set_progress(10, 0, "Launching ANCOM-BC2 (this may take several minutes)...", "running")
    shinyjs::disable("btn_ancombc2")

    rv$bg_ancom_start <- Sys.time()
    rv$bg_ancom <- callr::r_bg(
      function(ps, tax_level, fix_formula, rand_formula, group,
               p_adj_method, pseudo, pseudo_sens, prv_cut, lib_cut, s0_perc,
               struc_zero, neg_lb, alpha, n_cl,
               global, pairwise, dunnet, trend,
               iter_tol, iter_max, em_tol, em_max,
               mdfdr_method, mdfdr_B) {
        library(ANCOMBC)
        library(phyloseq)
        ancombc2(data = ps, tax_level = tax_level,
                 fix_formula = fix_formula,
                 rand_formula = rand_formula,
                 group = group,
                 p_adj_method = p_adj_method,
                 pseudo = pseudo, pseudo_sens = pseudo_sens,
                 prv_cut = prv_cut, lib_cut = lib_cut,
                 s0_perc = s0_perc,
                 struc_zero = struc_zero, neg_lb = neg_lb,
                 alpha = alpha, n_cl = n_cl, verbose = FALSE,
                 global = global, pairwise = pairwise,
                 dunnet = dunnet, trend = trend,
                 iter_control = list(tol = iter_tol, max_iter = iter_max, verbose = FALSE),
                 em_control = list(tol = em_tol, max_iter = em_max),
                 mdfdr_control = list(fwer_ctrl_method = mdfdr_method, B = mdfdr_B))
      },
      args = list(
        ps = rv$ps,
        tax_level = tax,
        fix_formula = fix_f,
        rand_formula = rand_formula,
        group = grp,
        p_adj_method = input$ancom_p_adj,
        pseudo = input$ancom_pseudo,
        pseudo_sens = input$ancom_pseudo_sens,
        prv_cut = input$ancom_prv_cut,
        lib_cut = input$ancom_lib_cut,
        s0_perc = input$ancom_s0_perc,
        struc_zero = input$ancom_struc_zero,
        neg_lb = input$ancom_neg_lb,
        alpha = input$ancom_alpha,
        n_cl = input$ancom_n_cl,
        global = input$ancom_global,
        pairwise = input$ancom_pairwise,
        dunnet = input$ancom_dunnet,
        trend = input$ancom_trend,
        iter_tol = input$ancom_iter_tol,
        iter_max = input$ancom_iter_max,
        em_tol = input$ancom_em_tol,
        em_max = input$ancom_em_max,
        mdfdr_method = input$ancom_mdfdr_method,
        mdfdr_B = input$ancom_mdfdr_B
      ),
      supervise = TRUE
    )
  })

  # Poll for ANCOM-BC2 completion
  observe({
    req(rv$bg_ancom)
    if (rv$bg_ancom$is_alive()) {
      elapsed <- as.numeric(difftime(Sys.time(), rv$bg_ancom_start, units = "secs"))
      elapsed_txt <- if (elapsed > 60) paste0(round(elapsed/60, 1), " min") else paste0(round(elapsed), "s")
      set_progress(10, 0, paste0("ANCOM-BC2 running... (", elapsed_txt, " elapsed)"), "running")
      invalidateLater(2000)
    } else {
      tryCatch({
        result <- rv$bg_ancom$get_result()
        rv$ancom_result <- result
        n_global <- if (!is.null(result$res_global)) sum(result$res_global$diff_abn, na.rm = TRUE) else 0
        add_log(10, paste("ANCOM-BC2 complete.", n_global, "globally differentially abundant taxa detected."), "success")
        set_progress(10, 100, "ANCOM-BC2 complete", "done")
        rv$completed_steps <- union(rv$completed_steps, 10)
        auto_save_session()
      }, error = function(e) {
        add_log(10, paste("ANCOM-BC2 error:", e$message), "error")
        set_progress(10, 0, paste("Error:", e$message), "error")
      })
      shinyjs::enable("btn_ancombc2")
      rv$bg_ancom <- NULL
    }
  })

  # ── Global test table ──
  output$ancom_global_dt <- renderDT({
    req(rv$ancom_result)
    res <- rv$ancom_result$res_global
    if (is.null(res) || nrow(res) == 0) {
      # Fallback: show primary result summary
      res <- rv$ancom_result$res
      if (is.null(res) || nrow(res) == 0) return(NULL)
    }
    num_cols <- sapply(res, is.numeric)
    res[num_cols] <- lapply(res[num_cols], round, 3)
    datatable(res, options = list(pageLength = 15, scrollX = TRUE, dom = "frtip"),
              rownames = FALSE, class = "compact")
  })

  # ── Pairwise directional test table ──
  output$ancom_pairwise_dt <- renderDT({
    req(rv$ancom_result)
    res <- rv$ancom_result$res_pair
    if (is.null(res) || nrow(res) == 0) {
      # Fallback: show primary result
      res <- rv$ancom_result$res
      if (is.null(res) || nrow(res) == 0) return(NULL)
    }
    num_cols <- sapply(res, is.numeric)
    res[num_cols] <- lapply(res[num_cols], round, 3)
    datatable(res, options = list(pageLength = 15, scrollX = TRUE, dom = "frtip"),
              rownames = FALSE, class = "compact")
  })

  # ── Volcano plot ──
  output$ancom_volcano_selector <- renderUI({
    req(rv$ancom_result)
    # Try res_pair first, then fall back to res (primary)
    res <- rv$ancom_result$res_pair
    if (is.null(res) || nrow(res) == 0) res <- rv$ancom_result$res
    if (is.null(res)) return(NULL)

    # Find lfc columns
    lfc_cols <- grep("^lfc_|^lfc\\.", colnames(res), value = TRUE)
    if (length(lfc_cols) == 0) return(NULL)
    comparisons <- gsub("^lfc_|^lfc\\.", "", lfc_cols)
    selectInput("ancom_volcano_comp", "Select Comparison",
                choices = comparisons, selected = comparisons[1])
  })

  make_volcano_plot <- function(result, comparison, alpha) {
    # Try res_pair first, then res
    res <- result$res_pair
    if (is.null(res) || nrow(res) == 0) res <- result$res
    if (is.null(res)) return(NULL)

    # Try different column naming patterns
    lfc_col <- NULL
    p_col <- NULL
    diff_col <- NULL
    for (prefix in c("lfc_", "lfc.")) {
      candidate <- paste0(prefix, comparison)
      if (candidate %in% colnames(res)) { lfc_col <- candidate; break }
    }
    for (prefix in c("q_", "q.")) {
      candidate <- paste0(prefix, comparison)
      if (candidate %in% colnames(res)) { p_col <- candidate; break }
    }
    for (prefix in c("diff_", "diff.")) {
      candidate <- paste0(prefix, comparison)
      if (candidate %in% colnames(res)) { diff_col <- candidate; break }
    }

    if (is.null(lfc_col) || is.null(p_col)) return(NULL)

    df <- data.frame(
      taxon = res$taxon,
      lfc = res[[lfc_col]],
      qval = res[[p_col]],
      stringsAsFactors = FALSE
    )
    if (!is.null(diff_col)) {
      df$diff <- res[[diff_col]]
    } else {
      df$diff <- df$qval < alpha
    }
    df <- df[complete.cases(df), ]
    df$neg_log10_q <- -log10(df$qval + 1e-300)
    df$Significant <- ifelse(df$diff, "Yes", "No")

    ggplot(df, aes(x = lfc, y = neg_log10_q, color = Significant)) +
      geom_point(size = 2.5, alpha = 0.7) +
      scale_color_manual(values = c("No" = "#5a6a80", "Yes" = "#f43f5e")) +
      geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "#f59e0b", linewidth = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#8899b0", linewidth = 0.3) +
      labs(title = paste("Volcano Plot -", comparison),
           x = "Log Fold Change", y = "-log10(q-value)") +
      theme_minimal(base_size = 14) +
      theme(
        plot.background = element_rect(fill = "#1a2332", color = NA),
        panel.background = element_rect(fill = "#1a2332", color = NA),
        legend.background = element_rect(fill = "#1a2332", color = NA),
        legend.key = element_rect(fill = "#1a2332", color = NA),
        text = element_text(color = "#e8ecf4"),
        axis.text = element_text(color = "#8899b0"),
        panel.grid.major = element_line(color = "#2a3a52"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold")
      )
  }

  output$ancom_volcano_plot <- renderPlot({
    req(rv$ancom_result, input$ancom_volcano_comp)
    make_volcano_plot(rv$ancom_result, input$ancom_volcano_comp, input$ancom_alpha)
  })

  # ── Heatmap of DA taxa ──
  make_ancom_heatmap <- function(result) {
    # Try res_pair first, then res
    res <- result$res_pair
    if (is.null(res) || nrow(res) == 0) res <- result$res
    if (is.null(res)) return(NULL)

    # Find all lfc and diff columns
    lfc_cols <- grep("^lfc_|^lfc\\.", colnames(res), value = TRUE)
    diff_cols <- grep("^diff_|^diff\\.", colnames(res), value = TRUE)

    if (length(lfc_cols) == 0) return(NULL)

    # Filter to only significant taxa (significant in at least one comparison)
    sig_mask <- apply(res[diff_cols], 1, function(x) any(x, na.rm = TRUE))
    if (sum(sig_mask) == 0) {
      plot.new()
      text(0.5, 0.5, "No differentially abundant taxa detected.", col = "#8899b0", cex = 1.4)
      return()
    }

    res_sig <- res[sig_mask, ]
    lfc_mat <- as.data.frame(res_sig[lfc_cols])
    rownames(lfc_mat) <- res_sig$taxon
    colnames(lfc_mat) <- gsub("^lfc_|^lfc\\.", "", colnames(lfc_mat))

    # Reshape for ggplot
    lfc_mat$taxon <- rownames(lfc_mat)
    df_long <- tidyr::pivot_longer(lfc_mat, cols = -taxon, names_to = "Comparison", values_to = "LFC")
    df_long <- df_long[complete.cases(df_long), ]

    ggplot(df_long, aes(x = Comparison, y = taxon, fill = LFC)) +
      geom_tile(color = "#2a3a52", linewidth = 0.3) +
      scale_fill_gradient2(low = "#3b82f6", mid = "#1a2332", high = "#f43f5e", midpoint = 0,
                           name = "Log FC") +
      labs(title = "Differentially Abundant Taxa - Log Fold Changes", x = "", y = "") +
      theme_minimal(base_size = 12) +
      theme(
        plot.background = element_rect(fill = "#1a2332", color = NA),
        panel.background = element_rect(fill = "#1a2332", color = NA),
        legend.background = element_rect(fill = "#1a2332", color = NA),
        legend.key = element_rect(fill = "#1a2332", color = NA),
        text = element_text(color = "#e8ecf4"),
        axis.text = element_text(color = "#8899b0"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank(),
        plot.title = element_text(size = 14, face = "bold")
      )
  }

  output$ancom_heatmap <- renderPlot({
    req(rv$ancom_result)
    make_ancom_heatmap(rv$ancom_result)
  })

  # ── ANCOM-BC2 download handlers ──
  output$dl_ancom_global <- downloadHandler(
    filename = function() paste0("ANCOMBC2_global_test_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$ancom_result)
      res <- rv$ancom_result$res_global
      if (is.null(res) || nrow(res) == 0) res <- rv$ancom_result$res
      if (!is.null(res)) write.csv(res, file, row.names = FALSE)
      else write.csv(data.frame(Note = "No results available"), file, row.names = FALSE)
    }
  )

  output$dl_ancom_pairwise <- downloadHandler(
    filename = function() paste0("ANCOMBC2_pairwise_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$ancom_result)
      res <- rv$ancom_result$res_pair
      if (is.null(res) || nrow(res) == 0) res <- rv$ancom_result$res
      if (!is.null(res)) write.csv(res, file, row.names = FALSE)
      else write.csv(data.frame(Note = "No results available"), file, row.names = FALSE)
    }
  )

  output$dl_ancom_volcano <- downloadHandler(
    filename = function() paste0("ANCOMBC2_volcano_", Sys.Date(), ".png"),
    content = function(file) {
      req(rv$ancom_result, input$ancom_volcano_comp)
      p <- make_volcano_plot(rv$ancom_result, input$ancom_volcano_comp, input$ancom_alpha)
      ggsave(file, plot = p, width = 10, height = 7, dpi = 300, bg = "#1a2332")
    }
  )

  output$dl_ancom_heatmap <- downloadHandler(
    filename = function() paste0("ANCOMBC2_heatmap_", Sys.Date(), ".png"),
    content = function(file) {
      req(rv$ancom_result)
      p <- make_ancom_heatmap(rv$ancom_result)
      ggsave(file, plot = p, width = 12, height = max(8, sum(rv$ancom_result$res_pair$diff_abn, na.rm = TRUE) * 0.3 + 4),
             dpi = 300, bg = "#1a2332")
    }
  )

}

# ── Run app ──────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
