############################################################
# Project configuration file
#
# Purpose:
#   Centralize all paths for reproducible analysis.
#
# Usage:
#   source("scripts/utils/project_config.R")
#
# Author:
#   Cross-species transcriptomics pipeline
#
############################################################


############################
# Project root
############################

PROJECT_ROOT <- normalizePath(
  ".",
  mustWork = FALSE
)


############################
# Data directories
############################

DATA_DIR <- file.path(
  PROJECT_ROOT,
  "data"
)


MOUSE_BULK_DIR <- file.path(
  DATA_DIR,
  "mouse_bulk"
)


MOUSE_COUNT_DIR <- file.path(
  MOUSE_BULK_DIR,
  "counts"
)


############################
# Result directories
############################

RESULT_DIR <- file.path(
  PROJECT_ROOT,
  "results"
)


MOUSE_BULK_RESULT_DIR <- file.path(
  RESULT_DIR,
  "mouse_bulk"
)


DESEQ2_QC_DIR <- file.path(
  MOUSE_BULK_RESULT_DIR,
  "deseq2_qc"
)


DEG_RESULT_DIR <- file.path(
  MOUSE_BULK_RESULT_DIR,
  "deseq2_deg"
)


GSEA_RESULT_DIR <- file.path(
  MOUSE_BULK_RESULT_DIR,
  "GSEA_GO_KEGG"
)



############################
# Create output directories
############################

dir.create(
  RESULT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  MOUSE_BULK_RESULT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  DESEQ2_QC_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  DEG_RESULT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  GSEA_RESULT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)



############################
# Analysis parameters
############################

DEG_LOG2FC_CUTOFF <- 1

DEG_FDR_CUTOFF <- 0.05


RANDOM_SEED <- 20260720



############################
# Mouse ECM/FSTL1 signature
############################

FSTL1_ECM_GENES <- c(

  "Fstl1",

  "Col1a1",
  "Col1a2",
  "Col3a1",

  "Postn",

  "Fn1",

  "Sparc"

)



############################
# Session information
############################

message(
  "Project root: ",
  PROJECT_ROOT
)

message(
  "Configuration loaded successfully."
)
