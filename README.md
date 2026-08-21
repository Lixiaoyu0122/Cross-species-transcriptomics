# Cross-species-transcriptomics
Code for Tanshinones attenuate early post-infarction cardiac fibrosis by suppressing FSTL1-associated ECM remodeling
Mouse bulk RNA-seq was analyzed to identify MI-induced and M-Tans-reversed transcriptional signatures. One-to-one mouse-human ortholog mapping was performed using Ensembl BioMart. The conserved signatures were projected onto human cardiac scRNA-seq and snRNA-seq datasets using Seurat module scoring. Spatial transcriptomic datasets were subsequently analyzed to determine the anatomical distribution of conserved fibroblast-associated remodeling signatures.
# Cross-species transcriptomic and spatial mapping of cardiac remodeling programs

This repository provides the computational workflow for mapping treatment-responsive fibroblast-associated remodeling signatures across human cardiac single-cell and spatial transcriptomic datasets.

## Overview

The workflow integrates:

1. Cross-omics identification of conserved fibroblast-associated remodeling signatures
2. Spatial transcriptomic localization of FSTL1-associated remodeling niches
3. Quantification and visualization of pathological cardiac regions

## Workflow

Mouse myocardial transcriptomic signatures were projected onto human cardiac datasets through gene-level correspondence and disease-associated remodeling signatures.

Human single-cell transcriptomic datasets were used to define fibroblast-associated states.

Spatial transcriptomic datasets were analyzed to determine anatomical localization of FSTL1-associated ECM remodeling signatures.

## Main outputs

- Cross-species remodeling signature comparison
- FSTL1-associated ECM enrichment analysis
- Spatial localization of pathological remodeling regions

## Data availability

Public human cardiac transcriptomic datasets were obtained from previously published studies. Raw sequencing data are not redistributed in this repository.
