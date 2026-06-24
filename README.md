# RNA-seq & Bioconductor Notes 2026

## Overview

This repository contains notes, practical exercises, and code developed during an RNA-seq and Bioconductor training module led by Dr. Leonardo Collado.

The material combines lecture notes, workshop exercises, and exploratory analyses covering modern transcriptomics workflows and the Bioconductor ecosystem.

---

## Topics Covered

### RNA-seq Analysis

* Experimental design
* Technical and biological variation
* Batch effects
* RNA quality assessment
* Read alignment and quantification
* Differential expression analysis
* Statistical modeling
* Interpretation of transcriptomic results

### Bioconductor Ecosystem

* Package development philosophy
* Reproducible workflows
* Vignettes and documentation
* Release management
* Dependency management
* Community standards

### Core Bioconductor Objects

* SummarizedExperiment
* GenomicRanges
* SingleCellExperiment

Topics include:

* Assays
* Row metadata
* Genomic coordinates
* Experimental metadata
* Data integration

### Interactive Data Exploration

Using the iSEE package:

* PCA visualization
* UMAP and t-SNE exploration
* Heatmaps
* Cluster inspection
* Interactive gene expression analysis

### Statistical Concepts

* Experimental design matrices
* Model specification
* Confounding variables
* Replication
* Multiple testing correction
* False discovery rate (FDR)

### Real-World RNA-seq Challenges

Discussion of practical issues including:

* Technical variability
* Sequencing artifacts
* RNA degradation
* Mapping ambiguities
* Biological noise
* Statistical uncertainty

---

## Guest Lectures

The notes include summaries and key concepts from talks covering:

### Regulatory Genomics and Immunology

* Lupus genetics
* B-cell biology
* Regulatory variants
* eQTL studies
* ATAC-seq

### Single-cell Transcriptomics

* Single-cell RNA-seq
* Single-nucleus RNA-seq
* Alzheimer's disease research
* Cell-type specific analyses

### Bioinformatics Career and Community Topics

* Bioconductor governance
* Open-source scientific software
* Reproducibility in computational biology

---

## Repository Structure

```text
.
├── notas/
├── figuras/
├── processed-data/
├── R/
└── rnaseq_2026_notas.Rproj
```

### R Scripts

The repository includes practical examples covering:

```text
01-notas.R
04_SummarizeData.R
05_ejercicio_iSEE.R
06_recount.R
07_statistic_models.R
08_ExploreModelMatrix.R
```

These scripts demonstrate common RNA-seq and Bioconductor workflows using real datasets and standard analysis pipelines.

---

## Technologies

* R
* Bioconductor
* SummarizedExperiment
* GenomicRanges
* iSEE
* recount3
* Quarto
* GitHub Pages

---

## Skills Demonstrated

* RNA-seq analysis
* Differential expression analysis
* Statistical modeling
* Bioconductor workflows
* Data visualization
* Reproducible research
* Scientific computing in R
* Transcriptomics

---

## Purpose

This repository serves as a personal knowledge base and practical reference for transcriptomics, Bioconductor workflows, and RNA-seq data analysis.
