# Single cell RNA seq analysis pipeline in python
Table of Contents
-[Single‑Cell RNA‑seq Processing Pipeline](#single_cell_rna_seq_processing_pipeline) \
-[Why this pipeline exists](#why_this_pipeline_exists) \
-[What the pipeline does](#what_the_pipeline_does) \
-[Citations](#citations) \

## Single‑Cell RNA‑seq Processing Pipeline
This repository contains a modular, end‑to‑end [single‑cell RNA‑seq analysis pipeline](scripts/single_cell_analysis.py) built on Scanpy. It’s designed for real lab workflows: multiple samples, batch effects, messy QC, and the need to keep everything reproducible and well‑organized. Very easy to use just update the [single_cell_config](scripts/single_cell_config.yaml).
Every major step — QC, normalization, HVG selection, batch correction, embedding, clustering, and marker detection — is wrapped in a clean function so you can reuse or swap modules without touching the rest of the pipeline.

## Why this pipeline exists
In most labs, single‑cell analysis ends up scattered across notebooks, half‑finished scripts, and ad‑hoc plots. This pipeline aims to solve that by providing:
A repeatable workflow you can run on any new dataset
Automatic figure saving so nothing gets lost
Logging for traceability (critical when working on HPC or collaborating)
Modular functions that can be reused in other projects
Interchangeable batch correction backends depending on your experiment
It’s meant to feel like a lightweight framework rather than a one‑off script.

## What the pipeline does
1. Load and merge multiple 10x samples
The pipeline accepts any number of 10x .h5 files, assigns each a sample name, ensures unique gene/cell identifiers, and concatenates them into a single AnnData object with a batch column.
This mirrors how most wet‑lab experiments are structured: multiple replicates, conditions, or donors.

2. Quality control with visual diagnostics
QC is intentionally opinionated but adjustable.
The pipeline computes:
Mitochondrial percentage
Total counts
Number of detected genes

It then generates [violin](figures/violin-1.png) and [scatter](figures/scatter-1.png) plots grouped by batch — the same plots you’d show in a lab meeting to justify filtering thresholds.
Cells are filtered by:
Mitochondrial content
Minimum/maximum gene count
Gene detection across cells

This step ensures downstream clustering reflects biology, not damaged cells or ambient RNA.

3. Normalization and highly variable gene selection
The workflow uses:
Library‑size normalization
Log1p transform
HVG selection using the Seurat v3 flavor (batch‑aware)
Only HVGs are retained for dimensionality reduction.
This keeps the pipeline fast and focuses the analysis on biologically informative genes.

4. Optional batch correction
Different experiments need different strategies, so the pipeline supports three modes:

No correction — useful for clean datasets or benchmarking
Harmony — robust for donor‑ or batch‑driven variation
BBKNN — neighbor‑graph correction that plays nicely with Scanpy
Harmony includes careful metadata cleaning to avoid common index‑related crashes.

5. PCA, neighbors, UMAP, and Leiden clustering
Once the representation is chosen (PCA or Harmony), the pipeline:

Builds a neighbor graph
Computes UMAP
Runs Leiden clustering
[UMAPs](figures/umap-1.png) are automatically saved and colored by:
Batch
Sample
Cluster

These are the plots collaborators usually want first.

6. Marker gene identification
The pipeline runs rank_genes_groups (Wilcoxon) and produces:

Summary marker plots
[A small heatmap of top markers per cluster](figures/heatmap-1.png)

This gives a quick sense of cluster identity and helps guide downstream annotation.

7. Final output

This file contains everything needed for downstream work:
Normalized counts
HVGs
PCA / Harmony / BBKNN embeddings
UMAP coordinates
[Leiden clusters](figures/rank_genes_groups_leiden-1.png)
Marker gene statistics

## Citations
If you use this pipeline in a publication, please cite the underlying tools so collaborators and reviewers can trace the computational steps.

Scanpy
[Wolf FA, Angerer P, Theis FJ. SCANPY: large‑scale single‑cell gene expression data analysis. Genome Biology (2018).](https://doi.org/10.1186/s13059-017-1382-0)

AnnData
[Virshup I, Rybakov S, Theis FJ, Angerer P. Anndata: Annotated data.](https://doi.org/10.5281/zenodo.3357167)  

Harmony (if used)
[Korsunsky I et al. Fast, sensitive, and accurate integration of single‑cell data with Harmony. Nature Methods (2019).](https://doi.org/10.1038/s41592-019-0619-0)

BBKNN (if used)
[Polański K et al. BBKNN: Fast batch alignment of single cell transcriptomes. Bioinformatics (2020).](https://doi.org/10.1093/bioinformatics/btz625)

Leiden Clustering
[Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well‑connected communities. Scientific Reports (2019).](https://doi.org/10.1038/s41598-019-41695-z)

UMAP
[McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction.](https://arxiv.org/abs/1802.03426)  

10x Genomics Format
[10x Genomics. Single Cell Gene Expression.](https://www.10xgenomics.com/)  
