# Spatial transcriptomics analysis (Seurat, R)

This repository contains R scripts for analyzing spatial transcriptomics data using Seurat, including:
- building a Seurat object from spatial transcriptomics input files
- quality control and normalization
- clustering and downstream analyses
- spatial projection/visualization functions to overlay results on tissue images

## Files
- `analysis.R`: end-to-end analysis workflow
- `spatial_plotting_functions.R`: function for spatial visualization
- `R/`: Helper functions sourced by the analysis pipeline
- `results/`: Output directory for generated tables and figures (created at runtime)

## Notes
This repository is intended as a portfolio example to demonstrate R coding skills.
Large raw data files are not included.
