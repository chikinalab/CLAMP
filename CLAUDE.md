# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

PLIER2 is an R/Bioconductor package for Pathway-Level Information ExtractoR 2, designed for prior-informed latent variable decomposition of high-dimensional transcriptomic data. It integrates curated gene sets to learn biologically interpretable latent variables and supports file-backed matrices for large datasets.

## Development Setup

### Local Development Environment
Use the conda environment for development:
```bash
conda env create -f envs/plier2.yaml
conda activate plier2
```

Then install the local package in R:
```r
library(remotes)
REPO_PATH <- "~/path/to/PLIER2"  # adjust to your path
remotes::install_local(REPO_PATH, force = TRUE, dependencies = FALSE)
library(PLIER2)
```

## Common Commands

### Package Development
- **Build package**: `R CMD build .` (requires R in PATH or use from conda environment)
- **Check package**: `R CMD check PLIER2_*.tar.gz`
- **Install locally**: `R CMD INSTALL .` or use `remotes::install_local()` from R
- **Run tests**: `testthat::test_check("PLIER2")` or `devtools::test()`
- **Generate documentation**: `devtools::document()` or `roxygen2::roxygenise()`

### Testing
- All tests are in `tests/testthat/` directory
- Run specific test file: `testthat::test_file("tests/testthat/test-filename.R")`
- Test runner is configured in `tests/testthat.R`

### Documentation
- Package website built with pkgdown: `pkgdown::build_site()`
- Vignettes are in `vignettes/` directory
- Main vignette: `vignettes/get_started.Rmd`

## Package Architecture

### Core Functions Structure
The package is organized into two main R files:

**`R/utils.R`** - Data preprocessing and utility functions:
- `preprocessPLIER2()` / `preprocessPLIER2FBM()` - Data filtering and preprocessing
- `zscorePLIER2()` / `zscorePLIER2FBM()` - Z-score normalization  
- `cpmPLIER2()` / `cpmPLIER2FBM()` - Counts per million normalization
- `getGMT()` - Download and parse GMT pathway files
- `gmtListToSparseMat()` - Convert GMT to sparse matrix format

**`R/solvers.R`** - Core PLIER algorithms:
- `PLIERbase()` - Base PLIER decomposition algorithm
- `PLIERfull()` - Full PLIER workflow with pathway integration
- `projectPLIER()` - Project new data onto existing PLIER model
- `getChat()` - Generate pathway matrix (Chat)
- `getMatchedPathwayMat()` - Match pathways to gene features
- `num.pc()` - Determine optimal number of principal components

### FBM Support
The package extensively supports Filebacked Big Matrix (FBM) objects from the `bigstatsr` package for handling large datasets:
- Most functions have both regular matrix and FBM variants (suffix `FBM`)
- FBM functions include parallel processing support via `ncores` parameter
- Use `mat_mult()` helper for matrix multiplication that handles both regular and FBM objects

### Key Dependencies
- **Core computation**: `bigstatsr`, `Matrix`, `rsvd`, `irlba`, `glmnet`
- **Data manipulation**: `dplyr`, `data.table`
- **Visualization**: `ggplot2`, `ggrepel`
- **Bioconductor**: Various annotation packages (suggested)
- **File I/O**: `hdf5r` (suggested for HDF5 support)

### Data Objects
- `dataWholeBlood` - Example whole blood gene expression dataset
- `majorCellTypes` - Major cell type markers for deconvolution

## Development Notes

### Package Standards
- This is a Bioconductor package following Bioconductor guidelines
- Uses roxygen2 for documentation (RoxygenNote: 7.3.2)
- testthat edition 3 for testing
- Supports both in-memory and file-backed matrix operations

### Large Dataset Focus  
The package is specifically designed to handle large transcriptomic datasets (tens of thousands of samples) through:
- File-backed matrix support via `bigstatsr`
- Chunked processing for memory efficiency
- Parallel computation support
- Integration with public resources like recount3 and ARCHS4