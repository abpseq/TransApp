# TransAPP RNA-seq Analysis Pipeline

TransAPP is an interactive Shiny application developed for the analysis and visualization of RNA-seq data. It integrates differential expression analysis, quality control, and multiple visualization options for users handling transcriptomic datasets.

## Features

- DE analysis with multiple dataset support
- PCA, heatmaps, and volcano plot visualization
- GUI via R Shiny
- Integrated support for `.rda` datasets

## Installation

1. Install R and RStudio.
2. Install required R packages:
```r
install.packages(c("shiny", "Biobase", "limma", "edgeR", "gplots", "heatmap3", "ggplot2", "shinythemes"))
```

## Launch

Clone the repository and run:
```r
shiny::runApp("app")
```

## Usage

- Upload `.rda` datasets from `example_data/`
- Explore DE genes and visualization tools
- Refer to `docs/` for installation and usage guides

## Example Data

- `GSE51808.rda`, `GSE42589.rda`: Sample expression sets
- `process.eset.R`: Script to pre-process data into ExpressionSet objects

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
