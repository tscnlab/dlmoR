# Residual heatmap explorer (Shiny)

This directory contains a Shiny application used to interactively explore the residual surface produced by the `dlmoR` hockey-stick DLMO estimation procedure.

The app visualizes the residual grid returned by `dlmoR::calculate_dlmo()`, identifies low-residual clusters using 4-connected component labeling, and summarizes cluster timing statistics in absolute clock time and relative to the inferred DLMO.

## Input
The app expects a precomputed object, returned by `dlmoR::calculate_dlmo()` (e.g., `dlmo_result`) to exist in the R session. This object must be the output of:

```r
dlmoR::calculate_dlmo(...)
```

## Outputs
During use, the app writes publication-ready SVG figures to the `results/` directory:

- `results/residual_heatmap.svg`
- `results/offset_violin.svg`

These files are overwritten when plots are re-rendered.

## Run
From the root of the repository:

```r
source("code/shiny/residual_explorer/app.R")
```
