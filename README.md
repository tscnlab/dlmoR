<img src="https://github.com/tscnlab/Templates/blob/main/logo/logo_with_text-01.png" width="400"/>

dlmoR: Dim-Light Melatonin Onset Estimation
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

## **`dlmoR`**

**`dlmoR`** is an R package that implements the hockey-stick method
(Danilenko et al., 2014) for estimating dim light melatonin onset (DLMO)
—a key circadian phase marker in chronobiology and sleep research.

The hockey-stick algorithm models melatonin rise as a piecewise
linear-parabolic curve, with the DLMO time point determined as the
inflection point of the piecewise curve that best fits the melatonin
profile, in a least-squares sense. This approach provides a more
objective and robust estimate of DLMO compared to traditional
threshold-based methods, which can be limited by variability in baseline
melatonin levels and subjectivity in visual estimation (Benloucif et
al., 2008; Kennaway, 2023; Glacet et al., 2023).

A Windows-based executable of this algorithm was previously released
(Danilenko & Verevkin, 2020), but its closed-source format limits
flexibility. The original software does not allow modification or
inspection of the underlying algorithm, lacks an API for integration
into analytical workflows, and requires manual operation, making batch
processing inefficient.

By bringing this method into the R-programming environment, **`dlmoR`**
provides an open-source, transparent, and scriptable alternative. It
enables reproducible and batchable DLMO estimation, allowing users to
efficiently analyze multiple melatonin time-series within
high-throughput workflows.

## Installation

You can install the latest development version of **`dlmoR`** from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("tscnlab/dlmoR")
```

## Example

This is a simple example that illustrates how to use **`dlmoR`**:

``` r
library(dlmoR)

# Load the sample melatonin profile data included in the package
filename <- system.file("extdata/sample_melatonin_profile.csv", package = "dlmoR")
# Calculate the DLMO using the sample data and a threshold of 3 pg/mL (default = 2.3 pg/mL)
sample_dlmo <- calculate_dlmo(file_path = filename, threshold = 3, fine_flag = TRUE)
```

``` r
# Load the package
library(dlmoR)

# Load sample data (provided with the package)
filename <- system.file("extdata/sample_melatonin_profile.csv", package = "dlmoR")

# Compute DLMO with a threshold of 3 pg/mL
dlmo_result <- calculate_dlmo(file_path = filename, threshold = 3, fine_flag = TRUE)

# View estimated DLMO inflection point (decimal-hours) and melatonin concentration (pg/mL) at this time
dlmo_result$ip$inflection_point
#        x     y
36   20.26667 0.1

# View estimated DLMO timestamp
dlmo_result$dlmo_time
20:16:00

# Plot the melatonin profile with DLMO detection
print(dlmo_result$dlmoplotcoarse)
print(dlmo_result$dlmoplotfine)
```

**Dim-Light Melatonin Onset (DLMO) Plot (coarse grid search view)**  
![DLMO Coarse Plot](man/figures/dlmo_example_coarse_plot.png)

**Dim-Light Melatonin Onset (reduced fine grid search view)**  
![DLMO Fine Plot](man/figures/dlmo_example_fine_plot.png)

## **Expected Runtime**

The runtime of `calculate_dlmo()` depends on the selected search method:

- **Coarse- and fine-grid DLMO search (`fine_flag = TRUE`, default)**:
  ~10 minutes per melatonin profile.  
- **Coarse-grid DLMO search only (`fine_flag = FALSE`)**: ~1 minute per
  melatonin profile.

For **batch processing**, runtime scales approximately linearly with the
number of profiles analyzed.
