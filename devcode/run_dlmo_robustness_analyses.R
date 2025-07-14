# -------------------------------------------------------------------------
# Wrapper Script: Run DLMO Robustness Analyses
#
# This script sets up the environment, loads the dlmoR package, and
# sequentially runs two analysis scripts:
#   1. deletion_analysis_multiple_pll.R
#   2. resampling_analysis_pll.R
#
# -------------------------------------------------------------------------

# ---- Load the dlmoR package (assumes installed or in dev environment)
message("Loading dlmoR package...")
if (!requireNamespace("dlmoR", quietly = TRUE)) {
  stop("The dlmoR package must be installed or available in the environment.")
}
library(dlmoR)


# ---- Run deletion analysis script ----
message("Starting deletion analysis...")
source("deletion_analysis_multiple_pll.R")
message("Deletion analysis complete.")

# ---- Run resampling analysis script ----
message("Starting resampling analysis...")
source("resampling_analysis_pll.R")
message("Resampling analysis complete.")

# ---- 5. Done ----
message("All DLMO analyses completed successfully.")
