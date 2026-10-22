# ---- Proposed patch: R/validate_melatonin_profile.R --------------------------
# 1) check_data_length(): return TRUE/FALSE indicating
#    whether the profile meets the minimum data-point requirement.
#
# 2) check_rise_threshold(): determine whether the profile contains a valid
#    rise above the threshold, and return TRUE/FALSE.
#    - Identify consecutive runs of points above the threshold. A profile
#      has a rise if at least one run:

#      (a) contains at least min_increase_points (default 2) or more consecutive points, or
#      (b) extends to the end of the profile (no point after it falls back
#          below the threshold).
#    - A single point that rises above the threshold and then falls back
#      below it (a "blip") does not count as a rise. Multiple blips do
#      not add up to a rise either. If the profile contains only isolated blips,
#      it is classified as having "no dynamic part", regardless of the
#      number of blips.
#    - Rationale: hockeystickexe returns inconsistent results for two separate blips
#      (it fits the earlier blip, the later blip, or reports "no dynamic part"
#      depending on the profile). This rule provides one consistent answer.
#      The choice is open for discussion.
#
# 3) preprocess_profile(): run both validation checks on the original profile,
#    before trimming terminal points below the threshold. If either check fails,
#    stop with the messages described in the paper:
#    "insufficient number of data points", or "no dynamic part".

check_data_length <- function(profile, min_points = 3) {
  length(profile) >= min_points
}

check_rise_threshold <- function(profile, threshold = 2.3, min_increase_points = 2) {
  above <- profile$melatonin >= threshold
  if (!any(above, na.rm = TRUE)) return(FALSE)
  runs   <- rle(above)
  ends   <- cumsum(runs$lengths)                  # last row index of each run
  is_up  <- runs$values
  at_end <- ends == length(above)                 # run reaches the end of the profile
  any(is_up & (runs$lengths >= min_increase_points | at_end))
}

preprocess_profile <- function(profile, threshold = 2.3, min_increase_points = 2, min_points = 3) {
  if (!check_data_length(profile$melatonin, min_points)) {
    stop("insufficient number of data points", call. = FALSE)
  }
  if (!check_rise_threshold(profile, threshold = threshold, min_increase_points = min_increase_points)) {
    stop("no dynamic part", call. = FALSE)
  }
  trim_end_below_threshold(profile, threshold = threshold)
}