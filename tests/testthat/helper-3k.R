# Build a profile sampled every `step_h` hours starting at `start`.
make_profile <- function(values, step_h = 0.5, start = "2024-01-01 16:00:00") {
  data.frame(
    datetime = as.POSIXct(start, tz = "UTC") + (seq_along(values) - 1) * step_h * 3600,
    melatonin = values
  )
}

# Run calculate_dlmo(method = "3k") quietly (no plot, messages and warnings suppressed).
quiet_3k <- function(...) {
  suppressMessages(suppressWarnings(calculate_dlmo(method = "3k", ..., plot = FALSE)))
}

# Run calculate_dlmo(method = "3k") quietly and return the "no DLMO" error condition.
fail_3k <- function(...) {
  tryCatch(quiet_3k(...), dlmo_3k_error = function(e) e)
}
