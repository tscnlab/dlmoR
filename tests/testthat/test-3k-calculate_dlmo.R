# ---- building blocks ------------------------------------------------------------

test_that("compute_3k_threshold uses mean + k_sd * SD of the first n samples", {
  thr <- compute_3k_threshold(c(1.0, 1.2, 1.4, 5, 9), n_baseline = 3, k_sd = 2)
  expect_equal(thr$baseline_mean, 1.2)
  expect_equal(thr$baseline_sd, 0.2)
  expect_equal(thr$threshold, 1.6)
  expect_false(thr$floored)

  floored <- compute_3k_threshold(c(0, 0, 0, 1, 5), min_threshold = 0.9)
  expect_equal(floored$threshold_raw, 0)
  expect_equal(floored$threshold, 0.9)
  expect_true(floored$floored)
})

test_that("find_threshold_intervals interpolates crossings and marks open ends", {
  iv <- find_threshold_intervals(0:5, c(1, 1, 3, 4, 1, 5), threshold = 2)
  expect_equal(nrow(iv), 2)
  expect_equal(iv$start_h, c(1.5, 4.25))
  expect_equal(iv$end_h[1], 3 + 2 / 3)
  expect_equal(iv$open_end, c(FALSE, TRUE))
  expect_equal(nrow(find_threshold_intervals(0:3, c(1, 1, 1, 1), threshold = 2)), 0)
  # values equal to the threshold do not count as above
  expect_equal(nrow(find_threshold_intervals(0:3, c(2, 2, 2, 2), threshold = 2)), 0)
  # a sample exactly on the threshold ends the interval; the next one starts there
  tie <- find_threshold_intervals(0:6, c(1, 1, 3, 2, 4, 5, 6), threshold = 2)
  expect_equal(tie$start_h, c(1.5, 3))
  expect_equal(tie$end_h, c(3, 6))
})

test_that("select_dlmo_3k applies the duration rule and reports the status", {
  iv <- find_threshold_intervals(0:8, c(1, 1, 1, 4, 1, 3, 5, 6, 7), threshold = 2)
  sel <- select_dlmo_3k(iv, min_duration = 2)
  expect_equal(sel$status, "ok")
  expect_equal(sel$interval, 2)
  expect_equal(sel$n_transient, 1L)
  expect_equal(select_dlmo_3k(iv, min_duration = 0)$interval, 1)
  expect_equal(select_dlmo_3k(iv[0, ], 2)$status, "no_crossing")
})


# ---- 3k method: normal profiles -------------------------------------------

test_that("a clean rise gives the interpolated crossing time", {
  p <- make_profile(c(1.0, 1.2, 1.4, 1.3, 1.5, 3.0, 6, 10, 14, 18, 20))
  res <- quiet_3k(data = p)
  expect_s3_class(res, "dlmo_3k")
  expect_equal(res$status, "ok")
  expect_equal(res$threshold, 1.6)
  # crossing between 18:00 (1.5) and 18:30 (3.0): 18:00 + (0.1 / 1.5) * 30 min
  expect_equal(format(res$dlmo_datetime, "%H:%M:%S"), "18:02:00")
  expect_equal(as.character(res$dlmo_time), "18:02:00")
  expect_false(any(res$flags))
  expect_equal(res$prof$baseline, c(1, 1, 1, rep(0, 8)))
})

test_that("short excursions before the sustained rise are ignored and noted", {
  p <- make_profile(c(1, 1.2, 1.4, 1.3, 5, 1.2, 1.3, 1.4, 3, 6, 9, 12, 15, 18))
  msgs <- capture_messages(res <- suppressWarnings(calculate_dlmo(method = "3k", data = p, plot = FALSE)))
  expect_match(msgs, "shorter excursion", all = FALSE)
  expect_equal(format(res$dlmo_datetime, "%H:%M:%S"), "19:33:45")
  expect_true(res$flags[["transient_excursion"]])
  expect_equal(sum(res$intervals$selected), 1)
  expect_equal(which(res$intervals$selected), 2)
})

test_that("min_duration = 0 returns the first crossing", {
  p <- make_profile(c(1, 1.2, 1.4, 1.3, 5, 1.2, 1.3, 1.4, 3, 6, 9, 12, 15, 18))
  res <- quiet_3k(data = p, min_duration = 0)
  expect_equal(res$selected_interval, 1)
})

test_that("profiles crossing midnight keep the correct date", {
  p <- make_profile(c(1.0, 1.2, 1.4, 1.3, 1.5, 3.0, 6, 10, 14, 18, 20), start = "2024-01-01 22:00:00")
  res <- quiet_3k(data = p)
  expect_equal(format(res$dlmo_datetime, "%Y-%m-%d %H:%M:%S"), "2024-01-02 00:02:00")
})

test_that("clock times without dates are unwrapped at midnight", {
  p <- make_profile(c(1.0, 1.2, 1.4, 1.3, 1.5, 3.0, 6, 10, 14, 18, 20), start = "2024-01-01 22:00:00")
  only_time <- data.frame(time = hms::as_hms(p$datetime), melatonin = p$melatonin)
  res <- quiet_3k(data = only_time)
  expect_equal(as.character(res$dlmo_time), "00:02:00")
  expect_equal(res$status, "ok")
})

test_that("the sample file can be analysed from file_path", {
  f <- system.file("extdata/sample_melatonin_profile.csv", package = "dlmoR")
  # the file starts during the morning decline: inflated threshold, rise cut off
  e <- fail_3k(file_path = f)
  expect_s3_class(e, "dlmo_3k_error")
  expect_match(conditionMessage(e), "profile ends")
  expect_equal(e$result$status, "truncated")

  # starting in the afternoon trough gives a DLMO
  d <- readr::read_delim(f, delim = ";", show_col_types = FALSE)
  res <- quiet_3k(data = d[d$datetime >= as.POSIXct("2024-04-17 12:00:00", tz = "UTC"), ])
  expect_equal(res$status, "ok")
  expect_equal(format(res$dlmo_datetime, "%H:%M:%S"), "20:27:51")
})


# ---- 3k method: no DLMO -------------------------------------------------------

test_that("no crossing stops with an error that carries the partial result", {
  p <- make_profile(c(1, 1.2, 1.4, 1.1, 1.3, 1.2, 1.0, 1.1))
  expect_error(quiet_3k(data = p), "never rises", class = "dlmo_3k_error")
  e <- fail_3k(data = p)
  expect_equal(e$result$status, "no_crossing")
  expect_true(is.na(e$result$dlmo_datetime))
  expect_equal(e$result$threshold, compute_3k_threshold(p$melatonin)$threshold)
})

test_that("a rise that does not last min_duration stops as not_sustained", {
  p <- make_profile(c(1, 1.2, 1.4, 1.3, 3, 1.2, 1.3, 1.1, 1.2, 1.3))
  expect_error(quiet_3k(data = p), "never stays above", class = "dlmo_3k_error")
  expect_equal(fail_3k(data = p)$result$status, "not_sustained")
})

test_that("a rise cut off by the end of the data stops as truncated", {
  p <- make_profile(c(1, 1.2, 1.4, 1.3, 1.2, 1.3, 1.4, 1.5, 3, 6))
  expect_error(quiet_3k(data = p), "profile ends", class = "dlmo_3k_error")
  e <- fail_3k(data = p)
  expect_equal(e$result$status, "truncated")
  expect_true(tail(e$result$intervals$open_end, 1))
})

test_that("the no-DLMO error carries a plot when plot = TRUE", {
  p <- make_profile(c(1, 1.2, 1.4, 1.1, 1.3, 1.2, 1.0, 1.1))
  e <- tryCatch(suppressMessages(calculate_dlmo(method = "3k", data = p)), dlmo_3k_error = function(e) e)
  expect_s3_class(e$result$dlmoplot, "ggplot")
  expect_null(fail_3k(data = p)$result$dlmoplot)
})


# ---- 3k method: warnings about the baseline ------------------------------------

test_that("a baseline with SD 0 is flagged and does not count as the rise", {
  p <- make_profile(c(0, 0, 0, 0, 0.5, 3, 6, 9, 12, 15))
  res <- quiet_3k(data = p)
  expect_true(res$flags[["baseline_sd_zero"]])
  expect_equal(res$threshold, 0)
  # leaves 0 between 17:30 and 18:00, not at the first sample
  expect_equal(format(res$dlmo_datetime, "%H:%M:%S"), "17:30:00")

  floored <- quiet_3k(data = p, min_threshold = 0.9)
  expect_true(floored$flags[["threshold_floored"]])
  expect_false(floored$flags[["baseline_sd_zero"]])
  expect_equal(floored$threshold, 0.9)
  expect_equal(format(floored$dlmo_datetime, "%H:%M:%S"), "18:04:48")
})

test_that("baseline values at the assay floor are flagged", {
  p <- make_profile(c(0, 0.1, 0, 0.1, 0.4, 3, 6, 9, 12, 15))
  expect_true(quiet_3k(data = p)$flags[["baseline_at_floor"]])
  p2 <- make_profile(c(0.2, 0.25, 0.3, 0.3, 0.4, 3, 6, 9, 12, 15))
  expect_false(quiet_3k(data = p2)$flags[["baseline_at_floor"]])
  expect_true(quiet_3k(data = p2, lod = 0.3)$flags[["baseline_at_floor"]])
})

test_that("a long sampling gap around the DLMO is flagged", {
  p <- data.frame(
    datetime = as.POSIXct("2024-01-01 16:00:00", tz = "UTC") +
      c(0, 0.5, 1, 1.5, 2, 5, 5.5, 6, 6.5, 7) * 3600,
    melatonin = c(1, 1.2, 1.4, 1.3, 1.2, 8, 10, 12, 14, 16))
  res <- quiet_3k(data = p)
  expect_equal(res$status, "ok")
  expect_true(res$flags[["gap_near_dlmo"]])
  expect_false(quiet_3k(data = p, max_gap = 4)$flags[["gap_near_dlmo"]])
})


# ---- input handling ---------------------------------------------------------------------

test_that("missing values are removed with a warning", {
  p <- make_profile(c(1.0, 1.2, 1.4, NA, 1.5, 3.0, 6, 10, 14, 18, 20))
  expect_warning(res <- suppressMessages(calculate_dlmo(method = "3k", data = p, plot = FALSE)), "missing")
  expect_true(res$flags[["missing_values_removed"]])
  expect_equal(nrow(res$prof), 10)
})

test_that("unsorted rows are sorted and give the same result", {
  p <- make_profile(c(1.0, 1.2, 1.4, 1.3, 1.5, 3.0, 6, 10, 14, 18, 20))
  shuffled <- p[c(5, 1, 11, 3, 2, 4, 7, 6, 9, 8, 10), ]
  res <- quiet_3k(data = shuffled)
  expect_true(res$flags[["unsorted_input"]])
  expect_equal(res$dlmo_datetime, quiet_3k(data = p)$dlmo_datetime)
})

test_that("invalid input stops with an informative error", {
  p <- make_profile(c(1.0, 1.2, 1.4, 1.3, 1.5, 3.0, 6, 10))
  expect_error(calculate_dlmo(method = "3k"), "either `data` or `file_path`")
  expect_error(quiet_3k(data = p, n_baseline = 1), "n_baseline")
  expect_error(quiet_3k(data = p, n_baseline = 2.5), "whole number")
  expect_error(quiet_3k(data = p, k_sd = -1), "k_sd")
  expect_error(quiet_3k(data = p, min_duration = -1), "min_duration")
  expect_error(quiet_3k(data = p[1:4, ]), "At least 5")
  expect_error(quiet_3k(data = rbind(p, p[3, ])), "Duplicated")
  neg <- p; neg$melatonin[2] <- -1
  expect_error(quiet_3k(data = neg), "negative")
})


# ---- output -------------------------------------------------------------------------------

test_that("the result carries a plot and prints a summary", {
  p <- make_profile(c(1.0, 1.2, 1.4, 1.3, 1.5, 3.0, 6, 10, 14, 18, 20))
  res <- suppressMessages(calculate_dlmo(method = "3k", data = p))
  expect_s3_class(res$dlmoplot, "ggplot")
  expect_s3_class(plot_profile_3k(res), "ggplot")
  expect_output(print(res), "18:02:00")
  no <- fail_3k(data = make_profile(c(1, 1.2, 1.4, 1.1, 1.3, 1.2, 1.0, 1.1)))$result
  expect_s3_class(plot_profile_3k(no), "ggplot")
  expect_output(print(no), "not determined")
})
