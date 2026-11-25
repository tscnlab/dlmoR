# Tests for the hockey-stick DLMOoff:
# calculate_dlmo(method = "hockeystick", dlmo_off = TRUE), find_offset_limb()
# and the mapping of the mirrored fit back to real time.

sample_file <- function() {
  system.file("extdata/sample_melatonin_profile.csv", package = "dlmoR")
}

# Night profile: falling limb from a peak at 03:00, trough at 16:00, next rise at 21:00.
offset_profile <- function() {
  make_profile(c(20, 25, 30, 28, 24, 18, 12, 7, 4, 2.5, 1.8, 1.5, 1.4, 1.2, 1.3, 1.1,
                 1.0, 0.9, 1.0, 0.8, 0.7, 0.8, 0.6, 0.9, 0.5, 1.2, 4, 9, 15),
               step_h = 0.75, start = "2024-01-02 00:00:00")
}

test_that("mirroring is its own inverse", {
  p <- offset_profile()
  m <- dlmoR:::.mirror_profile(p)
  back <- dlmoR:::.mirror_profile(m)
  expect_equal(as.numeric(back$datetime), as.numeric(p$datetime))
  expect_equal(back$melatonin, p$melatonin)
  expect_equal(rev(m$melatonin), p$melatonin)
})

test_that("fit parameters are transformed exactly from the mirrored axis", {
  K <- 37.3
  par <- dlmoR:::.mirror_fit_params(0.05, list(a = 0.4, b = -2, c = 3), K)
  x <- c(8, 10.5, 13)
  # the real-time functions at x equal the mirrored functions at x' = K - x
  expect_equal(par$base * x, -0.05 * x)
  expect_equal(par$descending$a * x^2 + par$descending$b * x + par$descending$c,
               0.4 * (K - x)^2 - 2 * (K - x) + 3)
  lin <- dlmoR:::.mirror_fit_params(0.05, 3.2, K)
  expect_equal(lin$descending, -3.2)
})

test_that("find_offset_limb cuts from the nocturnal peak to the trough before the next rise", {
  p <- offset_profile()
  lb <- find_offset_limb(p, threshold = 2.3)
  expect_true(lb$info$ok)
  expect_equal(lb$info$start, p$datetime[3])           # peak 30
  expect_equal(lb$info$crossing, p$datetime[11])       # first sample <= 2.3
  expect_equal(lb$info$end, p$datetime[25])            # trough 0.5 before the rise
  expect_true(lb$info$next_rise_found)
  expect_equal(nrow(lb$limb), 23)
})

test_that("find_offset_limb ignores blips and picks the requested episode", {
  # two offsets: one at the start of the profile, one 24 h later; a blip in between
  v <- c(10, 6, 3, 1, 1, 1, 4, 1, 1, 1, 1, 1, 5, 9, 12, 8, 4, 1, 1, 1)
  p <- make_profile(v, step_h = 1.5)
  first <- find_offset_limb(p, threshold = 2.3, episode = "first")
  last <- find_offset_limb(p, threshold = 2.3, episode = "last")
  start <- find_offset_limb(p, threshold = 2.3, episode = "start")
  expect_equal(first$info$start, p$datetime[1])
  expect_equal(start$info$start, p$datetime[1])
  expect_equal(last$info$start, p$datetime[15])        # peak 12 of the second episode
  expect_false(last$info$next_rise_found)
  expect_equal(find_offset_limb(p, threshold = 2.3, episode = 2)$info$start, p$datetime[15])
  # the blip at sample 7 is not an episode
  expect_match(find_offset_limb(p, threshold = 2.3, episode = 3)$info$reason, "only 2 found")
})

test_that("find_offset_limb reports why no limb was found", {
  low <- make_profile(c(1, 1.2, 1.1, 0.9, 1.0))
  expect_false(find_offset_limb(low, threshold = 2.3)$info$ok)
  expect_match(find_offset_limb(low, threshold = 2.3)$info$reason, "no sustained descent")
  rise <- make_profile(c(1, 1.2, 3, 6, 9))
  expect_match(find_offset_limb(rise, threshold = 2.3, episode = "start")$info$reason,
               "no sustained descent|does not start above")
  expect_error(find_offset_limb(rise, episode = "middle"), "`episode` must be")
  expect_error(find_offset_limb(rise, threshold = "a"), "`threshold` must be")
})

test_that("DLMOoff on the sample profile is consistent", {
  res <- suppressMessages(calculate_dlmo(file_path = sample_file(), dlmo_off = TRUE,
                                         threshold = 2.3, fine_flag = FALSE))
  expect_named(res, c("prof", "limb", "prl", "roi", "threshold", "interval_limit", "episode",
                      "ip", "dlmoff", "mirrored", "dlmoffplotcoarse", "dlmoffplotfine"))
  d <- res$dlmoff$coarse
  # the DLMOoff lies on the falling limb, after the threshold crossing ...
  expect_true(d$datetime >= res$limb$start && d$datetime <= res$limb$end)
  # ... and is the mirrored onset mapped back to real time
  t_on <- decimal_to_posixct(res$mirrored$fit$ip$inflection_point_coarse$x,
                             res$mirrored$fit$ip$datetime_ref)
  t_back <- dlmoR:::.unmirror_time(t_on, res$mirrored$t_first, res$mirrored$t_last)
  expect_lte(abs(as.numeric(difftime(d$datetime, t_back, units = "secs"))), 1)
  # both fitted functions meet at the DLMOoff
  x0 <- res$ip$inflection_point_coarse$x
  fl <- d$fit_lines
  f_desc <- if (fl$descending$type == "parabolic") {
    fl$descending$a * x0^2 + fl$descending$b * x0 + fl$descending$c
  } else fl$descending$m * x0 + fl$descending$b
  expect_equal(fl$base$m * x0 + fl$base$b, d$fit_melatonin, tolerance = 1e-6)
  expect_equal(f_desc, d$fit_melatonin, tolerance = 1e-6)
  # the base line follows the descent (non-increasing in real time is not required),
  # and the descending fit falls at the DLMOoff
  slope_desc <- if (fl$descending$type == "parabolic") 2 * fl$descending$a * x0 + fl$descending$b else fl$descending$m
  expect_lt(slope_desc, 0)
  # segment flags: descending samples come before base samples
  expect_true(max(res$prof$datetime[res$prof$descending == 1]) <
                min(res$prof$datetime[res$prof$base == 1]))
  expect_equal(nrow(res$prof), 19)
  expect_s3_class(res$dlmoffplotcoarse, "ggplot")
  expect_null(res$dlmoffplotfine)
  expect_s3_class(plot_profile_dlmoff(res), "ggplot")
})

test_that("DLMOoff stops with a clear error when no falling limb exists", {
  rise <- make_profile(c(1.0, 1.2, 1.4, 1.3, 1.5, 3.0, 6, 10, 14, 18, 20))
  expect_error(calculate_dlmo(data = rise, dlmo_off = TRUE, fine_flag = FALSE),
               "No falling limb found")
  expect_error(calculate_dlmo(data = rise, dlmo_off = TRUE, n_baseline = 3),
               "not used by method = \"hockeystick\": `n_baseline`")
})
