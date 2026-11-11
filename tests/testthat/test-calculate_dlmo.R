# Tests for the calculate_dlmo() dispatcher: method choice, dlmo_off and the
# method-specific arguments passed through `...`.

sample_file <- function() {
  system.file("extdata/sample_melatonin_profile.csv", package = "dlmoR")
}

test_that("the default method is a single registered method", {
  default <- formals(calculate_dlmo)$method
  expect_identical(default, "hockeystick")
  expect_true(default %in% names(dlmoR:::.dlmo_method_registry()))
})

test_that("every registered calculation function takes data and file_path first", {
  for (m in dlmoR:::.dlmo_method_registry()) {
    for (f in Filter(Negate(is.null), m)) {
      expect_equal(names(formals(f))[1:2], c("data", "file_path"))
    }
  }
})

test_that("the default method is the hockey-stick method", {
  res_default <- suppressMessages(
    calculate_dlmo(file_path = sample_file(), threshold = 5, fine_flag = FALSE))
  res_hs <- suppressMessages(
    dlmoR:::.calculate_dlmo_hockeystick(file_path = sample_file(), threshold = 5,
                                        fine_flag = FALSE))
  expect_named(res_default, c("prof", "prl", "roi", "threshold", "interval_limit",
                              "ip", "dlmo", "dlmoplotcoarse", "dlmoplotfine"))
  expect_equal(res_default$dlmo, res_hs$dlmo)
  expect_equal(res_default$threshold, 5)
})

test_that("method = \"3k\" gives the 3k result", {
  p <- make_profile(c(1.0, 1.2, 1.4, 1.3, 1.5, 3.0, 6, 10, 14, 18, 20))
  res <- suppressMessages(calculate_dlmo(data = p, method = "3k", plot = FALSE))
  expect_s3_class(res, "dlmo_3k")
  expect_equal(as.character(res$dlmo_time), "18:02:00")
})

test_that("partial argument names still work, as in a normal call", {
  p <- make_profile(c(1.0, 1.2, 1.4, 1.3, 1.5, 3.0, 6, 10, 14, 18, 20))
  res <- suppressMessages(calculate_dlmo(data = p, method = "3k", n_base = 4, plot = FALSE))
  expect_equal(res$settings$n_baseline, 4)
})

test_that("invalid method, dlmo_off and arguments stop with clear errors", {
  p <- make_profile(c(1.0, 1.2, 1.4, 1.3, 1.5, 3.0, 6, 10, 14, 18, 20))
  expect_error(calculate_dlmo(), "either `data` or `file_path`")
  expect_error(calculate_dlmo(data = p, method = "cosinor"), "Unknown `method`")
  expect_error(calculate_dlmo(data = p, method = c("3k", "hockeystick")), "single string")
  expect_error(calculate_dlmo(data = p, method = 3), "single string")
  expect_error(calculate_dlmo(data = p, method = NA_character_), "single string")
  # old positional call: the third argument used to be `threshold`
  expect_error(calculate_dlmo(p, NULL, 3), "single string")
  expect_error(calculate_dlmo(data = p, dlmo_off = "yes"), "TRUE or FALSE")
  expect_error(calculate_dlmo(data = p, dlmo_off = NA), "TRUE or FALSE")
  expect_error(calculate_dlmo(data = p, method = "3k", threshold = 5),
               "not used by method = \"3k\": `threshold`")
  expect_error(calculate_dlmo(data = p, method = "hockeystick", k_sd = 2),
               "not used by method = \"hockeystick\": `k_sd`")
  expect_error(calculate_dlmo(data = p, NULL, "3k", FALSE, 4), "must be named")
})

test_that("dlmo_off = TRUE is not implemented yet for any method", {
  p <- make_profile(c(1.0, 1.2, 1.4, 1.3, 1.5, 3.0, 6, 10, 14, 18, 20))
  for (m in names(dlmoR:::.dlmo_method_registry())) {
    expect_error(calculate_dlmo(data = p, method = m, dlmo_off = TRUE),
                 "not implemented yet")
  }
})
