# Helpers of the hockey-stick DLMOoff: mirror a profile in time and map a
# hockey-stick fit of the mirrored profile back to real time.
#
# Mirroring: t' = t_first + (t_last - t). The mirrored profile covers the same
# clock interval as the original, so a falling limb becomes a rise that the
# onset pipeline can fit. The onset t'_on found on it maps back to
# DLMOoff = t_first + t_last - t'_on.
#
# dlmoR keeps times in decimal hours relative to the date of a reference time
# (posixct_to_decimal()). A mirrored decimal hour x' and the real decimal hour
# x of the same sample are related by x' = K - x for a constant K, so every
# fitted function is transformed exactly:
#   line      y = m'x' + b'          ->  y = (-m')x + (m'K + b')
#   parabola  y = a'x'^2 + b'x' + c' ->  y = a'x^2 + (-2a'K - b')x + (a'K^2 + b'K + c')


#' Mirror a Melatonin Profile in Time
#'
#' @param data Data frame with a POSIXct `datetime` column.
#' @return The mirrored data frame, sorted by time, with attributes `t_first`
#'   and `t_last` (the ends of the original profile).
#' @keywords internal
#' @noRd
.mirror_profile <- function(data) {
  t_first <- min(data$datetime)
  t_last <- max(data$datetime)
  out <- data
  out$datetime <- .unmirror_time(data$datetime, t_first, t_last)
  out <- out[order(out$datetime), , drop = FALSE]
  attr(out, "t_first") <- t_first
  attr(out, "t_last") <- t_last
  out
}

#' Map Times Between the Real and the Mirrored Axis
#'
#' The mapping is its own inverse, so it is used in both directions.
#' @keywords internal
#' @noRd
.unmirror_time <- function(t, t_first, t_last) {
  t_first + as.numeric(difftime(t_last, t, units = "secs"))
}

#' Transform Fit Parameters from the Mirrored to the Real Axis
#'
#' @param base_params Slope of the base line on the mirrored axis.
#' @param asc_params Slope (linear) or list `a`, `b`, `c` (parabolic) of the
#'   rise fit on the mirrored axis.
#' @param K Constant with x' = K - x.
#' @return list(base = slope, descending = slope or list(a, b, c)) on the real axis.
#' @keywords internal
#' @noRd
.mirror_fit_params <- function(base_params, asc_params, K) {
  desc <- if (is.list(asc_params) && all(c("a", "b", "c") %in% names(asc_params))) {
    a <- unname(asc_params$a); b <- unname(asc_params$b); c <- unname(asc_params$c)
    list(a = a, b = -2 * a * K - b, c = a * K^2 + b * K + c)
  } else {
    -unname(asc_params)
  }
  list(base = -unname(base_params), descending = desc)
}

#' Map a Hockey-Stick Fit of a Mirrored Profile Back to Real Time
#'
#' @param fit Result of the hockey-stick onset calculation on the mirrored
#'   profile.
#' @param data The whole (validated) input profile, in real time.
#' @param t_first,t_last Ends of the profile that was mirrored (the limb).
#' @param fine_flag Whether a fine search was run.
#' @return A list with `prof`, `prl`, `roi`, `ip` and `dlmoff` on the real
#'   time axis (see `calculate_dlmo()`, section on the hockey-stick DLMOoff).
#' @keywords internal
#' @noRd
.offset_to_real_time <- function(fit, data, t_first, t_last, fine_flag) {
  unm <- function(t) .unmirror_time(t, t_first, t_last)
  ipm <- fit$ip
  levels <- c("coarse", if (fine_flag && !is.null(ipm$inflection_point_fine)) "fine")
  sfx <- c(coarse = "coarse", fine = "fine")
  grid_name <- c(coarse = "grid_big", fine = "grid_small")
  res_name <- c(coarse = "res_big", fine = "res_small")

  # DLMOoff of the final search level; its date is the reference of the real axis
  final <- levels[length(levels)]
  ip_final <- ipm[[paste0("inflection_point_", final)]]
  datetime_ref <- unm(decimal_to_posixct(ip_final$x, ipm$datetime_ref))

  # K with x' = K - x, computed with dlmoR's own conversion
  k_const <- function(ref_m) {
    posixct_to_decimal(unm(t_first), ref_m) + posixct_to_decimal(t_first, datetime_ref)
  }
  K <- k_const(ipm$datetime_ref)

  ip <- list()
  dlmoff <- list(coarse = NULL, fine = NULL)
  for (lv in levels) {
    s <- sfx[[lv]]
    ip_m <- ipm[[paste0("inflection_point_", s)]]
    par <- .mirror_fit_params(ipm[[paste0("base_params_", s)]],
                              ipm[[paste0("ascending_params_", s)]], K)
    ip_r <- ip_m
    ip_r$x <- K - ip_m$x
    grid_m <- ipm[[grid_name[[lv]]]]
    grid_r <- grid_m
    if (!is.null(grid_m)) grid_r$x <- K - grid_m$x
    ip[[paste0("inflection_point_", s)]] <- ip_r
    ip[[paste0("base_params_", s)]] <- par$base
    ip[[paste0("descending_params_", s)]] <- par$descending
    ip[[grid_name[[lv]]]] <- grid_r
    ip[[res_name[[lv]]]] <- ipm[[res_name[[lv]]]]

    t_off <- decimal_to_posixct(ip_r$x, datetime_ref)
    fl <- process_fits(list(inflection_point = ip_r, base_params = par$base,
                            ascending_params = par$descending),
                       list(time = hms::as_hms(t_off), fit_melatonin = ip_r$y))
    dlmoff[[lv]] <- list(time = fl$time, datetime = t_off, fit_melatonin = fl$fit_melatonin,
                         fit_lines = list(base = fl$fit_lines$base,
                                          descending = fl$fit_lines$ascending))
  }
  ip$datetime_ref <- datetime_ref

  # segment flags of the mirrored fit, on the whole input profile
  pm <- fit$prof
  key_m <- round(as.numeric(unm(pm$datetime)))
  key_r <- round(as.numeric(data$datetime))
  idx <- match(key_r, key_m)
  flag <- function(col) {
    v <- if (col %in% names(pm)) pm[[col]][idx] else rep(NA_real_, length(idx))
    v[is.na(v)] <- 0
    v
  }
  key_limb <- round(as.numeric(c(t_first, t_last)))
  prof <- tidyr::as_tibble(data)
  if (!"time" %in% names(prof)) prof$time <- hms::as_hms(prof$datetime)
  prof$in_limb <- as.numeric(key_r >= key_limb[1] & key_r <= key_limb[2])
  prof$base <- flag("base")
  prof$descending <- flag("ascending")
  prof$intermediate <- flag("intermediate")

  # region of interest
  roi <- fit$roi
  roi_r <- list(x_start = unm(roi$x_end), x_end = unm(roi$x_start),
                y_min = roi$y_min, y_max = roi$y_max)

  # parallelogram (corners in decimal hours relative to datetime_ref)
  pll <- fit$prl
  pll_r <- NULL
  if (!is.null(pll) && !is.null(pll$corners)) {
    Kp <- k_const(pm$datetime[3])  # dlmoR's parallelogram reference (see parallelogram_fit())
    pll_r <- pll
    pll_r$corners <- lapply(pll$corners, function(cn) c(Kp - cn[1], cn[2]))
    pll_r$pll_datetime_0 <- unm(pll$pll_datetime_1)
    pll_r$pll_datetime_1 <- unm(pll$pll_datetime_0)
    pll_r$pll_slope <- -pll$pll_slope
    pll_r$datetime_ref <- datetime_ref
  }

  list(prof = prof, prl = pll_r, roi = roi_r, ip = ip, dlmoff = dlmoff)
}
