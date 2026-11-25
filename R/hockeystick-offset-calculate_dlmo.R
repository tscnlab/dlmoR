#' @section Hockey-stick DLMOoff (`method = "hockeystick"`, `dlmo_off = TRUE`):
#' The dim-light melatonin offset is the end of the nocturnal melatonin
#' secretion. Following Danilenko et al. (2014), the profile is mirrored in
#' time so that the falling limb becomes a rise, and the mirrored profile is
#' fitted exactly like a DLMO with the hockey-stick method; the inflection point
#' is then mapped back to real time. Before mirroring, the falling limb of one
#' offset is cut out of the profile with [find_offset_limb()]: from the
#' nocturnal peak to the trough before the next rise. Without this cut, samples
#' of other days (e.g. the first morning of a 40-h profile) are left on the
#' mirrored profile after the rise and can be taken into the base segment.
#'
#' **Arguments** (passed through `...`):
#' - `threshold`, `interval_limit`, `fine_flag`: as for the onset (see the
#'   hockey-stick section above); they are applied to the mirrored limb.
#' - `episode`: which offset to use if the profile has several: `"last"`
#'   (default), `"first"`, `"start"` (a decline already under way at the first
#'   sample) or an integer `k`. See [find_offset_limb()].
#' - `min_run`: minimum number of consecutive samples above the threshold for an
#'   episode or a rise (default `2`). See [find_offset_limb()].
#'
#' **Value**: a list with the same structure as the onset result, on the real
#' time axis, with "ascending" replaced by "descending":
#' - **`prof`**: the whole input profile with `datetime`, `melatonin`, `time`,
#'   `in_limb` (1 = part of the falling limb that was fitted), `base`,
#'   `descending` and `intermediate` (0/1 segment flags).
#' - **`limb`**: the `info` list of [find_offset_limb()] (`start`, `end`,
#'   `crossing`, `n`, ...).
#' - **`prl`**, **`roi`**: parallelogram and region of interest, in real time.
#' - **`threshold`**, **`interval_limit`**, **`episode`**: settings used.
#' - **`ip`**: inflection point(s) and fit parameters, `inflection_point_coarse`,
#'   `base_params_coarse`, `descending_params_coarse`, `grid_big`, `res_big`
#'   (and the `_fine` / `grid_small` / `res_small` counterparts if
#'   `fine_flag = TRUE`), in decimal hours relative to the date of
#'   `ip$datetime_ref` (the DLMOoff).
#' - **`dlmoff`**: list with the `coarse` and (if `fine_flag = TRUE`) `fine`
#'   DLMOoff; each holds `time` (hh:mm:ss), `datetime` (POSIXct; profiles can
#'   span several days), `fit_melatonin` and `fit_lines` (`base` and
#'   `descending`, with x in decimal hours since 00:00 of the DLMOoff date).
#' - **`mirrored`**: the mirrored limb (`profile`), the ends of the limb
#'   (`t_first`, `t_last`) and the onset fit of the mirrored limb (`fit`, without
#'   plots), for checking.
#' - **`dlmoffplotcoarse`**, **`dlmoffplotfine`**: `ggplot` objects made with
#'   [plot_profile_dlmoff()], in the same format as the DLMO plots.
#'
#' The calculation stops with an error when no falling limb is found (the
#' reason is given) or when the hockey-stick fit of the mirrored limb fails.
#'
#' @references
#' Danilenko KV, Verevkin EG, Antyufeev VS, Wirz-Justice A, Cajochen C (2014).
#' The hockey-stick method to estimate evening dim light melatonin onset (DLMO)
#' in humans. *Chronobiol Int* 31(3):349-355.
#'
#' @examples
#'
#' # ---- Hockey-stick DLMOoff ----
#' # The sample profile starts with the morning decline of melatonin
#' filename <- system.file("extdata/sample_melatonin_profile.csv", package = "dlmoR")
#' off <- calculate_dlmo(file_path = filename, dlmo_off = TRUE, threshold = 2.3,
#'                       fine_flag = FALSE)
#' off$dlmoff$coarse$datetime
#' off$dlmoffplotcoarse
#'
#' @rdname calculate_dlmo
#' @order 4
#' @usage NULL
.calculate_dlmoff_hockeystick <- function(data = NULL, file_path = NULL, threshold = 2.3,
                                          interval_limit = lubridate::hours(2),
                                          fine_flag = TRUE, episode = "last", min_run = 2) {
  if (is.null(data) && is.null(file_path)) {
    stop("You must provide either `data` or `file_path`.")
  }
  if (!is.null(file_path)) {
    message("Loading data from file: ", file_path)
    data <- .read_melatonin_data(file_path)
  }
  data <- validate_df_structure(data)
  if (!"datetime" %in% names(data)) {
    stop("DLMOoff needs a `datetime` column (POSIXct).", call. = FALSE)
  }
  data <- data[order(data$datetime), , drop = FALSE]

  # 1. cut out the falling limb
  lb <- find_offset_limb(data, threshold = threshold, episode = episode, min_run = min_run)
  if (!isTRUE(lb$info$ok)) {
    stop("No falling limb found for DLMOoff at threshold = ", threshold, ": ",
         lb$info$reason, ".", call. = FALSE)
  }

  # 2. mirror the limb and fit it like a DLMO
  mir <- .mirror_profile(lb$limb[, c("datetime", "melatonin")])
  t_first <- attr(mir, "t_first")
  t_last <- attr(mir, "t_last")
  fit <- tryCatch(
    suppressMessages(.calculate_dlmo_hockeystick(
      data = as.data.frame(mir), threshold = threshold,
      interval_limit = interval_limit, fine_flag = fine_flag)),
    error = function(e) {
      stop("Hockey-stick fit of the mirrored falling limb (",
           format(t_first, "%Y-%m-%d %H:%M"), " to ", format(t_last, "%Y-%m-%d %H:%M"),
           ") failed: ", conditionMessage(e), call. = FALSE)
    })
  fit$dlmoplotcoarse <- NULL
  fit$dlmoplotfine <- NULL

  # 3. back to real time
  rt <- .offset_to_real_time(fit, data, t_first, t_last, fine_flag)
  res <- list(
    prof = rt$prof,
    limb = lb$info,
    prl = rt$prl,
    roi = rt$roi,
    threshold = threshold,
    interval_limit = fit$interval_limit,
    episode = episode,
    ip = rt$ip,
    dlmoff = rt$dlmoff,
    mirrored = list(profile = mir, t_first = t_first, t_last = t_last, fit = fit)
  )
  # list() keeps a NULL fine plot as an element, as in the onset result
  c(res, list(
    dlmoffplotcoarse = plot_profile_dlmoff(res, plot_coarse = TRUE),
    dlmoffplotfine = if (!is.null(res$dlmoff$fine)) plot_profile_dlmoff(res, plot_coarse = FALSE)
  ))
}
