#' Find the Falling Limb of a Melatonin Offset
#'
#' Cuts the falling limb of one melatonin offset out of a profile. This is the
#' first step of the hockey-stick DLMOoff calculation
#' (`calculate_dlmo(method = "hockeystick", dlmo_off = TRUE)`): the limb is
#' mirrored in time and fitted like a DLMO. Only melatonin values and the
#' threshold are used (no sleep or wake times).
#'
#' The rules mirror the way the hockey-stick onset defines a rise:
#' 1. **Episode.** A nocturnal episode is a run of at least `min_run`
#'    consecutive samples above `threshold` (or a run that starts at the first
#'    sample) followed by a sample at or below `threshold`. Single-sample
#'    excursions ("blips") are ignored, as for the onset.
#' 2. **Start of the limb** is the maximum of the chosen episode (the
#'    nocturnal peak; the later one if tied), so the rising part of the night is
#'    removed.
#' 3. **End of the limb** is the local minimum reached by walking back from the
#'    first sample of the next sustained rise above `threshold` (a run of at
#'    least `min_run` samples, or one that reaches the last sample) while the
#'    values keep decreasing, so the next rising limb is removed. Without a
#'    later rise the limb runs to the end of the profile.
#'
#' @param data A data frame with columns `datetime` (POSIXct) and `melatonin`
#'   (numeric), sorted or unsorted.
#' @param threshold Numeric melatonin threshold (pg/mL), as for the onset.
#' @param episode Which episode to use when the profile contains several:
#'   `"last"` (default; the most recent offset), `"first"`, `"start"` (only an
#'   episode that begins at the first sample, i.e. a decline already under way
#'   when sampling starts; useful for the first morning of a multi-day
#'   profile), or an integer `k` (the `k`-th episode in time order).
#' @param min_run Integer. Minimum number of consecutive samples above the
#'   threshold for an episode or a later rise (default `2`).
#'
#' @return A list with
#' - **`limb`**: tibble with the samples of the falling limb (`datetime`,
#'   `melatonin` and any other columns of `data`), or `NULL` if no limb was
#'   found;
#' - **`info`**: list with `ok` (logical), `reason` (why no limb was found,
#'   otherwise `NA`), `start`, `end` (POSIXct ends of the limb), `crossing` (first
#'   sample at or below the threshold after the episode), `n`, `n_above`, `peak`,
#'   `trough` and `next_rise_found` (logical).
#'
#' @examples
#' filename <- system.file("extdata/sample_melatonin_profile.csv", package = "dlmoR")
#' sample_profile <- readr::read_delim(filename, delim = ";", show_col_types = FALSE)
#' lb <- find_offset_limb(sample_profile, threshold = 2.3)
#' lb$info[c("start", "end", "n")]
#'
#' @seealso [calculate_dlmo()] (`method = "hockeystick"`, `dlmo_off = TRUE`)
#' @export
find_offset_limb <- function(data, threshold = 2.3, episode = "last", min_run = 2) {
  if (!is.data.frame(data) || !all(c("datetime", "melatonin") %in% names(data))) {
    stop("`data` must be a data frame with columns `datetime` and `melatonin`.", call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1 || is.na(threshold)) {
    stop("`threshold` must be a single number.", call. = FALSE)
  }
  ok_episode <- (is.character(episode) && length(episode) == 1 &&
                   episode %in% c("last", "first", "start")) ||
    (is.numeric(episode) && length(episode) == 1 && !is.na(episode) &&
       episode >= 1 && episode == round(episode))
  if (!ok_episode) {
    stop("`episode` must be \"last\", \"first\", \"start\" or a positive whole number.",
         call. = FALSE)
  }
  if (!is.numeric(min_run) || length(min_run) != 1 || is.na(min_run) || min_run < 1) {
    stop("`min_run` must be a positive whole number.", call. = FALSE)
  }

  data <- tidyr::as_tibble(data)
  data <- data[!is.na(data$datetime) & !is.na(data$melatonin), , drop = FALSE]
  data <- data[order(data$datetime), , drop = FALSE]
  n <- nrow(data)
  mel <- data$melatonin

  no_limb <- function(reason) {
    list(limb = NULL,
         info = list(ok = FALSE, reason = reason, start = NA, end = NA, crossing = NA,
                     n = 0L, n_above = 0L, peak = NA_real_, trough = NA_real_,
                     next_rise_found = NA))
  }
  if (n < 3) return(no_limb("insufficient number of data points"))

  above <- mel > threshold
  r <- rle(above)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  nr <- length(r$values)

  # 1. candidate episodes: a sustained run above threshold (or one touching the
  #    first sample) that is followed by a sample at or below threshold
  cand <- which(r$values & c(!r$values[-1], FALSE) &
                  (r$lengths >= min_run | starts == 1))
  if (length(cand) == 0) return(no_limb("no sustained descent below threshold"))

  if (identical(episode, "start")) {
    cand <- cand[starts[cand] == 1]
    if (length(cand) == 0) {
      return(no_limb("profile does not start above threshold (no decline at profile start)"))
    }
    k_ep <- cand[1]
  } else if (identical(episode, "last")) {
    k_ep <- max(cand)
  } else if (identical(episode, "first")) {
    k_ep <- min(cand)
  } else {
    if (episode > length(cand)) {
      return(no_limb(paste0("episode ", episode, " requested, but only ",
                            length(cand), " found")))
    }
    k_ep <- cand[episode]
  }
  i_cross <- ends[k_ep] + 1

  # 2. start of the limb: peak of the episode (ties: the later sample)
  ep_idx <- starts[k_ep]:ends[k_ep]
  i_start <- ep_idx[max(which(mel[ep_idx] == max(mel[ep_idx])))]

  # 3. end of the limb: local minimum before the next sustained rise
  k_rise <- which(r$values & starts > i_cross & (r$lengths >= min_run | ends == n))
  if (length(k_rise) == 0) {
    i_end <- n
  } else {
    m <- starts[k_rise[1]] - 1
    while (m - 1 >= i_cross && mel[m - 1] < mel[m]) m <- m - 1
    i_end <- m
  }

  limb <- data[i_start:i_end, , drop = FALSE]
  list(limb = limb,
       info = list(ok = TRUE, reason = NA_character_,
                   start = data$datetime[i_start], end = data$datetime[i_end],
                   crossing = data$datetime[i_cross],
                   n = nrow(limb), n_above = sum(limb$melatonin > threshold),
                   peak = mel[i_start], trough = mel[i_end],
                   next_rise_found = length(k_rise) > 0))
}
