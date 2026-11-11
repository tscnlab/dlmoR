#' Calculate Dim-Light Melatonin Onset (DLMO)
#'
#' Calculates the DLMO from a time series of melatonin concentrations. This is
#' the single entry point for all DLMO methods in `dlmoR`: choose the method
#' with `method`, and pass the settings of that method as further named
#' arguments. The data can be provided directly as a data frame or as a .csv
#' file.
#'
#' @param data A data frame with columns `datetime` (POSIXct) and `melatonin`
#'   (numeric). If `NULL`, use the `file_path` parameter to load data.
#' @param file_path A string specifying the path to a CSV file containing the
#'   data. The file must have two columns: `datetime` (POSIXct) and `melatonin`
#'   (numeric).
#' @param method The DLMO method: `"hockeystick"` (default, the hockey-stick
#'   method) or `"3k"` (the 3k threshold method). See the sections below.
#' @param dlmo_off Logical. If `FALSE` (default), the dim-light melatonin
#'   onset (DLMO) is calculated. If `TRUE`, the dim-light melatonin offset
#'   (DLMOoff) is calculated instead. DLMOoff is not implemented yet for any
#'   method, so `dlmo_off = TRUE` currently stops with an error.
#' @param ... Named settings of the chosen method, for example `threshold` and
#'   `fine_flag` for `method = "hockeystick"`, or `n_baseline` and `k_sd` for
#'   `method = "3k"`. The arguments of each method are listed in the sections
#'   below. Arguments that the chosen method does not use give an error.
#'
#' @return The result of the chosen method; see the "Value" part of the
#'   section of each method below.
#'
#' @seealso [plot_raw_profile()] to look at a profile before choosing a method,
#'   and [run_residual_explorer()] for hockey-stick results.
#'
#' @order 1
#' @export
calculate_dlmo <- function(data = NULL, file_path = NULL,
                           method = "hockeystick", dlmo_off = FALSE, ...) {

  if (is.null(data) && is.null(file_path)) {
    stop("You must provide either `data` or `file_path`.", call. = FALSE)
  }
  registry <- .dlmo_method_registry()
  if (!is.character(method) || length(method) != 1 || is.na(method)) {
    stop("`method` must be a single string, one of: ", .quote_names(names(registry)), ".",
         call. = FALSE)
  }
  if (!method %in% names(registry)) {
    stop("Unknown `method` = ", .quote_names(method), ". Use one of: ",
         .quote_names(names(registry)), ".", call. = FALSE)
  }
  if (!is.logical(dlmo_off) || length(dlmo_off) != 1 || is.na(dlmo_off)) {
    stop("`dlmo_off` must be TRUE or FALSE.", call. = FALSE)
  }

  marker <- if (dlmo_off) "offset" else "onset"
  method_fun <- registry[[method]][[marker]]
  if (is.null(method_fun)) {
    stop("DLMOoff (`dlmo_off = TRUE`) is not implemented yet for method = \"",
         method, "\".", call. = FALSE)
  }
  arg_names <- names(list(...))
  if (is.null(arg_names)) arg_names <- rep("", ...length())
  .check_method_args(arg_names, method_fun, method)

  method_fun(data = data, file_path = file_path, ...)
}


#' DLMO Method Registry
#'
#' Lists the calculation function of every method, for the onset (DLMO) and
#' the offset (DLMOoff). To add a method, write its calculation function in
#' files named `<method>-*.R` and add an entry here; to add an offset
#' calculation, replace the `NULL` with the offset function (kept in
#' `<method>-offset-*.R`). Every calculation function takes `data` and
#' `file_path` as its first two arguments.
#'
#' The names of the entries are the accepted values of the `method` argument
#' of [calculate_dlmo()]; its default (`"hockeystick"`) must be one of them.
#'
#' @return A named list with one element per method, each a list with the
#'   functions `onset` and `offset` (`NULL` if not implemented).
#' @keywords internal
#' @noRd
.dlmo_method_registry <- function() {
  list(
    hockeystick = list(onset = .calculate_dlmo_hockeystick, offset = NULL),
    `3k`        = list(onset = .calculate_dlmo_3k,          offset = NULL)
  )
}


#' Check the Method-Specific Arguments Given to calculate_dlmo()
#'
#' @param arg_names Names of the arguments in `...` (`""` for unnamed ones).
#' @param method_fun The calculation function of the method.
#' @param method The method name, for the error message.
#' @return `TRUE`, invisibly; stops if an argument is unnamed or not used by
#'   the method.
#' @keywords internal
#' @noRd
.check_method_args <- function(arg_names, method_fun, method) {
  if (length(arg_names) == 0) return(invisible(TRUE))
  if (any(arg_names == "")) {
    stop("Method settings must be named, for example `threshold = 5`.", call. = FALSE)
  }
  valid <- setdiff(names(formals(method_fun)), c("data", "file_path"))
  # partial matching as in a normal R function call (e.g. `fine` for `fine_flag`)
  unknown <- arg_names[is.na(pmatch(arg_names, valid, duplicates.ok = TRUE))]
  if (length(unknown) > 0) {
    stop("Argument(s) not used by method = \"", method, "\": ",
         paste0("`", unknown, "`", collapse = ", "), ". Arguments of this method: ",
         paste0("`", valid, "`", collapse = ", "), ".", call. = FALSE)
  }
  invisible(TRUE)
}

.quote_names <- function(x) paste0("\"", x, "\"", collapse = ", ")
