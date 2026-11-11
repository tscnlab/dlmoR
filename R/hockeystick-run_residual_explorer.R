#' Launch the residual explorer Shiny app
#'
#' Launches the interactive residual explorer for inspecting
#' residual landscapes generated during DLMO estimation.
#'
#' The app provides interactive visualization of residual
#' surfaces, low-residual clusters, radial profiles,
#' steepness metrics, and reverse gradient ascent paths
#' used to evaluate DLMO solution stability and uniqueness.
#'
#' @param dlmo_result A result object returned by
#'   [calculate_dlmo()]. The object can have any name
#'   in the user's environment and will be made available
#'   internally to the Shiny app as `dlmo_result`.
#'
#' @details
#' The residual explorer is intended for interactive
#' exploration of the residual grid and cluster structure
#' underlying a DLMO estimate.
#'
#' Example workflow:
#'
#' \preformatted{
#' res <- calculate_dlmo(...)
#' run_residual_explorer(res)
#' }
#'
#' @return Launches an interactive Shiny application for exploring residual landscape resulting from hockey-stick model fits.
#'
#' @export
run_residual_explorer <- function(dlmo_result) {

  # Validate input
  if (missing(dlmo_result)) {
    stop(
      "Please provide a result from calculate_dlmo().\n\n",
      "Example:\n",
      "res <- calculate_dlmo(...)\n",
      "run_residual_explorer(res)",
      call. = FALSE
    )
  }

  # Ensure shiny is installed
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop(
      "Package 'shiny' is required to run the residual explorer.\n",
      "Install it with:\n",
      "install.packages('shiny')",
      call. = FALSE
    )
  }

  # Locate installed app
  app_dir <- system.file(
    "shiny",
    "residual_explorer",
    package = "dlmoR"
  )

  if (!nzchar(app_dir)) {
    stop(
      "Could not find the residual_explorer Shiny app.",
      call. = FALSE
    )
  }

  # Put object into global env temporarily for app access
  assign("dlmo_result", dlmo_result, envir = .GlobalEnv)

  # Launch app
  shiny::runApp(app_dir)
}
