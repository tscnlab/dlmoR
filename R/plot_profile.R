#' Add DLMO Fit Lines to a Melatonin Profile Plot
#'
#' This function overlays the estimated DLMO fit lines (base and ascending phase) onto a melatonin profile plot.
#' The fit can be either **linear** or **parabolic**, depending on the `dlmoFit` results.
#'
#' @param plot A ggplot object. The base plot of the melatonin profile.
#' @param profile_data A dataframe containing melatonin concentration data, with required columns:
#'   \itemize{
#'     \item `datetime` (POSIXct) – Time of measurement.
#'     \item `melatonin` (numeric) – Melatonin concentration levels.
#'     \item `base` (binary) – Indicator column for base segment points.
#'     \item `ascending` (binary) – Indicator column for ascending segment points.
#'   }
#' @param dlmoFit A list containing the estimated DLMO inflection point and fit parameters, including:
#'   \itemize{
#'     \item `inflection_point$x` (decimal hours) – The estimated x-coordinate of the inflection point.
#'     \item `inflection_point$y` (numeric) – The estimated y-coordinate of the inflection point.
#'     \item `base_params` (numeric) – Slope of the linear base fit.
#'     \item `ascending_params` (list) – Either a single slope (for linear) or coefficients `a, b, c` (for parabolic).
#'   }
#'
#' @return A ggplot object with the fitted lines added.
#' @details
#' - The **base fit** is always **linear** and is drawn from the start of the base segment to the inflection point.
#' - The **ascending fit** can be either **linear** or **parabolic**:
#'   \itemize{
#'     \item **Linear fit:** A simple straight line extending from the inflection point to the last ascending point.
#'     \item **Parabolic fit:** A nonlinear curve extending from the inflection point to the last ascending point.
#'   }
#' - Uses `posixct_to_decimal()` and `decimal_to_posixct()` to convert between time formats.
#'
#' @examples
#' \dontrun{
#'   plot <- ggplot2::ggplot(profile_data, ggplot2::aes(x = datetime, y = melatonin)) +
#'       ggplot2::geom_point()
#'   dlmoFit <- list(inflection_point = list(x = 20.5, y = 3.2),
#'                   base_params = 0.1,
#'                   ascending_params = list(a = -0.03, b = 2.1, c = 0.5))
#'   plot_fit(plot, profile_data, dlmoFit)
#' }
#' @export
plot_fit <- function(plot, profile_data, dlmoFit) {
  # Convert base and ascending timepoints to decimal hours
  xstart_num <- posixct_to_decimal(dplyr::filter(profile_data, base == 1)$datetime[1], profile_data$datetime)
  xend_num <- posixct_to_decimal(tail(dplyr::filter(profile_data, ascending == 1)$datetime, n = 1), profile_data$datetime)

  # Convert the inflection point x back to POSIXct for plotting
  ipx_posix <- decimal_to_posixct(dlmoFit$inflection_point$x, profile_data$datetime)

  # --- BASE SEGMENT FIT ---
  # Always linear: y = m * (x - poi_x) + poi_y
  plot <- plot +
    ggplot2::geom_segment(
      x = dplyr::filter(profile_data, base == 1)$datetime[1],  # Start from first base point
      y = dlmoFit$base_params * (xstart_num - dlmoFit$inflection_point$x) + dlmoFit$inflection_point$y,  # Compute y-start using linear equation
      xend = ipx_posix,  # End at the inflection point
      yend = dlmoFit$inflection_point$y,
      color = "darkgray",
      size = 1
    )

  # --- ASCENDING SEGMENT FIT ---
  if (length(dlmoFit$ascending_params) == 1) {  # If ascending fit is linear
    plot <- plot +
      ggplot2::geom_segment(
        x = ipx_posix,  # Start from inflection point
        y = dlmoFit$inflection_point$y,
        xend = dplyr::filter(profile_data, ascending == 1)$datetime[length(dplyr::filter(profile_data, ascending == 1)$datetime)],  # Last ascending point
        yend = dplyr::filter(profile_data, ascending == 1)$melatonin[length(dplyr::filter(profile_data, ascending == 1)$melatonin)],  # Last ascending melatonin value
        color = "darkgray",
        size = 1
      )
  } else {  # If ascending fit is parabolic
    plot <- plot +
      ggplot2::geom_function(
        fun = function(x) {
          # Convert x (POSIXct) to decimal hours
          x_numeric <- posixct_to_decimal(x, profile_data$datetime)

          # Evaluate the parabolic function y = ax^2 + bx + c
          y <- dlmoFit$ascending_params$a * x_numeric^2 +
            dlmoFit$ascending_params$b * x_numeric +
            dlmoFit$ascending_params$c
          return(y)
        },
        color = 'darkgray',
        size = 1,
        n = 1000,
        xlim = c(
          decimal_to_posixct(dlmoFit$inflection_point$x, profile_data$datetime),  # Inflection point
          decimal_to_posixct(xend_num, profile_data$datetime)  # Last ascending point
        )
      )
  }

  # Format x-axis with time labels
  ggplot2::scale_x_datetime(
    limits = c(decimal_to_posixct(dlmoFit$inflection_point$x, profile_data$datetime),
               decimal_to_posixct(xend_num, profile_data$datetime)),
    date_labels = "%H:%M",  # Display hour and minute
    date_breaks = "1 hour"   # Set breaks at every hour
  )

  return(plot)
}



#' Add DLMO Inflection Point to the Plot
#'
#' This function adds a **DLMO inflection point** as a distinct marker to an existing melatonin profile plot.
#' The inflection point is represented by a **filled deep pink square**.
#'
#' @param plot A ggplot object. The base plot of the melatonin profile.
#' @param profile_data A dataframe containing melatonin concentration data, with required columns:
#'   \itemize{
#'     \item `datetime` (POSIXct) – Time of measurement.
#'   }
#' @param dlmoip A list containing the estimated **inflection point** coordinates:
#'   \itemize{
#'     \item `x` (numeric) – The estimated inflection point time (in decimal hours).
#'     \item `y` (numeric) – The estimated melatonin concentration at the inflection point.
#'   }
#'
#' @return A ggplot object with the inflection point added.
#'
#' @details
#' - The **inflection point** represents the DLMO time stamp, at which melatonin concentration begins to rise significantly.
#' - The function **converts** the inflection point `x` from **decimal hours** back to **POSIXct** for proper placement.
#' - The point is **styled** as a **deep pink** filled square.
#'
#' @examples
#' \dontrun{
#'   plot <- ggplot2::ggplot(profile_data, ggplot2::aes(x = datetime, y = melatonin)) +
#'       ggplot2::geom_point()
#'   dlmoip <- list(x = 20.5, y = 3.2)  # Inflection point at 20.5 decimal hours, y = 3.2 pg/mL
#'   plot_ip(plot, profile_data, dlmoip)
#' }
#' @export
plot_ip <- function(plot, profile_data, dlmoip) {
  plot <- plot +
    ggplot2::geom_point(
      x = decimal_to_posixct(dlmoip$x, profile_data$datetime),  # Convert decimal hours to POSIXct for x-axis
      y = dlmoip$y,  # Inflection point melatonin concentration
      color = "deeppink4",  # Outline color
      fill = "deeppink4",  # Fill color
      size = 3,  # Point size
      shape = 23  # Square shape with fill
    )

  return(plot)
}



#' Add Region of Interest (ROI) to a Plot
#'
#' This function overlays the **Region of Interest (ROI)** on an existing melatonin profile plot.
#' It adds a **horizontal segment** to mark the ROI range and optionally a **shaded rectangle** to highlight the full area that is parsed during the DLMO search.
#'
#' @param plot A ggplot object. The base plot to which the ROI should be added.
#' @param roi A list containing the **Region of Interest (ROI) coordinates**:
#'   \itemize{
#'     \item `x_start` (POSIXct) – The starting time of the ROI.
#'     \item `x_end` (POSIXct) – The ending time of the ROI.
#'     \item `y_min` (numeric) – The minimum melatonin concentration in the ROI.
#'     \item `y_max` (numeric) – The maximum melatonin concentration in the ROI.
#'   }
#' @param roi_line_only Logical. If `TRUE`, only the **horizontal segment** is plotted.
#'   If `FALSE`, a **semi-transparent rectangle** is added to indicate the **full ROI region**.
#'
#' @return A ggplot object with the ROI overlay.
#'
#' @details
#' - The **horizontal segment** is plotted at `y = -0.2` to clearly mark the time range of the ROI.
#' - If `roi_line_only = FALSE`, the **rectangle spans** from `y_min` to `y_max` with a light **orchid3** fill.
#'
#' @examples
#' \dontrun{
#'   plot <- ggplot2::ggplot(profile_data, ggplot2::aes(x = datetime, y = melatonin)) +
#'       ggplot2::geom_point()
#'   roi <- list(x_start = as.POSIXct("2024-04-16 18:00:00"),
#'               x_end = as.POSIXct("2024-04-16 22:00:00"),
#'               y_min = 0.2, y_max = 2.5)
#'   plot_roi(plot, roi, roi_line_only = FALSE)
#' }
#' @export
plot_roi <- function(plot, roi, roi_line_only) {
  # Add the ROI horizontal segment to indicate the time range
  plot <- plot +
    ggplot2::geom_segment(
      x = roi$x_start, xend = roi$x_end,  # Horizontal ROI range
      y = -0.2, yend = -0.2,  # Fixed y-position for clarity
      color = "orchid3", size = 1.5  # Color and thickness of the line
    )

  # Optionally add a shaded rectangle for the full ROI region
  if (!roi_line_only) {
    plot <- plot +
      ggplot2::geom_rect(
        xmin = roi$x_start, xmax = roi$x_end,  # Time range of ROI
        ymin = roi$y_min, ymax = roi$y_max,  # Melatonin concentration range
        fill = "orchid3", alpha = 0.01  # Light transparent shading
      )
  }

  return(plot)
}



#' Add ROI Heatmap to a Plot
#'
#' This function overlays a **heatmap of residuals** on a melatonin profile plot, visualizing
#' the **goodness of fit** across the **Region of Interest (ROI)** search grid.
#'
#' @param plot A ggplot object. The base plot to which the heatmap should be added.
#' @param data A dataframe containing the melatonin profile with `datetime` and `melatonin` values.
#' @param roi_grid A dataframe containing **roi grid points** (ROI search area).
#' @param residuals A numeric vector of residuals corresponding to `roi_grid`.
#'
#' @return A ggplot object with the residual heatmap overlay.
#'
#' @details
#' - The **color gradient** represents residual values, with lower residuals indicating **better fit**.
#' - Residuals are **log-transformed** for visualization.
#'
#' @examples
#' \dontrun{
#'   plot <- ggplot2::ggplot(profile_data, ggplot2::aes(x = datetime, y = melatonin)) +
#'       ggplot2::geom_point()
#'   plot <- plot_roi_heatmap(plot, data = profile_data, roi_grid = grid, residuals = res)
#' }
#' @export
plot_roi_heatmap <- function(plot, data = NULL, roi_grid = NULL, residuals = NULL) {
  # Ensure at least one grid is selected for visualization

    # Convert x-coordinates from decimal hours to POSIXct timestamps
    dt_roi_grid <- data.frame(
      x = decimal_to_posixct(roi_grid$x, data$datetime),
      y = roi_grid$y
    )

    # Add coarse grid heatmap to plot
    plot <- plot + ggplot2::geom_point(
      data = dt_roi_grid, ggplot2::aes(x = x, y = y, color = log(residuals)),
      size = 1.25
    ) +
      ggplot2::scale_color_gradient(low = "cyan", high = "deeppink")  # Color gradient from blue to pink


  return(plot)
}



#' Overlay a Parallelogram on a Melatonin Profile Plot
#'
#' This function adds a **parallelogram visualization** to a melatonin profile plot.
#' The best-fit parallelogram encompasses all ascending points with the smallest area parallelogram
#' possible. The ratio of parallelogram latera edge to long diagonal is used as a criterion
#' for trimming the ascending region to ensure that only a sufficiently fast enough melatonin
#' rise is taken into consideration when fitting the profile and searching for the DLMO point.
#'
#' @param plot A ggplot object. The base plot to which the parallelogram will be added.
#' @param profile_data A dataframe containing melatonin concentration data with `datetime` values.
#' @param pll_result A list containing the **optimized parallelogram parameters**, including:
#'   \itemize{
#'     \item `pll_datetime_0`: POSIXct timestamp for the left boundary.
#'     \item `pll_datetime_1`: POSIXct timestamp for the right boundary.
#'     \item `pll_slope`: Numeric. The slope of the parallelogram's edges.
#'     \item `corners`: A list of four corner coordinates (`ll`, `lr`, `ur`, `ul`).
#'   }
#' @return A ggplot object with the parallelogram overlay.
#'
#' @details
#' - The **parallelogram** is generated using the **optimized fit** from `parallelogram_fit()`.
#' - The **edges** of the parallelogram indicate a **bounded region of interest (ROI)**.
#' - Two **diagonal dashed lines** highlight the **shape constraints** used during optimization.
#' - Uses **red fill** with **transparency (alpha = 0.1)** to avoid obscuring data.
#'
#' @examples
#' \dontrun{
#'   plot <- ggplot2::ggplot(profile_data, ggplot2::aes(x = datetime, y = melatonin)) +
#'       ggplot2::geom_point()
#'   plot <- plot_parallelogram(plot, profile_data, pll_result)
#' }
#' @export
plot_parallelogram <- function(plot, profile_data, pll_result) {
  # Ensure parallelogram results are provided
  if (is.null(pll_result)) {
    stop("pll_result must be provided to plot the parallelogram.")
  }

  # Extract optimized parallelogram parameters
  x0_posix <- pll_result$pll_datetime_0  # Left boundary
  x1_posix <- pll_result$pll_datetime_1  # Right boundary
  slope <- pll_result$pll_slope          # Edge slope
  corners <- pll_result$corners          # Corner coordinates

  # Convert corner x-values from numeric (decimal hours) to POSIXct timestamps
  corners_datetime <- lapply(corners, function(corner) {
    list(datetime = decimal_to_posixct(corner[1], profile_data$datetime),
         melatonin = corner[2])  # Preserve melatonin concentration
  })

  # Create a dataframe for the parallelogram
  parallelogram_df <- do.call(rbind, lapply(corners_datetime, as.data.frame))

  # Identify diagonal points (1 → 3 and 2 → 4) for visualization
  diagonal_1 <- parallelogram_df[c(1, 3), ]  # Connect lower-left to upper-right
  diagonal_2 <- parallelogram_df[c(2, 4), ]  # Connect lower-right to upper-left

  # Add parallelogram overlay and diagonals to the plot
  plot <- plot +
    ggplot2::geom_polygon(
      data = parallelogram_df,
      ggplot2::aes(x = .data$datetime, y = .data$melatonin),
      fill = "red", alpha = 0.1  # Semi-transparent red fill
    ) +
    ggplot2::geom_line(
      data = diagonal_1,
      ggplot2::aes(x = .data$datetime, y = .data$melatonin),
      color = "red", linetype = "dashed", size = 0.5  # Dashed diagonal line 1
    ) +
    ggplot2::geom_line(
      data = diagonal_2,
      ggplot2::aes(x = .data$datetime, y = .data$melatonin),
      color = "red", linetype = "dashed", size = 0.5  # Dashed diagonal line 2
    )

  return(plot)
}


#' Generate a Melatonin Profile Plot
#'
#' This function creates a plot of melatonin concentration over time with various
#' optional overlays, including **DLMO fit lines, DLMO time stamp, threshold lines,**
#' **melatonin profile segment highlights, best-fit parallelogram, and ROI residual heatmaps**.
#'
#' @param profile_data A dataframe containing melatonin concentration data with `datetime` values.
#' @param show_threshold Logical. If `TRUE`, adds a **horizontal threshold line** at `threshold` value.
#' @param threshold Numeric. The threshold value for melatonin concentration (default `2.3 pg/mL`).
#' @param show_segments Logical. If `TRUE`, highlights **base, intermediate, and ascending segments**.
#' @param show_parallelogram Logical. If `TRUE`, overlays a **parallelogram** on the plot.
#' @param pll_result List containing parallelogram parameters, including:
#'   \itemize{
#'     \item `pll_datetime_0`: POSIXct timestamp for the left boundary.
#'     \item `pll_datetime_1`: POSIXct timestamp for the right boundary.
#'     \item `pll_slope`: Numeric. The slope of the parallelogram edges.
#'     \item `corners`: List of corner coordinates (`ll`, `lr`, `ur`, `ul`).
#'   }
#' @param show_roi Logical. If `TRUE`, overlays the **region of interest (ROI)**.
#' @param roi_line_only Logical. If `TRUE`, only plots the **ROI boundaries**, not the shaded area.
#' @param roi List containing ROI boundaries (`x_start`, `x_end`, `y_min`, `y_max`).
#' @param show_dlmoIP Logical. If `TRUE`, marks the **DLMO inflection point**.
#' @param dlmoFit List containing **DLMO fit results**, including:
#'   \itemize{
#'     \item `inflection_point`: A list with `x` (decimal hours) and `y` (melatonin level).
#'     \item `base_params`: Parameters of the **base segment fit**.
#'     \item `ascending_params`: Parameters of the **ascending segment fit**.
#'     \item `grid`: search grid for inflection point.
#'     \item `res`: Residuals from grid search.
#'   }
#' @param show_fit Logical. If `TRUE`, overlays **DLMO fit lines**.
#' @param show_roi_heatmap Logical. If `TRUE`, adds a **heatmap** for **ROI residuals**.
#' @return A `ggplot2` object with the melatonin profile and optional overlays.
#'
#' @details
#' - This function **plots melatonin concentration** over time with flexible overlays.
#' - The **DLMO fit** is visualized as a **piecewise-linear or parabolic** fit.
#' - The **inflection point** is highlighted using a **pink marker**.
#' - A **parallelogram** can be overlaid to show **truncated ascending segments**.
#' - The **ROI heatmap** provides a **residuals-based visualization** of the inflection search.
#'
#' @examples
#' \dontrun{
#'   plot <- plot_profile(profile_data, show_threshold = TRUE, show_segments = TRUE, show_dlmoIP = TRUE)
#'   print(plot)
#' }
#' @export
plot_profile <- function(profile_data, show_threshold = TRUE, threshold = 2.3,
                         show_segments = TRUE, show_parallelogram = FALSE, pll_result = NULL,
                         show_roi = FALSE, roi_line_only = TRUE, roi = NULL,
                         show_dlmoIP = TRUE, dlmo = NULL, dlmoFit = NULL, show_fit = FALSE,
                         show_roi_heatmap = FALSE, plot_coarse = FALSE) {

  # Function to safely replace NULL values
  ensure_non_null <- function(value, fallback = "NA") {
    if (is.null(value)) fallback else value
  }

  # Define title and coefficients conditionally
  dlmo_values <- if (!is.null(dlmoFit)) {
    if (is.null(dlmo$fine) || plot_coarse){
      # If no fine fit exists, use COARSE fit
      dlmo_time <- dlmo$coarse$time
      plot_title <- paste("Coarse Fit - DLMO time:", as.character(dlmo_time))
      a1 <- dlmo$coarse$fit_lines$base$m
      b1 <- dlmo$coarse$fit_lines$base$b
      a2 <- dlmo$coarse$fit_lines$ascending$m
      b2 <- dlmo$coarse$fit_lines$ascending$b
      c2 <- NULL  # Coarse fit is always linear
    } else {       # If fine fit exists
      if(!plot_coarse){ #if plot type is fine
        dlmo_time <- dlmo$fine$time
        plot_title <- paste("Fine Fit - DLMO time:", as.character(dlmo_time))
        a1 <- dlmo$fine$fit_lines$base$m
        b1 <- dlmo$fine$fit_lines$base$m
        if (length(dlmo$fine$fit_lines$ascending$params) == 1) {
          # Linear fine fit
          a2 <- dlmo$fine$fit_lines$ascending$m
          b2 <- dlmo$fine$fit_lines$ascending$b
          c2 <- NULL  # No quadratic term
        } else {
          # Parabolic fine fit
          a2 <- dlmo$fine$fit_lines$ascending$a
          b2 <- dlmo$fine$fit_lines$ascending$b
          c2 <- dlmo$fine$fit_lines$ascending$c
        }
      }
    }

    list(plot_title = plot_title, a1 = a1, b1 = b1, a2 = a2, b2 = b2, c2 = c2)
  } else {
    list(plot_title = "Melatonin profile", a1 = NULL, b1 = NULL, a2 = NULL, b2 = NULL, c2 = NULL)
  }

  # Unpack the returned values
  plot_title <- dlmo_values$plot_title
  a1 <- dlmo_values$a1
  b1 <- dlmo_values$b1
  a2 <- dlmo_values$a2
  b2 <- dlmo_values$b2
  c2 <- dlmo_values$c2



  # Define subtitle conditionally
  if (!is.null(a2) && !is.null(b2)) {
    if (is.null(c2)) {  # If c2 is NULL, it's a linear fit
      subtitle_text <- bquote(
        bold(f[e]) == .(ensure_non_null(a1)) * x[t] + .(ensure_non_null(b1)) * "," ~ "\n" ~
          bold(f[l]) == .(ensure_non_null(a2)) * x[t] + .(ensure_non_null(b2))
      )
    } else {  # If c2 exists, it's a parabolic fit
      subtitle_text <- bquote(
        bold(f[e]) == .(ensure_non_null(a1)) * x[t] + .(ensure_non_null(b1)) * "," ~ "\n" ~
          bold(f[l]) == .(ensure_non_null(a2)) * x[t]^2 + .(ensure_non_null(b2)) * x[t] + .(ensure_non_null(c2))
      )
    }
  } else {
    subtitle_text <- "No valid fit available"
  }


  # Initialize base plot
  plot <- ggplot2::ggplot(profile_data, ggplot2::aes(x = .data$datetime, y = .data$melatonin)) +
    ggplot2::geom_point(color = 'grey', size = 2) +
    ggplot2::geom_line(color = 'grey', linetype = "dotted", size = 1) +
    ggplot2::scale_x_datetime(labels = scales::date_format("%H:%M"), date_breaks = "2 hours") +

    ggplot2::labs(
      title = plot_title,
      subtitle = subtitle_text,
      x = "Local time [hh:mm]",
      y = "Melatonin concentration [pg/mL]"
    ) +

    ggplot2::theme_minimal() +

    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 12, hjust = 0)
    ) +

    # Corrected legend position
    ggplot2::theme(legend.position = "none") +

    ggplot2::scale_linetype_manual(values = c("Full Profile" = "dotted"))


  # Overlay Parallelogram if enabled
  if (show_parallelogram) {
    plot <- plot_parallelogram(plot, profile_data, pll_result)
  }

  # Overlay Region of Interest (ROI)
  if (show_roi) {
    plot <- plot_roi(plot, roi, roi_line_only)
  }

  # Add ROI heatmap if enabled
  if (show_roi_heatmap) {
    plot <- plot_roi_heatmap(plot, data = profile_data, roi_grid = dlmoFit$grid,residuals = dlmoFit$res)
  }

  # Re-add full profile line to ensure clarity
  plot <- plot + ggplot2::geom_line(color = 'grey', linetype = "dotted", size = 1)

  # Overlay DLMO Fit if enabled
  if (show_fit) {
    plot <- plot_fit(plot, profile_data, dlmoFit)
  }

  # Mark Inflection Point (DLMO) if enabled
  if (show_dlmoIP) {
    plot <- plot_ip(plot, profile_data, dlmoFit$inflection_point)
  }

  # Add threshold line if enabled
  if (show_threshold) {
    plot <- plot + ggplot2::geom_hline(ggplot2::aes(yintercept = threshold), color = "burlywood3", size = 1)
  }

  # Highlight Base and Ascending Segments
  if (show_segments) {
    plot <- plot +
      ggplot2::geom_point(
        data = dplyr::filter(profile_data, .data$base == 1),
        ggplot2::aes(x = .data$datetime, y = .data$melatonin, color = "Base Segment"),
        color = '#56B4E9',
        size = 2
      ) +
      ggplot2::geom_point(
        data = dplyr::filter(profile_data, .data$ascending == 1),
        ggplot2::aes(x = .data$datetime, y = .data$melatonin, color = "Ascending Segment"),
        color = 'lightgreen',
        size = 2
      )

    # Highlight Intermediate Segment if present
    if ("intermediate" %in% colnames(profile_data)) {
      plot <- plot +
        ggplot2::geom_point(
          data = dplyr::filter(profile_data, .data$intermediate == 1),
          ggplot2::aes(x = .data$datetime, y = .data$melatonin, color = "Intermediate Segment"),
          color = 'darkgoldenrod1',
          size = 2
        )
    }
  }

  # Return the final plot
  return(plot)
}
