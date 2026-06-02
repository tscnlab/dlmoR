# =============================================================================
# Residual heatmap explorer (Shiny) for dlmoR DLMO fits
#
# This Shiny app is included in the publication repository for transparency.
# It will be migrated into the dlmoR package in a future release.
#
# Purpose
#   Interactive exploration of the residual surface produced by the dlmoR
#   hockey-stick fitting procedure. The app:
#     1) visualizes the residual grid as a heatmap,
#     2) thresholds low-residual points and labels 4-connected clusters,
#     3) reports cluster summary statistics (nominal time and relative to DLMO),
#     4) visualizes time-offset distributions by cluster.
#
# Required input object
#   This script assumes an object named `dlmo_result` exists in the R session
#   before launching the app. `dlmo_result` must be the output of:
#     dlmoR::calculate_dlmo(...)
#
# Output behavior
#   This app writes publication-ready SVG figures to the `results/` directory
#   as a side effect when plots are rendered. The following files are created
#   or overwritten during use:
#     - results/residual_heatmap.svg
#     - results/offset_violin.svg
#
#   This behavior reflects the workflow used to generate figures for the
#   associated publication.
#
# =============================================================================

library(shiny)
library(viridis)
library(ggplot2)
library(dplyr)
library(tidyr)

# -----------------------------------------------------------------------------
# Helper: 4-connected cluster labeling on a logical matrix
# -----------------------------------------------------------------------------
label_clusters <- function(binary_mask) {
  nr <- nrow(binary_mask)
  nc <- ncol(binary_mask)

  visited <- matrix(FALSE, nr, nc)
  labels <- matrix(0, nr, nc)
  label_id <- 0

  flood_fill <- function(i, j, label_id) {
    stack <- list(c(i, j))
    while (length(stack) > 0) {
      loc <- stack[[1]]
      stack <- stack[-1]

      x <- loc[1]
      y <- loc[2]

      if (x < 1 || x > nr || y < 1 || y > nc) next
      if (visited[x, y]) next
      if (!binary_mask[x, y]) next

      visited[x, y] <<- TRUE
      labels[x, y] <<- label_id

      stack <- c(
        stack,
        list(c(x - 1, y)),
        list(c(x + 1, y)),
        list(c(x, y - 1)),
        list(c(x, y + 1))
      )
    }
  }

  for (i in 1:nr) {
    for (j in 1:nc) {
      if (binary_mask[i, j] && !visited[i, j]) {
        label_id <- label_id + 1
        flood_fill(i, j, label_id)
      }
    }
  }

  list(label_matrix = labels, n_clusters = label_id)
}

# -----------------------------------------------------------------------------
# Shiny app definition
# -----------------------------------------------------------------------------
shinyApp(
  ui = fluidPage(
    titlePanel(HTML("<b><i>dlmoR</i></b>: Residual heatmap explorer")),

    sidebarLayout(
      sidebarPanel(
        sliderInput(
          "top_percent",
          "Clustering threshold (% lowest residuals):",
          min = 0.1, max = 100, value = 5, step = 0.5
        ),
        numericInput(
          "top_percent_num",
          "Enter threshold %:",
          value = 5, min = 0.1, max = 100, step = 0.5
        ),
        numericInput("min_cluster_size", "Minimum cluster size:", value = 10),
        uiOutput("residual_range"),

        tags$div(
          style = "margin-top: 25px; text-align: left;",
          tags$b("If you use this app or package, please cite:"),
          tags$hr(style = "border-top: 1px solid black; margin-top: 10px; margin-bottom: 10px;"),

          HTML("Thalji, Salma M., and Manuel Spitschan (2026). <i>dlmoR: An open-source R package for the dim-light melatonin onset (DLMO) hockey-stick method</i>. R package version 2.0.0. <a href='https://github.com/tscnlab/dlmoR'>https://github.com/tscnlab/dlmoR</a>"),
          tags$br(), tags$br(),

          HTML(
            "Thalji, S. M., &amp; Spitschan, M. (2026).
            <i><a href='https://doi.org/10.1177/07487304251389994' target='_blank'>
            dlmoR: An Open-Source R Package for the Dim-Light Melatonin Onset (DLMO) Hockey-Stick Method.
            </a></i>
            <i>Journal of Biological Rhythms</i>, 41(3), 301–323."
                    ),
          tags$br(), tags$br(),

          tags$details(
            tags$summary(
              tags$span(
                "Show BibTeX",
                style = "cursor: pointer; background-color: #f0f0f0; border: 1px solid #ccc; border-radius: 4px; padding: 4px 8px; font-weight: bold;"
              )
            ),
            tags$pre("
@Manual{thalji2026dlmor,
  title = {dlmoR: An open-source R package for the dim-light melatonin onset (DLMO) hockey-stick method},
  author = {Salma M. Thalji and Manuel Spitschan},
  year = {2026},
  note = {R package version 2.0.0},
  url = {https://github.com/tscnlab/dlmoR}
}

@Article{thalji2026dlmor-paper,
  title = {dlmoR: An open-source R package for the dim-light melatonin onset (DLMO) hockey-stick method},
  author = {Salma M. Thalji and Manuel Spitschan},
  journal = {Journal of Biological Rhythms},
  year = {2026},
  doi = {10.1177/07487304251389994}
}
")
          )
        )
      ),

      mainPanel(
        plotOutput("residPlot"),
        div(style = "margin-top:30px;"),
        plotOutput("offsetViolinPlot"),
        div(style = "margin-top:30px;"),
        h4(tags$b("DLMO time:"), dlmo_result$dlmo$fine$time),
        div(style = "margin-top:30px;"),
        h4(tags$b("Cluster stats: nominal")),
        tableOutput("clusterNominalStats"),
        div(style = "margin-top:30px;"),
        h4(tags$b("Cluster stats: time difference relative to DLMO")),
        tableOutput("clusterRelativeStats")
      )
    )
  ),

  server = function(input, output, session) {

    # -------------------------------------------------------------------------
    # Convert dlmo_result grid output into a matrix (residual surface)
    # -------------------------------------------------------------------------
    get_matrix <- reactive({
      req(dlmo_result)

      grid_df <- dlmo_result$ip$grid_small
      x_vals <- sort(unique(grid_df$x))
      y_vals <- sort(unique(grid_df$y))

      mat <- matrix(NA, nrow = length(y_vals), ncol = length(x_vals))

      for (i in seq_along(dlmo_result$ip$res_small)) {
        x <- grid_df$x[i]
        y <- grid_df$y[i]
        row <- which(y_vals == y)
        col <- which(x_vals == x)
        mat[row, col] <- dlmo_result$ip$res_small[i]
      }

      list(mat = mat, x = x_vals, y = y_vals)
    })

    # -------------------------------------------------------------------------
    # Cluster labeling on the residual surface
    # -------------------------------------------------------------------------
    full_clusters <- reactive({
      req(input$top_percent, input$min_cluster_size)

      m <- get_matrix()
      mat <- m$mat

      threshold <- quantile(mat, input$top_percent / 100, na.rm = TRUE)
      mask <- mat <= threshold

      labeled <- label_clusters(mask)
      labels <- labeled$label_matrix

      sizes <- table(labels[labels > 0])
      keep <- as.integer(names(sizes[sizes >= input$min_cluster_size]))
      labels[!labels %in% keep] <- NA

      dlmo_x <- dlmo_result$ip$inflection_point_fine$x
      dlmo_y <- dlmo_result$ip$inflection_point_fine$y

      col_idx <- which.min(abs(m$x - dlmo_x))
      row_idx <- which.min(abs(m$y - dlmo_y))

      dlmo_cluster <- labels[row_idx, col_idx]

      if (!is.na(dlmo_cluster)) {
        relabels <- setdiff(sort(unique(labels[!is.na(labels)])), dlmo_cluster)
        new_labels <- labels

        new_labels[labels == dlmo_cluster] <- 1
        label_counter <- 2

        for (l in relabels) {
          new_labels[labels == l] <- label_counter
          label_counter <- label_counter + 1
        }

        labels <- new_labels
      }

      list(labels = labels, x = m$x, y = m$y, mat = m$mat)
    })

    # Keep numeric input synchronized with slider
    observeEvent(input$top_percent, {
      updateNumericInput(session, "top_percent_num", value = input$top_percent)
    })

    # Keep slider synchronized with numeric input
    observeEvent(input$top_percent_num, {
      val <- input$top_percent_num
      if (!is.null(val) && val >= 0.1 && val <= 100) {
        updateSliderInput(session, "top_percent", value = val)
      }
    })

    # -------------------------------------------------------------------------
    # Display residual value range represented by the selected percentile
    # -------------------------------------------------------------------------
    output$residual_range <- renderUI({
      m <- get_matrix()
      mat <- m$mat

      if (all(is.na(mat))) return(NULL)

      q <- quantile(mat, input$top_percent / 100, na.rm = TRUE)
      min_val <- min(mat, na.rm = TRUE)

      HTML(paste0(
        "<b>Residual values in lowest ", input$top_percent, "%:</b><br>",
        round(min_val, 3), " to ", round(q, 3),
        " (", sum(mat <= q, na.rm = TRUE), " points)"
      ))
    })

    # -------------------------------------------------------------------------
    # Residual heatmap with cluster outlines and DLMO marker
    # -------------------------------------------------------------------------
    output$residPlot <- renderPlot({
      fc <- full_clusters()

      df <- expand.grid(x = fc$y, y = fc$x)
      df$resid <- as.vector(fc$mat)
      df$cluster <- as.vector(fc$labels)

      # Convert decimal-hour x-axis (time) to POSIXct for plotting
      df$y <- decimal_to_posixct(df$y, dlmo_result$prof$datetime[3])

      # Convex hull polygons per cluster (requires at least 3 points)
      hulls <- df %>%
        filter(!is.na(cluster)) %>%
        group_by(cluster) %>%
        filter(n() >= 3) %>%
        slice(chull(y, x)) %>%
        ungroup()

      # Cluster centroids for label placement
      centroids <- df %>%
        filter(!is.na(cluster)) %>%
        group_by(cluster) %>%
        summarise(
          y = mean(y, na.rm = TRUE),
          x = mean(x, na.rm = TRUE),
          .groups = "drop"
        )

      hulls$cluster <- factor(hulls$cluster)
      centroids$cluster <- factor(centroids$cluster)

      p <- ggplot(df, aes(x = y, y = x)) +
        geom_raster(aes(fill = resid)) +
        geom_polygon(
          data = hulls,
          aes(x = y, y = x, group = cluster, color = factor(cluster)),
          fill = NA,
          linewidth = 1,
          inherit.aes = FALSE
        ) +
        geom_text(
          data = centroids,
          aes(x = y, y = x, label = cluster, color = factor(cluster)),
          inherit.aes = FALSE,
          size = 5,
          fontface = "bold"
        ) +
        geom_point(
          data = data.frame(
            x = decimal_to_posixct(
              dlmo_result$ip$inflection_point_fine$x,
              dlmo_result$prof$datetime[3]
            ),
            y = dlmo_result$ip$inflection_point_fine$y,
            shape_label = "DLMO estimate"
          ),
          aes(x = x, y = y, shape = shape_label),
          color = "black",
          fill = "white",
          stroke = 0.5,
          size = 4
        ) +
        scale_shape_manual(name = "", values = c("DLMO estimate" = 23)) +
        scale_fill_gradientn(
          colours = viridis::viridis(256),
          name = "Residual values",
          limits = range(df$resid, na.rm = TRUE),
          breaks = pretty(range(df$resid, na.rm = TRUE), 5),
          labels = scales::number_format(accuracy = 0.01)
        ) +
        scale_color_brewer(palette = "Set3") +
        labs(
          x = "Time [hh:mm]",
          y = "Melatonin concentration [pg/mL]"
        ) +
        guides(color = "none") +
        theme_minimal(base_size = 14) +
        theme(
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x  = element_text(size = 14),
          axis.text.y  = element_text(size = 14)
        )

      # Export side effect: writes the current plot to disk
      if (!dir.exists("results")) dir.create("results", recursive = TRUE)
      ggsave(file.path("results", "residual_heatmap.svg"), p, width = 7, height = 5, dpi = 300)

      p
    })

    # -------------------------------------------------------------------------
    # Cluster stats: nominal time summaries (clock time)
    # -------------------------------------------------------------------------
    output$clusterNominalStats <- renderTable({
      fc <- full_clusters()
      labels <- fc$labels
      mat <- fc$mat
      x_vals <- fc$x

      total_points <- sum(!is.na(mat))

      clusters <- sort(unique(as.integer(labels[!is.na(labels)])))
      if (length(clusters) == 0) return(NULL)

      stats <- lapply(clusters, function(cl) {
        coords <- which(labels == cl, arr.ind = TRUE)
        times <- x_vals[coords[, 2]]
        size <- length(times)

        data.frame(
          Cluster = cl,
          Size = size,
          "% of grid" = round(100 * size / total_points, 2),
          "Mean time" = format(decimal_to_posixct(mean(as.numeric(times)), dlmo_result$prof$datetime[3]), "%H:%M"),
          "Median time" = format(decimal_to_posixct(median(as.numeric(times)), dlmo_result$prof$datetime[3]), "%H:%M"),
          "Min time" = format(decimal_to_posixct(min(as.numeric(times)), dlmo_result$prof$datetime[3]), "%H:%M"),
          "Max time" = format(decimal_to_posixct(max(as.numeric(times)), dlmo_result$prof$datetime[3]), "%H:%M"),
          check.names = FALSE
        )
      })

      do.call(rbind, stats)
    })

    # -------------------------------------------------------------------------
    # Cluster stats: time differences relative to reference point
    # -------------------------------------------------------------------------
    output$clusterRelativeStats <- renderTable({
      fc <- full_clusters()
      labels <- fc$labels
      mat <- fc$mat
      x_vals <- fc$x

      cluster1_coords <- which(labels == 1, arr.ind = TRUE)
      if (nrow(cluster1_coords) == 0) return(NULL)

      resids <- mat[cluster1_coords]
      dlmo_idx <- which.min(resids)
      dlmo_time <- x_vals[cluster1_coords[dlmo_idx, 2]]

      clusters <- sort(unique(as.integer(labels[!is.na(labels)])))

      stats <- lapply(clusters, function(cl) {
        coords <- which(labels == cl, arr.ind = TRUE)
        times <- x_vals[coords[, 2]]
        diffs_min <- (times - dlmo_time) * 60

        data.frame(
          Cluster = cl,
          "Mean [min]" = round(mean(diffs_min), 2),
          "Median [min]" = round(median(diffs_min), 2),
          "Min [min]" = round(min(diffs_min), 2),
          "Max [min]" = round(max(diffs_min), 2),
          check.names = FALSE
        )
      })

      do.call(rbind, stats)
    })

    # -------------------------------------------------------------------------
    # Violin plot of time offsets by cluster (relative to inferred reference)
    # -------------------------------------------------------------------------
    output$offsetViolinPlot <- renderPlot({
      req(dlmo_result)

      fc <- full_clusters()
      labels <- fc$labels
      mat <- fc$mat
      x_vals <- fc$x

      cluster1_coords <- which(labels == 1, arr.ind = TRUE)
      cluster1_resids <- sapply(seq_len(nrow(cluster1_coords)), function(i) {
        mat[cluster1_coords[i, 1], cluster1_coords[i, 2]]
      })
      min_idx <- which.min(cluster1_resids)
      dlmo_time <- x_vals[cluster1_coords[min_idx, 2]]

      all_coords <- which(!is.na(labels), arr.ind = TRUE)
      all_times <- sapply(seq_len(nrow(all_coords)), function(i) x_vals[all_coords[i, 2]])
      cluster_ids <- sapply(seq_len(nrow(all_coords)), function(i) labels[all_coords[i, 1], all_coords[i, 2]])

      time_diffs <- (all_times - dlmo_time) * 60
      df <- data.frame(TimeDifference = time_diffs, Cluster = factor(cluster_ids))

      p <- ggplot(df, aes(x = TimeDifference, y = Cluster, fill = Cluster)) +
        geom_violin(trim = FALSE) +
        geom_boxplot(width = 0.1, fill = "white") +
        geom_vline(aes(xintercept = 0, color = "DLMO"), linewidth = 1) +
        scale_color_manual(name = NULL, values = c("DLMO" = "black"), guide = guide_legend(order = 1)) +
        scale_fill_brewer(name = "Cluster (mean shown)", palette = "Set3", guide = guide_legend(order = 2)) +
        stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
        scale_x_continuous(
          name = "Time [hh:mm]",
          labels = function(x) {
            dlmo_anchor <- decimal_to_posixct(dlmo_time, dlmo_result$prof$datetime[3])
            format(dlmo_anchor + lubridate::minutes(round(x)), "%H:%M")
          },
          sec.axis = sec_axis(~ ., name = "Time difference relative to reference [minutes]")
        ) +
        labs(y = "Cluster") +
        theme_minimal(base_size = 14) +
        theme(
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x  = element_text(size = 14),
          axis.text.y  = element_text(size = 14)
        )

      # Export side effect: writes the current plot to disk
      if (!dir.exists("results")) dir.create("results", recursive = TRUE)
      ggsave(file.path("results", "offset_violin.svg"), p, width = 6, height = 5, dpi = 300)

      p
    })
  }
)
