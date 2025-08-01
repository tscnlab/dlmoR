library(shiny)
library(viridis)
library(ggplot2)
library(dplyr)
library(tidyr)
# ---- Label clusters (4-connected) ----
label_clusters <- function(binary_mask) {
  nr <- nrow(binary_mask)
  nc <- ncol(binary_mask)
  visited <- matrix(FALSE, nr, nc)
  labels <- matrix(0, nr, nc)
  label_id <- 0

  flood_fill <- function(i, j, label_id) {
    stack <- list(c(i, j))
    while (length(stack) > 0) {
      loc <- stack[[1]]; stack <- stack[-1]
      x <- loc[1]; y <- loc[2]
      if (x < 1 || x > nr || y < 1 || y > nc || visited[x, y] || !binary_mask[x, y]) next
      visited[x, y] <<- TRUE
      labels[x, y] <<- label_id
      stack <- c(stack, list(c(x-1, y)), list(c(x+1, y)), list(c(x, y-1)), list(c(x, y+1)))
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
  return(list(label_matrix = labels, n_clusters = label_id))
}

# ---- Helper to convert decimal hours to HH:MM:SS ----
decimal_to_hms <- function(time_decimal) {
  h <- floor(time_decimal)
  m <- floor((time_decimal - h) * 60)
  s <- round((((time_decimal - h) * 60) - m) * 60)
  sprintf("%02d:%02d:%02d", h, m, s)
}

# ---- Shiny App ----
shinyApp(
  ui = fluidPage(
    titlePanel("DLMO Residual Cluster Viewer"),
    sidebarLayout(
      sidebarPanel(
        sliderInput("top_percent", "Clustering threshold (% lowest residuals):",
                    min = 0.1, max = 100, value = 5, step = 0.5),
        numericInput("min_cluster_size", "Minimum cluster size:", value = 10),
        selectInput("selected_cluster", "Select cluster:", choices = NULL)
      ),
      mainPanel(
        plotOutput("residPlot"),
        plotOutput("offsetViolinPlot"),
        tableOutput("offsetSummaryTable"),
        tableOutput("offsetSummaryByCluster")

      )

    )
  ),

  server = function(input, output, session) {

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

      list(labels = labels, x = m$x, y = m$y)
    })

    observe({
      labs <- sort(unique(as.integer(full_clusters()$labels)))
      labs <- labs[!is.na(labs)]
      updateSelectInput(session, "selected_cluster", choices = labs)
    })

    cluster_time_diff_stats <- reactive({
      m <- get_matrix()
      mat <- m$mat
      x_vals <- m$x
      labels <- full_clusters()$labels

      # Get all coordinates in cluster 1
      cluster1_coords <- which(labels == 1, arr.ind = TRUE)
      if (nrow(cluster1_coords) == 0) return(NULL)

      # Find point in cluster 1 with the minimum residual
      cluster1_residuals <- apply(cluster1_coords, 1, function(row) mat[row[1], row[2]])
      min_idx <- which.min(cluster1_residuals)
      dlmo_j <- cluster1_coords[min_idx, 2]
      dlmo_time <- x_vals[dlmo_j]

      # Time difference from DLMO to all points in any cluster
      all_coords <- which(!is.na(labels), arr.ind = TRUE)
      all_times <- x_vals[all_coords[, 2]]
      time_diffs <- all_times - dlmo_time

      # Summary statistics
      data.frame(
        Statistic = c("Mean", "Median", "SD", "Min", "Max"),
        `Time Difference (hrs)` = round(c(
          mean(time_diffs),
          median(time_diffs),
          sd(time_diffs),
          min(time_diffs),
          max(time_diffs)
        ), 4)
      )
    })

    output$offsetSummaryTable <- renderTable({
      req(dlmo_result)
      m <- get_matrix()
      labels <- full_clusters()$labels
      mat <- m$mat
      x_vals <- m$x

      cluster1_coords <- which(labels == 1, arr.ind = TRUE)
      cluster1_resids <- sapply(1:nrow(cluster1_coords), function(i) mat[cluster1_coords[i,1], cluster1_coords[i,2]])
      min_idx <- which.min(cluster1_resids)
      dlmo_time <- x_vals[cluster1_coords[min_idx, 2]]

      all_coords <- which(!is.na(labels), arr.ind = TRUE)
      all_times <- sapply(1:nrow(all_coords), function(i) x_vals[all_coords[i, 2]])
      time_diffs_hr <- all_times - dlmo_time
      time_diffs_min <- time_diffs_hr * 60

      data.frame(
        Statistic = c("Mean", "Median", "Standard Deviation", "Minimum", "Maximum"),
        `Time Difference (hrs)` = round(c(
          mean(time_diffs_hr),
          median(time_diffs_hr),
          sd(time_diffs_hr),
          min(time_diffs_hr),
          max(time_diffs_hr)
        ), 4),
        `Time Difference (min)` = round(c(
          mean(time_diffs_min),
          median(time_diffs_min),
          sd(time_diffs_min),
          min(time_diffs_min),
          max(time_diffs_min)
        ), 2),
        check.names = FALSE
      )
    })


    output$offsetViolinPlot <- renderPlot({
      req(dlmo_result)
      m <- get_matrix()
      labels <- full_clusters()$labels
      mat <- m$mat
      x_vals <- m$x

      cluster1_coords <- which(labels == 1, arr.ind = TRUE)
      cluster1_resids <- sapply(1:nrow(cluster1_coords), function(i) mat[cluster1_coords[i,1], cluster1_coords[i,2]])
      min_idx <- which.min(cluster1_resids)
      dlmo_time <- x_vals[cluster1_coords[min_idx, 2]]

      all_coords <- which(!is.na(labels), arr.ind = TRUE)
      all_times <- sapply(1:nrow(all_coords), function(i) x_vals[all_coords[i, 2]])
      cluster_ids <- sapply(1:nrow(all_coords), function(i) labels[all_coords[i,1], all_coords[i,2]])
      time_diffs <- (all_times - dlmo_time) * 60

      df <- data.frame(TimeDifference = time_diffs, Cluster = factor(cluster_ids))

      stats_df <- df %>%
        group_by(Cluster) %>%
        summarise(
          mean = mean(TimeDifference),
          median = median(TimeDifference),
          sd = sd(TimeDifference),
          min = min(TimeDifference),
          max = max(TimeDifference)
        )

      # Add a column for shape labels
      stats_long <- stats_df %>%
        pivot_longer(cols = c(mean, median, min, max), names_to = "stat", values_to = "TimeDifference") %>%
        mutate(shape = case_when(
          stat == "mean" ~ "Mean",
          stat == "median" ~ "Median",
          stat == "min" ~ "Min",
          stat == "max" ~ "Max"
        ))

      ggplot(df, aes(x = TimeDifference, y = Cluster)) +
        geom_violin(aes(fill = Cluster), trim = FALSE, alpha = 0.3, color = "black", show.legend = FALSE) +
        geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", alpha = 0.4, show.legend = FALSE) +
        geom_point(data = stats_long, aes(x = TimeDifference, y = Cluster, shape = shape), size = 3, color = "black", inherit.aes = FALSE) +
        geom_errorbar(data = stats_df, aes(y = Cluster, xmin = mean - sd, xmax = mean + sd, linetype = "±1 SD"),
                      width = 0.2, linewidth = 0.7, color = "black", inherit.aes = FALSE) +
        geom_vline(xintercept = 0, color = "cyan3", linewidth = 1) +
        annotate("text", x = .5, y = 1, label = "DLMO point", vjust = -10, color = "black", fontface = "bold", size = 3.5) +
        scale_shape_manual(
          name = "Summary Stats",
          values = c("Mean" = 18, "Median" = 17, "Min" = 3, "Max" = 4)
        ) +
        scale_linetype_manual(
          name = "Summary Stats",
          values = c("±1 SD" = "dashed")
        ) +
        labs(title = "Time Differences from DLMO to All Cluster Points",
             x = "Time Difference (minutes)", y = "Cluster") +
        theme_minimal() +
        theme(legend.position = "right")
    })





    output$offsetSummaryByCluster <- renderTable({
      req(dlmo_result)
      m <- get_matrix()
      labels <- full_clusters()$labels
      mat <- m$mat
      x_vals <- m$x

      cluster1_coords <- which(labels == 1, arr.ind = TRUE)
      cluster1_resids <- sapply(1:nrow(cluster1_coords), function(i) mat[cluster1_coords[i,1], cluster1_coords[i,2]])
      min_idx <- which.min(cluster1_resids)
      dlmo_time <- x_vals[cluster1_coords[min_idx, 2]]

      all_coords <- which(!is.na(labels), arr.ind = TRUE)
      all_times <- sapply(1:nrow(all_coords), function(i) x_vals[all_coords[i, 2]])
      cluster_ids <- sapply(1:nrow(all_coords), function(i) labels[all_coords[i,1], all_coords[i,2]])
      time_diffs_min <- (all_times - dlmo_time) * 60

      df <- data.frame(Cluster = factor(cluster_ids), Minutes = time_diffs_min)

      df_long <- df %>%
        group_by(Cluster) %>%
        summarise(
          Mean = round(mean(Minutes), 2),
          Median = round(median(Minutes), 2),
          `Standard Deviation` = round(sd(Minutes), 2),
          Minimum = round(min(Minutes), 2),
          Maximum = round(max(Minutes), 2)
        ) %>%
        pivot_longer(-Cluster, names_to = "Statistic", values_to = "Time Difference (min)") %>%
        pivot_wider(names_from = Cluster, values_from = `Time Difference (min)`)

      df_long
    })



    radial_profile <- reactive({
      req(input$selected_cluster)
      fc <- full_clusters()
      selected <- as.integer(input$selected_cluster)
      if (is.na(selected)) return(NULL)

      labels <- fc$labels
      coords <- which(labels == selected, arr.ind = TRUE)
      if (nrow(coords) < 1) return(NULL)

      mat <- get_matrix()$mat
      x_vals <- get_matrix()$x
      z_vals <- apply(coords, 1, function(row) mat[row[1], row[2]])
      min_idx <- which.min(z_vals)
      min_i <- coords[min_idx, 1]
      min_j <- coords[min_idx, 2]

      dists <- sqrt((coords[,1] - min_i)^2 + (coords[,2] - min_j)^2)
      max_radius <- floor(max(dists))
      radii <- 0:max_radius

      means <- sapply(radii, function(r) {
        ring_mask <- outer(1:nrow(mat), 1:ncol(mat), function(i, j) {
          d <- sqrt((i - min_i)^2 + (j - min_j)^2)
          d >= r - 0.5 & d < r + 0.5
        })
        mean(mat[ring_mask], na.rm = TRUE)
      })

      cluster_times <- sapply(1:nrow(coords), function(i) x_vals[coords[i, 2]])
      residuals <- sapply(1:nrow(coords), function(i) mat[coords[i, 1], coords[i, 2]])
      min_resid <- min(residuals, na.rm = TRUE)
      in_range <- which(residuals <= min_resid + 0.1 * abs(min_resid))
      time_inrange <- cluster_times[in_range]
      time_bounds <- range(time_inrange, na.rm = TRUE)
      time_width <- diff(time_bounds)
      time_center <- x_vals[min_j]
      offset_lower <- round(time_bounds[1] - time_center, 3)
      offset_upper <- round(time_bounds[2] - time_center, 3)

      instantaneous_slopes <- c(NA, diff(means))
      laplacian <- c(NA, diff(means, differences = 2), NA)

      label_prefix <- ifelse(selected == 1, "DLMO", "Minimum point")

      # Contrast = mean(edge residuals) - min(residual in cluster)
      edge_vals <- apply(coords, 1, function(row) {
        i <- row[1]; j <- row[2]
        neighbors <- expand.grid(i = (i-1):(i+1), j = (j-1):(j+1))
        neighbors <- subset(neighbors, i >= 1 & i <= nrow(mat) & j >= 1 & j <= ncol(mat))
        any(sapply(1:nrow(neighbors), function(k) {
          ni <- neighbors[k, 1]; nj <- neighbors[k, 2]
          is.na(labels[ni, nj]) || labels[ni, nj] != selected
        }))
      })
      outer_vals <- z_vals[edge_vals]
      contrast_val <- if (length(outer_vals) > 0) mean(outer_vals, na.rm = TRUE) - z_vals[min_idx] else NA


      list(
        radii = radii,
        means = means,
        slope = coef(lm(means ~ radii))[2],
        slope_inner = ifelse(length(radii) > 3, coef(lm(means[1:3] ~ radii[1:3]))[2], NA),
        contrast = contrast_val,
        unique_score = 1 / max(1, length(which(diff(sign(diff(means))) == -2))),
        slice_width_10 = sum(abs(means - means[1]) <= 0.1 * abs(means[1])),
        inst_slopes = instantaneous_slopes,
        laplacian = laplacian,
        min_val = means[1],
        time_bounds = time_bounds,
        time_width = time_width,
        offset_lower = offset_lower,
        offset_upper = offset_upper,
        dlmo_time = time_center,
        dlmo_time_hms = decimal_to_hms(time_center),
        bound_left_hms = decimal_to_hms(time_bounds[1]),
        bound_right_hms = decimal_to_hms(time_bounds[2]),
        offset_lower_hms = decimal_to_hms(time_center + offset_lower),
        offset_upper_hms = decimal_to_hms(time_center + offset_upper),
        size = nrow(coords),
        mean_resid = mean(z_vals, na.rm = TRUE),
        label_prefix = label_prefix
      )
    })

    output$timeDiffStats <- renderTable({
      cluster_time_diff_stats()
    })

    output$residPlot <- renderPlot({
      req(input$selected_cluster)
      m <- get_matrix()
      mat <- m$mat
      x <- m$x
      y <- m$y
      lbls <- full_clusters()$labels

      image(x, y, t(mat), col = viridis(100), xlab = "Time (decimal hours)", ylab = "Melatonin",
            main = "Residuals + Clusters")
      contour(x, y, t(mat), add = TRUE, col = "white", lwd = 0.5)

      for (k in sort(unique(na.omit(as.integer(lbls))))) {
        coords <- which(lbls == k, arr.ind = TRUE)
        cx <- x[coords[, 2]]
        cy <- y[coords[, 1]]
        if (length(cx) >= 3) {
          h <- chull(cx, cy)
          lines(c(cx[h], cx[h[1]]), c(cy[h], cy[h[1]]), col = "red", lty = 2, lwd = 1.5)
          text(mean(cx), mean(cy), labels = k, col = "white", font = 2)
        }
      }

      dlmo_x <- dlmo_result$ip$inflection_point_fine$x
      dlmo_y <- dlmo_result$ip$inflection_point_fine$y
      points(dlmo_x, dlmo_y, col = "cyan", pch = 18, cex = 2.5)
    })

    output$radialPlot <- renderPlot({
      rp <- radial_profile()
      if (is.null(rp)) return()
      plot(rp$radii, rp$means, type = "b", pch = 16, col = "blue",
           main = paste("Radial Profile: Cluster", input$selected_cluster),
           xlab = "Radius", ylab = "Mean residual")
      lines(rp$radii, rp$inst_slopes, col = "orange", lty = 3)
      if (as.integer(input$selected_cluster) == 1) {
        points(0, rp$means[1], col = "cyan", pch = 18, cex = 2)
        text(1, rp$means[1], labels = paste0("DLMO residual = ", round(rp$means[1], 3)), col = "cyan", pos = 4)
      }
    })

    output$instSlopePlot <- renderPlot({
      rp <- radial_profile()
      if (is.null(rp)) return()
      plot(rp$radii, rp$inst_slopes, type = "l", lwd = 2, col = "darkgreen",
           main = "Instantaneous Slope of Radial Profile",
           xlab = "Radius", ylab = "Instantaneous Slope")
      abline(h = 0, col = "gray", lty = 2)
    })

    output$laplacianPlot <- renderPlot({
      rp <- radial_profile()
      if (is.null(rp)) return()
      plot(rp$radii, rp$laplacian, type = "l", lwd = 2, col = "purple",
           main = "Laplacian of Radial Profile",
           xlab = "Radius", ylab = "Second Derivative (Laplacian)")
      abline(h = 0, col = "gray", lty = 2)
    })

    output$metricsTable <- renderTable({
      rp <- radial_profile()
      if (is.null(rp)) return()
      label <- rp$label_prefix
      data.frame(
        Metric = c("Pixels in Cluster", "Mean Residual", "Radial Slope (full)", "Radial Slope (inner)",
                   "Contrast (outer - center)", "Uniqueness Score", "Slice Width @10%",
                   paste(label, "Estimate (Time)"), paste(label, "Time (HH:MM:SS)"),
                   "Left Bound of Time Near Min", "Right Bound of Time Near Min",
                   "Left Bound (HH:MM:SS)", "Right Bound (HH:MM:SS)",
                   "Offset from Center (lower)", "Offset from Center (upper)",
                   "Offset Lower (HH:MM:SS)", "Offset Upper (HH:MM:SS)"),
        Value = c(round(rp$size, 3), round(rp$mean_resid, 3), round(rp$slope, 3), round(rp$slope_inner, 3),
                  round(rp$contrast, 3), round(rp$unique_score, 3), round(rp$slice_width_10, 3),
                  round(rp$dlmo_time, 3), rp$dlmo_time_hms,
                  round(rp$time_bounds[1], 3), round(rp$time_bounds[2], 3),
                  rp$bound_left_hms, rp$bound_right_hms,
                  rp$offset_lower, rp$offset_upper,
                  rp$offset_lower_hms, rp$offset_upper_hms),
        Interpretation = c("How many pixels in cluster", "Average residual across cluster",
                           "Overall radial rise from center", "How sharp the center is (steepness)",
                           "Residual difference between cluster edges and minimum",
                           "1 / number of local minima in full residuals",
                           "# of timepoints near minimum within 10% error",
                           paste("Time of", label), paste("Formatted", label),
                           "Earliest time with residual within 10% of minimum",
                           "Latest time with residual within 10% of minimum",
                           "Formatted lower bound time", "Formatted upper bound time",
                           "Lower time offset from center (hrs)",
                           "Upper time offset from center (hrs)",
                           "Formatted lower offset", "Formatted upper offset"),
        stringsAsFactors = FALSE
      )
    })
  }
)
