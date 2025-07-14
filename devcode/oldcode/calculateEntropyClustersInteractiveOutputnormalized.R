# Load viridis for perceptual colormap
if (!requireNamespace("viridis", quietly = TRUE)) install.packages("viridis")
library(viridis)

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

# ---- Main analysis and plot function ----
analyze_and_plot_residual_clusters <- function(resid_vector, grid_df,
                                               top_percent = 0.05,
                                               min_cluster_size = 10,
                                               low = TRUE,
                                               dlmo_x = NULL,
                                               dlmo_y = NULL) {
  x_vals <- sort(unique(grid_df$x))
  y_vals <- sort(unique(grid_df$y))
  nrow_grid <- length(y_vals)
  ncol_grid <- length(x_vals)

  # Build residual matrix (row = y, col = x)
  resid_matrix <- matrix(NA, nrow = nrow_grid, ncol = ncol_grid)
  for (i in seq_along(resid_vector)) {
    row <- match(grid_df$y[i], y_vals)
    col <- match(grid_df$x[i], x_vals)
    resid_matrix[row, col] <- resid_vector[i]
  }

  # Normalize
  resid_matrix <- resid_matrix / sum(resid_matrix, na.rm = TRUE)

  # Entropy
  p <- as.vector(resid_matrix)
  p <- p[!is.na(p) & p > 0]
  entropy <- -sum(p * log(p))
  total_area <- sum(!is.na(resid_matrix))
  max_entropy <- log(total_area)
  normalized_entropy <- entropy / max_entropy

  # Thresholding
  threshold <- quantile(resid_matrix, if (low) top_percent else 1 - top_percent, na.rm = TRUE)
  high_mask <- if (low) resid_matrix <= threshold else resid_matrix >= threshold

  # Label clusters
  label_info <- label_clusters(high_mask)
  label_matrix <- label_info$label_matrix

  # Count pixels in each cluster
  cluster_sizes <- table(label_matrix[label_matrix > 0])

  # Filter clusters by size
  filtered_clusters <- as.integer(names(cluster_sizes[cluster_sizes >= min_cluster_size]))
  filtered_mask <- matrix(label_matrix %in% filtered_clusters, nrow = nrow_grid, ncol = ncol_grid)

  # Update label matrix to remove small clusters and reindex
  label_matrix[!filtered_mask] <- 0
  cluster_sizes <- as.integer(cluster_sizes[names(cluster_sizes) %in% filtered_clusters])
  n_clusters <- length(filtered_clusters)

  # Reindex labels consecutively
  new_label_matrix <- matrix(0, nrow = nrow_grid, ncol = ncol_grid)
  for (new_id in seq_along(filtered_clusters)) {
    new_label_matrix[label_matrix == filtered_clusters[new_id]] <- new_id
  }
  label_matrix <- new_label_matrix

  # Compute centroids and areas
  centroids <- data.frame(cluster = integer(), decimal_hours = numeric(), melatonin = numeric(), area = numeric(), pct_area = numeric())
  for (k in seq_len(n_clusters)) {
    indices <- which(label_matrix == k, arr.ind = TRUE)
    if (!is.null(dim(indices)) && nrow(indices) > 0) {
      x_coords <- x_vals[indices[, 2]]
      y_coords <- y_vals[indices[, 1]]
      n_pixels <- length(x_coords)
      centroids <- rbind(centroids, data.frame(
        cluster = k,
        decimal_hours = mean(x_coords),
        melatonin = mean(y_coords),
        area = n_pixels,
        pct_area = round(100 * n_pixels / total_area, 2)
      ))
    }
  }

  # Use correct DLMO coordinates for plotting
  dlmo_x <- dlmo_result$ip$inflection_point_fine$x
  dlmo_y <- dlmo_result$ip$inflection_point_fine$y

  # ---- Plot heatmap ----
  image(
    x = x_vals,
    y = y_vals,
    z = t(resid_matrix),
    col = viridis(100, option = "D"),
    xlab = "Time (decimal hours)",
    ylab = "Melatonin",
    main = if (low) "Residual Surface with Clustered Low Residuals" else "High Residual Clusters"
  )

  contour(
    x = x_vals,
    y = y_vals,
    z = t(resid_matrix),
    add = TRUE,
    col = "grey50",
    drawlabels = FALSE
  )

  for (k in seq_len(n_clusters)) {
    mask <- label_matrix == k
    indices <- which(mask, arr.ind = TRUE)
    x_coords <- x_vals[indices[, 2]]
    y_coords <- y_vals[indices[, 1]]
    if (length(x_coords) >= 3) {
      hull <- chull(x_coords, y_coords)
      lines(
        c(x_coords[hull], x_coords[hull[1]]),
        c(y_coords[hull], y_coords[hull[1]]),
        col = "red", lwd = 2, lty = 2
      )
    }
    centroid <- centroids[centroids$cluster == k, ]
    text(centroid$decimal_hours, centroid$melatonin, labels = k, col = "red", font = 2, cex = 1.2)
  }

  # Add DLMO marker and label offset 45 degrees (top right)
  text_x <- dlmo_x + 0.1
  text_y <- dlmo_y + 0.1
  points(dlmo_x, dlmo_y, col = "cyan", pch = 18, cex = 2.5)  # filled diamond
  # text(
  #   text_x,
  #   text_y,
  #   labels = paste0("DLMO\n", round(dlmo_x, 2), "h, ", format(dlmo_result$dlmo$fine$time)),
  #   col = "cyan",
  #   font = 2,
  #   cex = 1.2
  # )

  # Print summary to console
  cat("\n---- Residual Cluster Summary ----\n")
  cat("Entropy:", round(entropy, 6), "(normalized:", round(normalized_entropy, 6), ")\n")
  cat("Modality (clusters):", n_clusters, "\n")
  cat("Top", round(top_percent * 100, 2), "% threshold value:", signif(threshold, 3), "\n")
  cat("DLMO time:", format(dlmo_result$dlmo$fine$time), "\n")
  cat("DLMO melatonin:", round(dlmo_y, 3), "\n")
  cat("\nCluster centroids:\n")
  print(centroids)

  return(list(
    entropy = entropy,
    normalized_entropy = normalized_entropy,
    modality = n_clusters,
    cluster_sizes = cluster_sizes,
    centroids = centroids,
    top_percent_threshold = threshold,
    dlmo_x = dlmo_x,
    dlmo_y = dlmo_y,
    dlmo_time_str = dlmo_result$dlmo$fine$time
  ))
}

# ---- Launch Shiny App if interactive ----
if (interactive()) {
  if (!requireNamespace("shiny", quietly = TRUE)) install.packages("shiny")
  library(shiny)

  shinyApp(
    ui = fluidPage(
      titlePanel("Interactive Residual Cluster Viewer"),
      sidebarLayout(
        sidebarPanel(
          sliderInput("top_percent", "Threshold for clustering (% of residuals):",
                      min = 0.1, max = 20, value = 5, step = 0.5),
          numericInput("min_cluster_size", "Minimum cluster size (pixels):", value = 10),
          helpText("Thresholding selects the lowest X% of residuals to identify clusters.")
        ),
        mainPanel(
          plotOutput("residPlot"),
          uiOutput("residStats")
        )
      )
    ),
    server = function(input, output) {
      result <- reactive({
        analyze_and_plot_residual_clusters(
          resid_vector = dlmo_result$ip$res_small,
          grid_df = dlmo_result$ip$grid_small,
          top_percent = input$top_percent / 100,
          min_cluster_size = input$min_cluster_size,
          low = TRUE
        )
      })

      output$residPlot <- renderPlot({ result() })

      output$residStats <- renderUI({
        r <- result()
        tagList(
          h4("Cluster Summary"),
          tags$p(strong("Entropy (spread):"), round(r$entropy, 6),
                 em(paste("(normalized:", round(r$normalized_entropy, 6), ") - higher = more diffuse residuals"))),
          tags$p(strong("Modality (clusters):"), r$modality,
                 em("- number of detected clusters")),
          tags$p(strong("Threshold value (residual):"), signif(r$top_percent_threshold, 3),
                 em("- cutoff for lowest X% of values")),
          tags$p(strong("DLMO estimate:"), paste0("Time = ", r$dlmo_time_str, ", Melatonin = ", round(r$dlmo_y, 2))),
          tags$hr(),
          h5("Cluster Centroids and Sizes"),
          renderTable(r$centroids)
        )
      })
    }
  )
}
