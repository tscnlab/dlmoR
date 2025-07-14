# Load required package
if (!requireNamespace("viridis", quietly = TRUE)) {
  install.packages("viridis")
}
library(viridis)

# ---- Helper: 4-connected component cluster labeling ----
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
      x <- loc[1]; y <- loc[2]
      if (x < 1 || x > nr || y < 1 || y > nc || visited[x, y] || !binary_mask[x, y]) next
      visited[x, y] <<- TRUE
      labels[x, y] <<- label_id
      stack <- c(stack, list(c(x - 1, y)), list(c(x + 1, y)), list(c(x, y - 1)), list(c(x, y + 1)))
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

# ---- Main Function: Analyze & Plot ----
analyze_and_plot_residual_clusters <- function(resid_vector, grid_df, top_percent = 0.05) {
  # Prepare grid
  x_vals <- sort(unique(grid_df$x))  # time
  y_vals <- sort(unique(grid_df$y))  # melatonin
  nrow_grid <- length(y_vals)
  ncol_grid <- length(x_vals)

  # Fill matrix
  resid_matrix <- matrix(NA, nrow = nrow_grid, ncol = ncol_grid)
  for (i in seq_along(resid_vector)) {
    row <- which(y_vals == grid_df$y[i])
    col <- which(x_vals == grid_df$x[i])
    resid_matrix[row, col] <- resid_vector[i]
  }

  # Normalize
  resid_matrix <- resid_matrix / sum(resid_matrix, na.rm = TRUE)

  # Entropy
  p <- as.vector(resid_matrix)
  p <- p[!is.na(p) & p > 0]
  entropy <- -sum(p * log(p))

  # Threshold for bottom X%
  threshold <- quantile(resid_matrix, probs = top_percent, na.rm = TRUE)
  low_mask <- resid_matrix <= threshold

  # Label low-residual clusters
  label_info <- label_clusters(low_mask)
  label_matrix <- label_info$label_matrix
  n_clusters <- label_info$n_clusters

  # Compute centroids
  centroids <- data.frame(cluster = integer(), x = numeric(), y = numeric())
  for (k in 1:n_clusters) {
    indices <- which(label_matrix == k, arr.ind = TRUE)
    if (nrow(indices) > 0) {
      x_coords <- x_vals[indices[, 2]]
      y_coords <- y_vals[indices[, 1]]
      centroids <- rbind(centroids, data.frame(
        cluster = k,
        x = mean(x_coords),
        y = mean(y_coords)
      ))
    }
  }

  # --- Plotting ---
  image(
    x_vals,
    y_vals,
    t(resid_matrix),
    col = viridis(100, direction = -1),  # reversed: dark = high residuals
    xlab = "Time (decimal hours)",
    ylab = "Melatonin",
    main = "Residual Surface with Clustered Low Residuals"
  )

  contour(
    x_vals,
    y_vals,
    t(resid_matrix),
    add = TRUE,
    col = "grey50",
    drawlabels = FALSE
  )

  for (k in 1:n_clusters) {
    idx <- which(label_matrix == k, arr.ind = TRUE)
    x_pts <- x_vals[idx[, 2]]
    y_pts <- y_vals[idx[, 1]]
    points(x_pts, y_pts, pch = 20, col = k + 1)
    text(mean(x_pts), mean(y_pts), labels = k, col = "black", font = 2, cex = 0.8)
  }

  # --- Return summary
  return(list(
    entropy = entropy,
    modality = n_clusters,
    top_percent_threshold = threshold,
    centroids = centroids
  ))
}

# ---- Run ----
result <- analyze_and_plot_residual_clusters(
  resid_vector = dlmo_result$ip$res_small,
  grid_df = dlmo_result$ip$grid_small,
  top_percent = 0.05
)

print(result)
