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
                                               low = TRUE) {
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

  # Update label matrix to remove small clusters
  label_matrix[!filtered_mask] <- 0
  cluster_sizes <- as.integer(cluster_sizes[names(cluster_sizes) %in% filtered_clusters])

  # Re-compute number of valid clusters
  n_clusters <- length(filtered_clusters)

  # Compute centroids and areas
  centroids <- data.frame(cluster = integer(), x = numeric(), y = numeric(), area = numeric())
  for (k in filtered_clusters) {
    indices <- which(label_matrix == k, arr.ind = TRUE)
    if (!is.null(dim(indices)) && nrow(indices) > 0) {
      x_coords <- x_vals[indices[, 2]]
      y_coords <- y_vals[indices[, 1]]
      centroids <- rbind(centroids, data.frame(
        cluster = k,
        x = mean(x_coords),
        y = mean(y_coords),
        area = length(x_coords)
      ))
    }
  }

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

  # Overlay residual contour
  contour(
    x = x_vals,
    y = y_vals,
    z = t(resid_matrix),
    add = TRUE,
    col = "grey50",
    drawlabels = FALSE
  )

  # Overlay cluster outlines using polygons
  for (k in filtered_clusters) {
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
    text(centroid$x, centroid$y, labels = k, col = "red", font = 2, cex = 1.2)
    # ---- Draw reverse gradient ascent arrow ----
    # ---- Draw reverse gradient ascent arrow ----
    i <- which.min(abs(y_vals - centroid$y))
    j <- which.min(abs(x_vals - centroid$x))
    z_current <- resid_matrix[i, j]
    best_z <- z_current
    best_i <- i
    best_j <- j
    for (di in -1:1) {
      for (dj in -1:1) {
        if (di == 0 && dj == 0) next
        ni <- i + di
        nj <- j + dj
        if (ni >= 1 && ni <= nrow(resid_matrix) &&
            nj >= 1 && nj <= ncol(resid_matrix)) {
          z_neighbor <- resid_matrix[ni, nj]
          if (!is.na(z_neighbor) && z_neighbor > best_z) {
            best_z <- z_neighbor
            best_i <- ni
            best_j <- nj
          }
        }
      }
    }
    # If uphill neighbor found, draw arrow
    if (best_i != i || best_j != j) {
      arrows(
        x0 = centroid$x,
        y0 = centroid$y,
        x1 = x_vals[best_j],
        y1 = y_vals[best_i],
        col = "blue",
        length = 0.1,
        lwd = 2
      )
    }


  }

  return(list(
    entropy = entropy,
    modality = n_clusters,
    cluster_sizes = cluster_sizes,
    centroids = centroids,
    top_percent_threshold = threshold
  ))
}

# ---- Run ----
result <- analyze_and_plot_residual_clusters(
  resid_vector = dlmo_result$ip$res_small,
  grid_df = dlmo_result$ip$grid_small,
  top_percent = 0.05,
  min_cluster_size = 10,
  low = TRUE
)

print(result)
