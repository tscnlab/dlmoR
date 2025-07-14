library(shiny)
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

# ---- Shiny App ----
shinyApp(
  ui = fluidPage(
    titlePanel("Residual Cluster Viewer"),
    sidebarLayout(
      sidebarPanel(
        sliderInput("top_percent", "Clustering threshold (% lowest residuals):",
                    min = 0.1, max = 100, value = 5, step = 0.5),
        numericInput("min_cluster_size", "Minimum cluster size:", value = 10),
        selectInput("selected_cluster", "Select cluster:", choices = NULL)
      ),
      mainPanel(
        plotOutput("residPlot"),
        plotOutput("radialPlot"),
        plotOutput("directionalPlot"),
        tableOutput("clusterStats")
      )
    )
  ),

  server = function(input, output, session) {

    get_matrix <- reactive({
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
      m <- get_matrix()
      mat <- m$mat
      threshold <- quantile(mat, input$top_percent / 100, na.rm = TRUE)
      mask <- mat <= threshold
      labeled <- label_clusters(mask)
      labels <- labeled$label_matrix
      sizes <- table(labels[labels > 0])
      keep <- as.integer(names(sizes[sizes >= input$min_cluster_size]))
      labels[!labels %in% keep] <- 0
      labels <- matrix(factor(labels, levels = keep), nrow = nrow(mat))
      list(labels = labels, x = m$x, y = m$y)
    })

    observe({
      labs <- sort(unique(as.integer(full_clusters()$labels)))
      labs <- labs[!is.na(labs)]
      updateSelectInput(session, "selected_cluster", choices = labs)
    })

    radial_profile <- reactive({
      fc <- full_clusters()
      selected <- as.integer(input$selected_cluster)
      if (is.null(selected)) return(NULL)

      labels <- fc$labels
      coords <- which(labels == selected, arr.ind = TRUE)
      if (nrow(coords) < 1) return(NULL)

      mat <- get_matrix()$mat
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

      list(
        radii = radii,
        means = means,
        slope = coef(lm(means ~ radii))[2],
        size = nrow(coords),
        mean_resid = mean(z_vals, na.rm = TRUE)
      )
    })

    directional_slices <- reactive({
      fc <- full_clusters()
      selected <- as.integer(input$selected_cluster)
      if (is.null(selected)) return(NULL)

      labels <- fc$labels
      coords <- which(labels == selected, arr.ind = TRUE)
      if (nrow(coords) < 1) return(NULL)

      m <- get_matrix()
      mat <- m$mat
      x_vals <- m$x
      y_vals <- m$y

      z_vals <- apply(coords, 1, function(row) mat[row[1], row[2]])
      min_idx <- which.min(z_vals)
      ci <- coords[min_idx, 1]
      cj <- coords[min_idx, 2]

      horizontal <- mat[ci, ]
      vertical <- mat[, cj]

      diag_idx <- -min(ci - 1, cj - 1):(min(nrow(mat) - ci, ncol(mat) - cj))
      diag_vals <- sapply(diag_idx, function(d) {
        i <- ci + d
        j <- cj + d
        if (i >= 1 && i <= nrow(mat) && j >= 1 && j <= ncol(mat)) mat[i, j] else NA
      })

      list(
        x_vals = x_vals,
        y_vals = y_vals,
        time_slice = horizontal,
        mel_slice = vertical,
        diag_vals = diag_vals,
        diag_x = diag_idx,
        center_i = ci,
        center_j = cj
      )
    })

    output$residPlot <- renderPlot({
      m <- get_matrix()
      mat <- m$mat
      x <- m$x
      y <- m$y
      lbls <- full_clusters()$labels

      image(x, y, t(mat), col = viridis(100), xlab = "Time", ylab = "Melatonin",
            main = "Residuals + Clusters")
      contour(x, y, t(mat), add = TRUE, col = "white", lwd = 0.5)

      for (k in na.omit(unique(as.integer(lbls)))) {
        coords <- which(lbls == k, arr.ind = TRUE)
        cx <- x[coords[, 2]]
        cy <- y[coords[, 1]]
        if (length(cx) >= 3) {
          h <- chull(cx, cy)
          lines(c(cx[h], cx[h[1]]), c(cy[h], cy[h[1]]), col = "red", lty = 2, lwd = 1.5)
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
      #abline(lm(rp$means ~ rp$radii), col = "red", lty = 2)
      points(0, rp$means[1], col = "cyan", pch = 18, cex = 2)
      text(1, rp$means[1], labels = paste0("DLMO residual = ", round(rp$means[1], 3)), col = "cyan", pos = 4)
    })

    output$directionalPlot <- renderPlot({
      ds <- directional_slices()
      if (is.null(ds)) return()
      ylim_range <- range(c(ds$time_slice, ds$mel_slice, ds$diag_vals), na.rm = TRUE)
      plot(ds$x_vals, ds$time_slice, type = "l", col = "red", lwd = 2,
           xlab = "Index", ylab = "Residual",
           main = "Directional Slices through Cluster Center",
           ylim = ylim_range)
      #lines(ds$y_vals, ds$mel_slice, col = "green", lwd = 2)
      #lines(ds$diag_x, ds$diag_vals, col = "blue", lwd = 2)
      legend("topright", legend = c("Horizontal (time)", "Vertical (melatonin)", "Diagonal"),
             col = c("red", "green", "blue"), lwd = 2, bty = "n")
    })

    output$clusterStats <- renderTable({
      rp <- radial_profile()
      if (is.null(rp)) return()
      data.frame(
        Pixels = rp$size,
        MeanResidual = round(rp$mean_resid, 4),
        RadialSlope = round(rp$slope, 4)
      )
    })
  }
)
