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
    titlePanel("Residual Cluster Viewer (DEBUG MODE)"),
    sidebarLayout(
      sidebarPanel(
        sliderInput("top_percent", "Threshold for clustering (% of residuals):",
                    min = 0.1, max = 20, value = 1, step = 0.5),
        numericInput("min_cluster_size", "Minimum cluster size (pixels):", value = 10),
        selectInput("selected_cluster", "Select cluster:", choices = NULL)
      ),
      mainPanel(
        plotOutput("residPlot"),
        plotOutput("radialPlot"),
        tableOutput("clusterStats")
      )
    )
  ),

  server = function(input, output, session) {

    get_matrix <- reactive({
      validate(need(exists("dlmo_result"), "dlmo_result not found in global environment"))
      grid_df <- dlmo_result$ip$grid_small
      res_vals <- dlmo_result$ip$res_small
      x_vals <- sort(unique(grid_df$x))
      y_vals <- sort(unique(grid_df$y))
      mat <- matrix(NA, nrow = length(y_vals), ncol = length(x_vals))

      for (i in seq_along(res_vals)) {
        x <- grid_df$x[i]
        y <- grid_df$y[i]
        row <- which(y_vals == y)
        col <- which(x_vals == x)
        mat[row, col] <- res_vals[i]
      }

      list(mat = mat, x = x_vals, y = y_vals)
    })

    full_clusters <- reactive({
      m <- get_matrix()
      resid_matrix <- m$mat
      threshold <- quantile(resid_matrix, input$top_percent / 100, na.rm = TRUE)
      mask <- resid_matrix <= threshold
      labeled <- label_clusters(mask)
      label_matrix <- labeled$label_matrix
      sizes <- table(label_matrix[label_matrix > 0])
      keep <- as.integer(names(sizes[sizes >= input$min_cluster_size]))
      label_matrix[!label_matrix %in% keep] <- 0
      label_matrix <- matrix(factor(label_matrix, levels = keep), nrow = nrow(resid_matrix))
      list(labels = label_matrix, x = m$x, y = m$y)
    })

    observe({
      labs <- sort(unique(as.integer(full_clusters()$labels)))
      labs <- labs[!is.na(labs)]
      updateSelectInput(session, "selected_cluster", choices = labs)
    })

    radial_profile <- reactive({
      fc <- full_clusters()  # force dependency
      selected <- input$selected_cluster
      if (is.null(selected) || selected == "") return(NULL)

      m <- get_matrix()
      resid_matrix <- m$mat
      x_vals <- m$x
      y_vals <- m$y

      cluster_mask <- fc$labels
      coords <- which(cluster_mask == as.integer(selected), arr.ind = TRUE)
      if (nrow(coords) < 1) return(NULL)

      # Fetch fresh residuals
      z_vals <- resid_matrix[coords]

      # Get new cluster minimum
      min_idx <- which.min(z_vals)
      min_i <- coords[min_idx, 1]
      min_j <- coords[min_idx, 2]

      # Compute new radial profile around the new center
      radii <- 1:5
      means <- sapply(radii, function(r) {
        ring_mask <- outer(1:nrow(resid_matrix), 1:ncol(resid_matrix), function(i, j) {
          d <- sqrt((i - min_i)^2 + (j - min_j)^2)
          d >= r - 0.5 & d < r + 0.5
        })
        mean(resid_matrix[ring_mask], na.rm = TRUE)
      })

      list(
        radii = radii,
        means = means,
        slope = coef(lm(means ~ radii))[2],
        size = nrow(coords),
        mean_resid = mean(z_vals, na.rm = TRUE),
        min_ij = c(min_i, min_j)
      )
    })


    output$residPlot <- renderPlot({
      m <- get_matrix()
      resid_matrix <- m$mat
      x_vals <- m$x
      y_vals <- m$y
      lbls <- full_clusters()$labels

      image(x_vals, y_vals, t(resid_matrix), col = viridis(100), xlab = "Time", ylab = "Melatonin",
            main = "Residuals + Cluster Boundaries")
      contour(x_vals, y_vals, t(resid_matrix), add = TRUE, col = "white", lwd = 0.7)

      for (k in na.omit(unique(as.integer(lbls)))) {
        indices <- which(lbls == k, arr.ind = TRUE)
        x <- x_vals[indices[, 2]]
        y <- y_vals[indices[, 1]]
        if (length(x) >= 3) {
          h <- chull(x, y)
          lines(c(x[h], x[h[1]]), c(y[h], y[h[1]]), col = "red", lty = 2, lwd = 2)
        }
      }
    })

    output$radialPlot <- renderPlot({
      rp <- radial_profile()
      if (is.null(rp)) return(NULL)
      plot(rp$radii, rp$means, type = "b", pch = 16, col = "blue",
           xlab = "Radius (pixels)", ylab = "Mean residual",
           main = paste("Radial Profile for Cluster", input$selected_cluster))
      abline(lm(rp$means ~ rp$radii), col = "red", lty = 2)
    })

    output$clusterStats <- renderTable({
      rp <- radial_profile()
      if (is.null(rp)) return(NULL)
      data.frame(
        Cluster = input$selected_cluster,
        Pixels = rp$size,
        MeanResidual = round(rp$mean_resid, 3),
        RadialSlope = round(rp$slope, 4)
      )
    })
  }
)
