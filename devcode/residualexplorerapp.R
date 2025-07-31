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
      if (x < 1 || x > nr || y < 1 || visited[x, y] || !binary_mask[x, y]) next
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
    titlePanel("DLMO residual heatmap explorer"),
    sidebarLayout(
      sidebarPanel(
        sliderInput("top_percent", "Clustering threshold (% lowest residuals):",
                    min = 0.1, max = 100, value = 5, step = 0.5),
        numericInput("min_cluster_size", "Minimum cluster size:", value = 10),
      ),
      mainPanel(
        plotOutput("residPlot"),
        plotOutput("offsetViolinPlot"),
        h3("Cluster Stats (Nominal)"),
        tableOutput("clusterNominalStats"),
        h3("Cluster Stats Relative to DLMO"),
        tableOutput("clusterRelativeStats")
      )
    )
  ),

  server = function(input, output, session) {

    # ---- Convert dlmo_result into matrix form ----
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

    # ---- Cluster labeling ----
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

    observe({
      labs <- sort(unique(as.integer(full_clusters()$labels)))
      labs <- labs[!is.na(labs)]
      updateSelectInput(session, "selected_cluster", choices = labs)
    })

    # ---- Residual heatmap ----
    output$residPlot <- renderPlot({
      fc <- full_clusters()
      df <- expand.grid(x = fc$y, y = fc$x)
      df$resid <- as.vector(fc$mat)
      df$cluster <- as.vector(fc$labels)

      ggplot(df, aes(x = y, y = x, fill = resid)) +
        geom_raster() +
        geom_point(data = df %>% filter(!is.na(cluster)), aes(color = factor(cluster)), size = 1) +
        scale_fill_viridis(option = "magma") +
        labs(fill = "Residual", color = "Cluster") +
        theme_minimal()
    })

    # ---- Nominal cluster stats ----
    output$clusterNominalStats <- renderTable({
      fc <- full_clusters()
      labels <- fc$labels
      mat <- fc$mat
      total_points <- sum(!is.na(mat))

      clusters <- sort(unique(as.integer(labels[!is.na(labels)])))
      if (length(clusters) == 0) return(NULL)

      stats <- lapply(clusters, function(cl) {
        coords <- which(labels == cl, arr.ind = TRUE)
        resids <- mat[coords]
        size <- length(resids)
        data.frame(
          Cluster = cl,
          Size = size,
          PercentOfGrid = round(100 * size / total_points, 2),
          MeanResidual = round(mean(resids), 4),
          MedianResidual = round(median(resids), 4),
          MinResidual = round(min(resids), 4),
          MaxResidual = round(max(resids), 4)
        )
      })

      do.call(rbind, stats)
    })

    # ---- Cluster stats relative to DLMO ----
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
        diffs_hr <- times - dlmo_time
        diffs_min <- diffs_hr * 60
        data.frame(
          Cluster = cl,
          MeanDiffHr = round(mean(diffs_hr), 4),
          MedianDiffHr = round(median(diffs_hr), 4),
          MinDiffHr = round(min(diffs_hr), 4),
          MaxDiffHr = round(max(diffs_hr), 4),
          MeanDiffMin = round(mean(diffs_min), 2),
          MedianDiffMin = round(median(diffs_min), 2)
        )
      })
      do.call(rbind, stats)
    })

    # ---- Violin plot of offsets ----
    output$offsetViolinPlot <- renderPlot({
      req(dlmo_result)
      fc <- full_clusters()
      labels <- fc$labels
      mat <- fc$mat
      x_vals <- fc$x

      cluster1_coords <- which(labels == 1, arr.ind = TRUE)
      cluster1_resids <- sapply(1:nrow(cluster1_coords), function(i) mat[cluster1_coords[i,1], cluster1_coords[i,2]])
      min_idx <- which.min(cluster1_resids)
      dlmo_time <- x_vals[cluster1_coords[min_idx, 2]]

      all_coords <- which(!is.na(labels), arr.ind = TRUE)
      all_times <- sapply(1:nrow(all_coords), function(i) x_vals[all_coords[i, 2]])
      cluster_ids <- sapply(1:nrow(all_coords), function(i) labels[all_coords[i,1], all_coords[i,2]])
      time_diffs <- (all_times - dlmo_time) * 60

      df <- data.frame(TimeDifference = time_diffs, Cluster = factor(cluster_ids))

      ggplot(df, aes(x = TimeDifference, y = Cluster, fill = Cluster)) +
        geom_violin(trim = FALSE) +
        geom_boxplot(width = 0.1, fill = "white") +
        labs(y = "Time Difference (minutes)", x = "Cluster") +
        theme_minimal()
    })
  }
)
