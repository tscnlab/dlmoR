source("/Users/salmathalji/Documents/dlmoR/devcode/bootstrap_modular.R")

# Example: assuming you already have a DLMO result from dlmoR
filename <- system.file("extdata/civibe_melatonin_FD207_day2.csv", package = "dlmoR")
sample_dlmo <- calculate_dlmo(file_path = filename)

# --- Setup output directory ---
output_dir <- "bootstrap_outputs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# --- Bootstrap all methods with error safety + skip existing ---
methods <- list(
  monte_carlo = list(method = "monte_carlo"),
  residual    = list(method = "residual"),
  wild_radem  = list(method = "wild", wild_type = "rademacher"),
  wild_norm   = list(method = "wild", wild_type = "normal"),
  hybrid      = list(method = "hybrid", time_sd = 5/60)
)

results_list <- list()
summary_status <- list()

for (m in names(methods)) {
  rds_file <- file.path(output_dir, paste0("bootstrap_result_", m, ".rds"))
  svg_file <- file.path(output_dir, paste0("bootstrap_plot_", m, ".svg"))

  if (file.exists(rds_file)) {
    cat(sprintf("⏩ Skipping %s (already exists)\n", m))
    summary_status[[m]] <- "skipped (already exists)"
    next
  }

  cat(sprintf("\n🚀 Running %s bootstrap...\n", m))

  tryCatch({
    boot_result <- do.call(dlmo_bootstrap, c(list(sample_dlmo, n_iter = 10), methods[[m]]))
    results_list[[m]] <- boot_result

    saveRDS(boot_result, file = rds_file)
    ggsave(
      svg_file,
      plot = plot_dlmo_bootstrap(boot_result, method_label = m, bw = 300),
      width = 5.5, height = 5, dpi = 300
    )

    cat(sprintf("✅ %s completed and saved.\n", m))
    summary_status[[m]] <- "success"

  }, error = function(e) {
    cat(sprintf("❌ Error in %s bootstrap: %s\n", m, e$message))
    summary_status[[m]] <- paste("failed:", e$message)
  })
}

# --- Final summary report ---
cat("\n📊 Bootstrap Summary Report:\n")
for (m in names(methods)) {
  status <- summary_status[[m]]
  if (grepl("success", status)) {
    cat(sprintf("  ✅ %-12s : %s\n", m, status))
  } else if (grepl("skipped", status)) {
    cat(sprintf("  ⏩ %-12s : %s\n", m, status))
  } else {
    cat(sprintf("  ❌ %-12s : %s\n", m, status))
  }
}
