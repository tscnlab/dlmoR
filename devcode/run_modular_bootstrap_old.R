source("/Users/salmathalji/Documents/dlmoR/devcode/bootstrap_modular.R")

# Example: assuming you already have a DLMO result from dlmoR
filename <- system.file("extdata/civibe_melatonin_FD207_day2.csv", package = "dlmoR")
sample_dlmo <- calculate_dlmo(file_path = filename)

# # Run bootstrap with residual resampling
# # boot_result <- dlmo_bootstrap(sample_dlmo, method = "residual", n_iter = 300)
# # boot_result <- dlmo_bootstrap(sample_dlmo, method = "wild", wild_type = "rademacher", n_iter = 300)
# boot_result <- dlmo_bootstrap(sample_dlmo, method = "hybrid", time_sd = 5/60, n_iter = 300)
#
# # Inspect results
# boot_result$mean
# boot_result$ci
# boot_result$bootstrap_values[1:5]
#
# # Plot results
# plot_dlmo_bootstrap(boot_result, method_label = "Hybrid Bootstrap", bw = 300)

# --- Bootstrap all methods with error safety + summary ---
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
  cat(sprintf("\n🚀 Running %s bootstrap...\n", m))

  tryCatch({
    boot_result <- do.call(dlmo_bootstrap, c(list(sample_dlmo, n_iter = 300), methods[[m]]))
    results_list[[m]] <- boot_result

    saveRDS(boot_result, file = paste0("bootstrap_result_", m, ".rds"))
    ggsave(
      paste0("bootstrap_plot_", m, ".svg"),
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
  } else {
    cat(sprintf("  ❌ %-12s : %s\n", m, status))
  }
}
