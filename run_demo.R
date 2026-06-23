#!/usr/bin/env Rscript

# ============================================================
# PhyloSOLIDvis - Demo Script
# ============================================================
# This script demonstrates the complete PhyloSOLIDvis pipeline
# using built-in demo data.
# ============================================================

# Load required packages
library(ggplot2)
library(PhyloSOLIDvis)

# Set paths to demo data
demo_path <- system.file("examples/input", package = "PhyloSOLIDvis")
output_path <- "demo_output"
annotation_file <- system.file("examples/annotation.txt", package = "PhyloSOLIDvis")

# Create output directory
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

# Print header
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PhyloSOLIDvis - Demo Run\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Demo data path:", demo_path, "\n")
cat("Output path:", output_path, "\n")
cat("Annotation file:", annotation_file, "\n")
cat(paste(rep("-", 70), collapse = ""), "\n")

# ============================================================
# Method 1: Run complete pipeline with one function
# ============================================================
cat("\n[Method 1] Running complete pipeline (run_all)...\n")
cat(paste(rep("-", 50), collapse = ""), "\n")

results <- run_all(
  inputpath = demo_path,
  outputpath = output_path,
  annotation_file = annotation_file,
  target_mut = "no",
  selected_mutlist = "all",
  manual_fp_file = "no",
  tip_label_offset = 6,
  tip_label_size = 2.5,
  tip_point_size = 0.5,
  heatmap_width = 0.3,
  heatmap_circos_offset = 0.05,
  flipping_point_size = 1.3,
  plot_height = 12,
  plot_width = 18,
  heatmap_colors = c("#D4E8F0", "#7D2224"),
  heatmap_na_color = "white",
  run_adjusted = TRUE,
  adjusted_tip_label_offset = 6,
  adjusted_flipping_point_size = 1.3,
  verbose = TRUE
)

# ============================================================
# Summary
# ============================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("DEMO COMPLETED SUCCESSFULLY!\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

cat("\nGenerated files:\n")
output_files <- list.files(output_path)
for (f in output_files) {
  file_size <- file.info(file.path(output_path, f))$size
  size_str <- if (file_size > 1024 * 1024) {
    paste0(round(file_size / 1024 / 1024, 2), " MB")
  } else if (file_size > 1024) {
    paste0(round(file_size / 1024, 2), " KB")
  } else {
    paste0(file_size, " B")
  }
  cat("  - ", f, " (", size_str, ")\n", sep = "")
}

cat("\nOutput directory:", output_path, "\n")

# ============================================================
# Optional: Show how to use individual functions
# ============================================================
cat("\n", paste(rep("-", 70), collapse = ""), "\n")
cat("Individual functions are also available:\n")
cat("\n")
cat("  # Step 1: Matrix ordering\n")
cat("  order_matrix(inputpath, outputpath)\n")
cat("\n")
cat("  # Step 2: Circular tree plot\n")
cat("  plot_circos(inputpath, outputpath, annotation_file)\n")
cat("\n")
cat("  # Step 3: Heatmap\n")
cat("  plot_heatmap(inputpath, outputpath)\n")
cat("\n")
cat("  # Step 4: Adjusted circos plot\n")
cat("  plot_circos(inputpath, outputpath, annotation_file,\n")
cat("              tip_label_offset = 6, flipping_point_size = 1.3)\n")

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("✅ Done!\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")