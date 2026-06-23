#' Run complete PhyloSOLIDvis pipeline
#'
#' This function runs the complete visualization pipeline in one go:
#' 1. Matrix ordering (order_matrix)
#' 2. Circular tree plot (plot_circos)
#' 3. Heatmap (plot_heatmap)
#' 4. Optional: Adjusted circular tree plot with custom parameters
#'
#' @param inputpath Path to PhyloSOLID phylo directory
#' @param outputpath Path to output directory
#' @param annotation_file Path to annotation file (optional)
#' @param target_mut Target mutation to highlight (default: "no")
#' @param selected_mutlist Selected mutations (default: "all")
#' @param manual_fp_file Manual false positive file (default: "no")
#' @param tip_label_offset Offset for tip labels (default: 10)
#' @param tip_label_size Size of tip labels (default: 2.5)
#' @param tip_point_size Size of tip points (default: 0.5)
#' @param heatmap_width Width of heatmap (default: 0.3)
#' @param heatmap_circos_offset Offset for heatmap (default: 0.04)
#' @param flipping_point_size Size of flipping points (default: 0.2)
#' @param plot_height Plot height (default: 12)
#' @param plot_width Plot width (default: 18)
#' @param heatmap_colors Colors for heatmap (default: c("#D4E8F0", "#7D2224"))
#' @param heatmap_na_color Color for NA values in heatmap (default: "white")
#' @param run_adjusted Logical. Whether to run adjusted plot with different parameters. Default: FALSE
#' @param adjusted_tip_label_offset Adjusted tip label offset (default: 6)
#' @param adjusted_flipping_point_size Adjusted flipping point size (default: 1.3)
#' @param verbose Print progress messages. Default: TRUE
#'
#' @return Invisibly returns a list containing all results
#' @export
#'
#' @examples
#' \dontrun{
#' results <- run_all(
#'   inputpath = "path/to/phylo/",
#'   outputpath = "path/to/output/"
#' )
#' }
run_all <- function(
    inputpath,
    outputpath,
    annotation_file = NULL,
    target_mut = "no",
    selected_mutlist = "all",
    manual_fp_file = "no",
    tip_label_offset = 10,
    tip_label_size = 2.5,
    tip_point_size = 0.5,
    heatmap_width = 0.3,
    heatmap_circos_offset = 0.04,
    flipping_point_size = 0.2,
    plot_height = 12,
    plot_width = 18,
    heatmap_colors = c("#D4E8F0", "#7D2224"),
    heatmap_na_color = "white",
    run_adjusted = FALSE,
    adjusted_tip_label_offset = 6,
    adjusted_flipping_point_size = 1.3,
    verbose = TRUE
) {
  
  # Start timing
  start_time <- Sys.time()
  
  if (verbose) {
    message("\n", paste(rep("=", 70), collapse = ""))
    message("PhyloSOLIDvis - Complete Pipeline")
    message(paste(rep("=", 70), collapse = ""))
    message("Start time: ", start_time)
    message("Input path: ", inputpath)
    message("Output path: ", outputpath)
    if (!is.null(annotation_file)) {
      message("Annotation file: ", annotation_file)
    }
    message(paste(rep("-", 70), collapse = ""))
  }
  
  # Create output directory
  if (!dir.exists(outputpath)) {
    dir.create(outputpath, recursive = TRUE)
    if (verbose) message("Created output directory: ", outputpath)
  }
  
  # Initialize results list
  results <- list()
  
  # ============================================================
  # Step 1: Matrix Ordering
  # ============================================================
  if (verbose) message("\n[Step 1/4] Matrix Ordering...")
  
  results$order <- order_matrix(
    inputpath = inputpath,
    outputpath = outputpath,
    verbose = verbose
  )
  
  if (verbose) message("✓ Matrix ordering completed")
  
  # ============================================================
  # Step 2: Circular Tree Plot
  # ============================================================
  if (verbose) message("\n[Step 2/4] Circular Tree Plot...")
  
  results$circos <- plot_circos(
    inputpath = inputpath,
    outputpath = outputpath,
    annotation_file = annotation_file,
    target_mut = target_mut,
    selected_mutlist = selected_mutlist,
    manual_fp_file = manual_fp_file,
    tip_label_offset = tip_label_offset,
    tip_label_size = tip_label_size,
    tip_point_size = tip_point_size,
    heatmap_width = heatmap_width,
    heatmap_circos_offset = heatmap_circos_offset,
    flipping_point_size = flipping_point_size,
    plot_height = plot_height,
    plot_width = plot_width,
    verbose = verbose
  )
  
  if (verbose) message("✓ Circular tree plot completed")
  
  # ============================================================
  # Step 3: Heatmap
  # ============================================================
  if (verbose) message("\n[Step 3/4] Heatmap...")
  
  results$heatmap <- plot_heatmap(
    inputpath = inputpath,
    outputpath = outputpath,
    ordered_metadata_file = file.path(outputpath, "ordered_metadata_for_heatmap.txt"),
    colors = heatmap_colors,
    na_color = heatmap_na_color,
    verbose = verbose
  )
  
  if (verbose) message("✓ Heatmap completed")
  
  # ============================================================
  # Step 4: Adjusted Circular Tree Plot (optional)
  # ============================================================
  if (run_adjusted) {
    if (verbose) message("\n[Step 4/4] Adjusted Circular Tree Plot...")
    
    results$circos_adjusted <- plot_circos(
      inputpath = inputpath,
      outputpath = outputpath,
      annotation_file = annotation_file,
      target_mut = target_mut,
      selected_mutlist = selected_mutlist,
      manual_fp_file = manual_fp_file,
      tip_label_offset = adjusted_tip_label_offset,
      tip_label_size = tip_label_size,
      tip_point_size = tip_point_size,
      heatmap_width = heatmap_width,
      heatmap_circos_offset = heatmap_circos_offset,
      flipping_point_size = adjusted_flipping_point_size,
      plot_height = plot_height,
      plot_width = plot_width,
      verbose = verbose
    )
    
    if (verbose) message("✓ Adjusted circular tree plot completed")
  } else {
    if (verbose) message("\n[Step 4/4] Skipped (run_adjusted = FALSE)")
    results$circos_adjusted <- NULL
  }
  
  # ============================================================
  # Summary
  # ============================================================
  end_time <- Sys.time()
  elapsed <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
  
  if (verbose) {
    message("\n", paste(rep("=", 70), collapse = ""))
    message("COMPLETED successfully!")
    message("Time elapsed: ", elapsed, " seconds")
    message("Output directory: ", outputpath)
    
    output_files <- list.files(outputpath)
    if (length(output_files) > 0) {
      message("\nGenerated files:")
      for (f in output_files) {
        message("  - ", f)
      }
    }
    message(paste(rep("=", 70), collapse = ""), "\n")
  }
  
  invisible(results)
}