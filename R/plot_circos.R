#' Create circular phylogenetic tree visualization
#'
#' Generates a publication-ready circular (circos-style) plot from PhyloSOLID output,
#' featuring mutation heatmaps, cell annotations, and genotype flipping markers.
#'
#' @param inputpath Path to directory containing PhyloSOLID output files
#' @param outputpath Path to output directory for saving plots
#' @param annotation_file Path to cell annotation file (TSV format)
#' @param target_mut Target mutation ID to highlight (default: "no")
#' @param conflict_muts Comma-separated conflict mutations (default: NULL)
#' @param selected_mutlist Mutations to include or "all" (default: "all")
#' @param adding_mut_file Additional mutation file path (default: "no")
#' @param manual_fp_file Manual false positive file path (default: "no")
#' @param barcode_file Barcode file for highlighting (default: "no")
#' @param barcode_name Barcode identifier (default: "no")
#' @param tip_label_offset Offset distance for tip labels (default: 10)
#' @param tip_label_size Font size for tip labels (default: 2.5)
#' @param tip_point_size Size of tip points (default: 0.5)
#' @param heatmap_width Width of heatmap track (default: 0.3)
#' @param heatmap_circos_offset Offset for heatmap (default: 0.04)
#' @param flipping_point_size Size of flipping marker points (default: 0.2)
#' @param plot_height Height of output plot in inches (default: 12)
#' @param plot_width Width of output plot in inches (default: 18)
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return Invisibly returns a list containing the plot, legend, data, parameters, and output files
#'
#' @examples
#' \dontrun{
#' result <- plot_circos(
#'   inputpath = "results/phylosolid/",
#'   outputpath = "figures/",
#'   annotation_file = "data/annotations.txt"
#' )
#' }
#'
#' @export
#' @author Qing Yang <qing.yang@nyu.edu>, Yonghe Xia <yonghe.xia@nyu.edu>
plot_circos <- function(
    inputpath,
    outputpath,
    annotation_file,
    target_mut = "no",
    conflict_muts = NULL,
    selected_mutlist = "all",
    adding_mut_file = "no",
    manual_fp_file = "no",
    barcode_file = "no",
    barcode_name = "no",
    tip_label_offset = 10,
    tip_label_size = 2.5,
    tip_point_size = 0.5,
    heatmap_width = 0.3,
    heatmap_circos_offset = 0.04,
    flipping_point_size = 0.2,
    plot_height = 12,
    plot_width = 18,
    verbose = TRUE
) {
  
  # Start timing
  start_time <- Sys.time()
  
  if (verbose) {
    message("\n", paste(rep("=", 60), collapse = ""))
    message("PhyloSOLIDvis - Circular Tree Visualization")
    message(paste(rep("=", 60), collapse = ""))
    message("Start time: ", start_time)
    message("Package version: ", packageVersion("PhyloSOLIDvis"))
  }
  
  # Create output directory
  if (!dir.exists(outputpath)) {
    dir.create(outputpath, recursive = TRUE)
    if (verbose) message("Created output directory: ", outputpath)
  }
  
  # Set file paths
  inputfile <- file.path(inputpath, "final_cleaned_I_full_withNA3_for_circosPlot.txt")
  cffile <- file.path(inputpath, "final_cleaned_M_full_basedPivots.filtered_sites_inferred.CFMatrix")
  features_file <- file.path(inputpath, "df_flipping_count_for_each_mut.txt")
  total_flipping_file <- file.path(inputpath, "df_total_flipping_count.txt")
  
  # Validate input files
  required_files <- c(inputfile, cffile, features_file, total_flipping_file, annotation_file)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    stop("Required files not found:\n", paste("  -", missing_files, collapse = "\n"))
  }
  
  # Read data
  if (verbose) message("\n[1/6] Reading input files...")
  
  raw_cf_table <- read.table(cffile, header = TRUE, sep = "\t", row.names = 1)
  if (verbose) message("  - CF matrix: ", nrow(raw_cf_table), " cells, ", ncol(raw_cf_table), " mutations")
  
  # Sort mutation matrix
  if (verbose) message("\n[2/6] Sorting mutation matrix...")
  cf_table <- sort_mutation_matrix(raw_cf_table)
  
  # Read input table
  raw_input_table <- read.table(inputfile, header = TRUE, sep = "\t", row.names = 1)
  input_table <- raw_input_table[match(rownames(cf_table), rownames(raw_input_table)), 
                                  match(colnames(cf_table), colnames(raw_input_table))]
  
  # Build tree
  if (verbose) message("\n[3/6] Building phylogenetic tree...")
  cffile_no0 <- .remove_zero_barcode(cffile)
  cf_tableno0 <- read.table(cffile_no0, header = TRUE, sep = "\t")
  
  tree_dat <- cf2treedata(cffile_no0)
  tree <- dat2tree(tree_dat)
  if (verbose) message("  - Tree: ", Ntip(tree), " tips, ", Nnode(tree), " nodes")
  
  # Process annotations
  if (verbose) message("\n[4/6] Processing annotations...")
  df_all_info <- read.table(annotation_file, header = TRUE, sep = "\t")
  annotation_results <- .process_annotations(df_all_info)
  
  # Sort cells and identify subclones
  if (verbose) message("\n[5/6] Identifying subclones...")
  cf_tableno0_for_sorting_mat <- cf_tableno0[, 2:ncol(cf_tableno0), drop = FALSE]
  rownames(cf_tableno0_for_sorting_mat) <- cf_tableno0[, 1]
  mat <- as.matrix(cf_tableno0_for_sorting_mat)
  
  cell_clone_structure <- .recursive_clone_partition(rownames(mat), mat, by_row = TRUE)
  sorted_cells <- .extract_ordering(cell_clone_structure)
  sorted_mat <- mat[sorted_cells, , drop = FALSE]
  
  cf_tableno0_sorted <- data.frame(cellIDxmutID = rownames(sorted_mat), 
                                    sorted_mat, 
                                    stringsAsFactors = FALSE,
                                    check.names = FALSE)
  subclones <- subclone_finder(cf_tableno0_sorted)
  if (verbose) message("  - Found ", length(unique(subclones$clone)), " clones")
  
  # Prepare leaf data
  leaf_df <- tree_dat[tree_dat$label %in% tree$tip.label, ]
  subclones_transformed <- tidyr::separate_rows(subclones, cells, sep = ",")
  leaf_df <- leaf_df %>% 
    dplyr::left_join(subclones_transformed, by = c("label" = "cells"))
  leaf_df$clone[is.na(leaf_df$clone)] <- "root_clone"
  
  if (target_mut != "no") {
    leaf_df <- .add_target_clone(leaf_df, cf_tableno0_sorted, target_mut)
    if (verbose) message("  - Targeting mutation: ", target_mut)
  }
  
  leaf_df <- .sort_leaf_df(leaf_df, cf_tableno0_sorted)
  clone_order_tree <- ape::rotateConstr(tree, leaf_df$label)
  
  # Create plot
  if (verbose) message("\n[6/6] Generating circular plot...")
  plot_result <- .create_circos_plot(
    tree = clone_order_tree,
    leaf_df = leaf_df,
    cf_tableno0 = cf_tableno0_sorted,
    input_table = input_table,
    cf_table = cf_table,
    df_all_info = df_all_info,
    annotation_results = annotation_results,
    target_mut = target_mut,
    selected_mutlist = selected_mutlist,
    manual_fp_file = manual_fp_file,
    tip_label_offset = tip_label_offset,
    tip_label_size = tip_label_size,
    tip_point_size = tip_point_size,
    heatmap_width = heatmap_width,
    heatmap_circos_offset = heatmap_circos_offset,
    flipping_point_size = flipping_point_size,
    outputpath = outputpath
  )
  
  # Save outputs
  .save_outputs(plot_result, outputpath, target_mut, plot_width, plot_height)
  
  # Summary
  end_time <- Sys.time()
  elapsed <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
  
  if (verbose) {
    message("\n", paste(rep("=", 60), collapse = ""))
    message("COMPLETED successfully!")
    message("Time elapsed: ", elapsed, " seconds")
    message("Output directory: ", outputpath)
    message(paste(rep("=", 60), collapse = ""), "\n")
  }
  
  invisible(plot_result)
}