#' Order mutation matrix for visualization
#'
#' This function performs the complex matrix ordering of cells and mutations
#' based on clonal relationships. It saves the sorted matrix and generates
#' the metadata needed for circos and heatmap visualizations.
#'
#' @param inputpath Path to the PhyloSOLID phylo directory containing:
#'   - final_cleaned_M_full_basedPivots.filtered_sites_inferred.CFMatrix
#'   - final_cleaned_I_full_withNA3_for_circosPlot.txt
#' @param outputpath Path to output directory for saving sorted files
#' @param verbose Print progress messages. Default: TRUE
#'
#' @return Invisibly returns a list containing:
#'   - cf_table: The sorted CF table
#'   - input_table: The sorted input table
#'   - cell_order: Order of cells
#'   - mutation_order: Order of mutations
#' @export
#'
#' @examples
#' \dontrun{
#' result <- order_matrix(
#'   inputpath = "path/to/phylo/",
#'   outputpath = "path/to/output/"
#' )
#' }
order_matrix <- function(
    inputpath,
    outputpath,
    verbose = TRUE
) {
  
  # Start timing
  start_time <- Sys.time()
  
  if (verbose) {
    message("\n", paste(rep("=", 60), collapse = ""))
    message("PhyloSOLIDvis - Matrix Ordering")
    message(paste(rep("=", 60), collapse = ""))
    message("Start time: ", start_time)
  }
  
  # Create output directory
  if (!dir.exists(outputpath)) {
    dir.create(outputpath, recursive = TRUE)
    if (verbose) message("Created output directory: ", outputpath)
  }
  
  # Set file paths
  cffile <- file.path(inputpath, "final_cleaned_M_full_basedPivots.filtered_sites_inferred.CFMatrix")
  inputfile <- file.path(inputpath, "final_cleaned_I_full_withNA3_for_circosPlot.txt")
  
  # Validate input files
  required_files <- c(cffile, inputfile)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    stop("Required files not found:\n", paste("  -", missing_files, collapse = "\n"))
  }
  
  # 1. Read and sort CF matrix
  if (verbose) message("\n[1/4] Reading and sorting CF matrix...")
  
  raw_cf_table <- read.table(cffile, header = TRUE, sep = "\t", row.names = 1)
  if (verbose) message("  - CF matrix: ", nrow(raw_cf_table), " cells, ", ncol(raw_cf_table), " mutations")
  
  # Sort mutation matrix (using existing function)
  cf_table <- sort_mutation_matrix(raw_cf_table)
  if (verbose) message("  - Matrix sorted")
  
  # 2. Read and align input table
  if (verbose) message("\n[2/4] Reading and aligning input table...")
  
  raw_input_table <- read.table(inputfile, header = TRUE, sep = "\t", row.names = 1)
  input_table <- raw_input_table[match(rownames(cf_table), rownames(raw_input_table)), 
                                  match(colnames(cf_table), colnames(raw_input_table))]
  if (verbose) message("  - Input table aligned")
  
  # 3. Remove zero barcodes and get order
  if (verbose) message("\n[3/4] Removing zero barcodes and getting order...")
  
  outputree_dat_zero_barcode <- find_zero_barcode(cf_table)
  
  # Process cffile
  cffile_no0 <- remove_zero_barcode(cffile)
  cf_tableno0 <- read.table(cffile_no0, header = TRUE, sep = "\t")
  if (verbose) message("  - Spots drawn: ", dim(cf_tableno0)[1])
  
  # Build tree to get cell order
  tree_dat <- converTree::cf2treedata(cffile_no0)
  tree <- dat2tree(tree_dat)
  
  # Get leaf order from tree
  leaf_df <- tree_dat[tree_dat$label %in% tree$tip.label, ]
  leaf_df <- leaf_df %>% left_join(
    data.frame(label = tree$tip.label, order = 1:length(tree$tip.label)),
    by = "label"
  )
  leaf_df <- leaf_df %>% arrange(order)
  tip_label_order <- leaf_df$label
  
  # 4. Generate and save ordered metadata
  if (verbose) message("\n[4/4] Generating ordered metadata...")
  
  # Get mutation order (from cf_table column names)
  mutation_order <- colnames(cf_table)
  
  # Generate ordered metadata for heatmap
  metadata_for_heatmap_temp <- cf_table
  metadata_for_heatmap_temp <- metadata_for_heatmap_temp[tip_label_order, mutation_order]
  metadata_for_heatmap_temp$cellIDxmutID <- rownames(metadata_for_heatmap_temp)
  
  # Add mutid_in_heatmap row
  metadata_for_heatmap_temp <- metadata_for_heatmap_temp[, c("cellIDxmutID", mutation_order)]
  metadata_for_heatmap_temp[] <- lapply(metadata_for_heatmap_temp, as.character)
  
  mutid_row <- c("mutid_in_heatmap", paste0("M", 1:length(mutation_order)))
  ordered_metadata_for_heatmap <- rbind(mutid_row, metadata_for_heatmap_temp)
  
  # Save ordered metadata
  ordered_metadata_file <- file.path(outputpath, "ordered_metadata_for_heatmap.txt")
  write.table(ordered_metadata_for_heatmap, file = ordered_metadata_file, 
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  if (verbose) message("  - Saved ordered metadata: ", ordered_metadata_file)
  
  # Save sorted CF matrix
  sorted_cf_file <- file.path(outputpath, "sorted_cf_matrix.txt")
  write.table(cf_table, file = sorted_cf_file, sep = "\t", quote = FALSE)
  if (verbose) message("  - Saved sorted CF matrix: ", sorted_cf_file)
  
  # Summary
  end_time <- Sys.time()
  if (verbose) {
    message("\n", paste(rep("=", 60), collapse = ""))
    message("COMPLETED successfully!")
    message("Time elapsed: ", round(as.numeric(difftime(end_time, start_time, units = "secs")), 2), " seconds")
    message("Output directory: ", outputpath)
    message("  - sorted_cf_matrix.txt")
    message("  - ordered_metadata_for_heatmap.txt")
    message(paste(rep("=", 60), collapse = ""), "\n")
  }
  
  invisible(list(
    cf_table = cf_table,
    input_table = input_table,
    cell_order = tip_label_order,
    mutation_order = mutation_order,
    sorted_files = c(sorted_cf_file, ordered_metadata_file)
  ))
}