#' Order mutation matrix for visualization
#'
#' This function performs the complex matrix ordering of cells and mutations
#' based on clonal relationships. It saves the sorted matrix and generates
#' the metadata needed for circos and heatmap visualizations.
#'
#' @param inputpath Path to the PhyloSOLID phylo directory. If data is provided, this can be NULL.
#' @param outputpath Path to output directory for saving sorted files
#' @param data A PhyloData object (optional). If provided, uses this instead of reading files.
#' @param verbose Print progress messages. Default: TRUE
#'
#' @return Invisibly returns a list containing:
#'   - cf_table: The sorted CF table
#'   - input_table: The sorted input table
#'   - cell_order: Order of cells
#'   - mutation_order: Order of mutations
#'   - data: PhyloData object if created
#' @export
#'
#' @examples
#' \dontrun{
#' # Traditional usage
#' result <- order_matrix(
#'   inputpath = "path/to/phylo/",
#'   outputpath = "path/to/output/"
#' )
#'
#' # Using prepared data
#' data <- prepare_data("path/to/phylo/")
#' result <- order_matrix(
#'   data = data,
#'   outputpath = "path/to/output/"
#' )
#' }
order_matrix <- function(
    inputpath = NULL,
    outputpath,
    data = NULL,
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
  
  # ==================== Get data ====================
  if (!is.null(data) && is_PhyloData(data)) {
    # Use provided PhyloData object
    if (verbose) message("Using provided PhyloData object")
    cf_table <- data$sorted_cf
    input_table <- data$input_table
    tree <- data$tree
    tree_data <- data$tree_data
    zero_barcodes <- data$zero_barcodes
    mutation_order <- colnames(cf_table)
    tip_label_order <- data$tree$tip.label
  } else if (!is.null(inputpath)) {
    # Read from files (original behavior)
    if (verbose) message("Reading from files: ", inputpath)
    
    cffile <- file.path(inputpath, "final_cleaned_M_full_basedPivots.filtered_sites_inferred.CFMatrix")
    inputfile <- file.path(inputpath, "final_cleaned_I_full_withNA3_for_circosPlot.txt")
    
    required_files <- c(cffile, inputfile)
    missing_files <- required_files[!file.exists(required_files)]
    if (length(missing_files) > 0) {
      stop("Required files not found:\n", paste("  -", missing_files, collapse = "\n"))
    }
    
    # Read and sort
    raw_cf_table <- read.table(cffile, header = TRUE, sep = "\t", row.names = 1)
    cf_table <- sort_mutation_matrix(raw_cf_table)
    
    raw_input_table <- read.table(inputfile, header = TRUE, sep = "\t", row.names = 1)
    input_table <- raw_input_table[match(rownames(cf_table), rownames(raw_input_table)), 
                                    match(colnames(cf_table), colnames(raw_input_table))]
    
    # Get zero barcodes
    zero_barcodes <- find_zero_barcode(cf_table)
    
    # Build tree
    cffile_no0 <- remove_zero_barcode(cffile)
    tree_data <- converTree::cf2treedata(cffile_no0)
    tree <- dat2tree(tree_data)
    tip_label_order <- tree$tip.label
    mutation_order <- colnames(cf_table)
  } else {
    stop("Either inputpath or data must be provided")
  }
  
  # ==================== Generate ordered metadata ====================
  if (verbose) message("\nGenerating ordered metadata...")
  
  # Get cell order from tree
  leaf_df <- tree_data[tree_data$label %in% tree$tip.label, ]
  leaf_df <- leaf_df %>% left_join(
    data.frame(label = tree$tip.label, order = 1:length(tree$tip.label)),
    by = "label"
  )
  leaf_df <- leaf_df %>% arrange(order)
  tip_label_order <- leaf_df$label
  
  # Generate ordered metadata for heatmap
  metadata_for_heatmap_temp <- cf_table
  metadata_for_heatmap_temp <- metadata_for_heatmap_temp[tip_label_order, mutation_order]
  metadata_for_heatmap_temp$cellIDxmutID <- rownames(metadata_for_heatmap_temp)
  
  metadata_for_heatmap_temp <- metadata_for_heatmap_temp[, c("cellIDxmutID", mutation_order)]
  metadata_for_heatmap_temp[] <- lapply(metadata_for_heatmap_temp, as.character)
  
  mutid_row <- c("mutid_in_heatmap", paste0("M", 1:length(mutation_order)))
  ordered_metadata_for_heatmap <- rbind(mutid_row, metadata_for_heatmap_temp)
  
  # Save files
  ordered_metadata_file <- file.path(outputpath, "ordered_metadata_for_heatmap.txt")
  write.table(ordered_metadata_for_heatmap, file = ordered_metadata_file, 
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  if (verbose) message("  - Saved ordered metadata: ", ordered_metadata_file)
  
  sorted_cf_file <- file.path(outputpath, "sorted_cf_matrix.txt")
  write.table(cf_table, file = sorted_cf_file, sep = "\t", quote = FALSE)
  if (verbose) message("  - Saved sorted CF matrix: ", sorted_cf_file)
  
  # ==================== Summary ====================
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