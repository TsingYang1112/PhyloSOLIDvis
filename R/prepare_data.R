#' Prepare phylogenetic data for visualization
#'
#' Reads and processes all input data once, creating a PhyloData object
#' that can be reused across visualization functions.
#'
#' @param inputpath Path to PhyloSOLID phylo directory
#' @param annotation_file Path to annotation file (optional)
#' @param verbose Print progress messages
#'
#' @return A PhyloData object containing all processed data
#' @export
#'
#' @examples
#' \dontrun{
#' data <- prepare_data(
#'   inputpath = "path/to/phylo/",
#'   annotation_file = "path/to/annotations.txt"
#' )
#' }
prepare_data <- function(
    inputpath,
    annotation_file = NULL,
    verbose = TRUE
) {
  if (verbose) {
    message("\n", paste(rep("=", 60), collapse = ""))
    message("Preparing Phylogenetic Data")
    message(paste(rep("=", 60), collapse = ""))
    message("Input path: ", inputpath)
    if (!is.null(annotation_file)) {
      message("Annotation file: ", annotation_file)
    }
  }
  
  # Set file paths
  cffile <- file.path(inputpath, "final_cleaned_M_full_basedPivots.filtered_sites_inferred.CFMatrix")
  inputfile <- file.path(inputpath, "final_cleaned_I_full_withNA3_for_circosPlot.txt")
  features_file <- file.path(inputpath, "df_flipping_count_for_each_mut.txt")
  total_flipping_file <- file.path(inputpath, "df_total_flipping_count.txt")
  clone_file <- file.path(inputpath, "df_barcode_clones_from_phylo_tree.csv")
  
  # Validate required files
  required_files <- c(cffile, inputfile, features_file, total_flipping_file)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    stop("Required files not found:\n", paste("  -", missing_files, collapse = "\n"))
  }
  
  # ==================== 1. Read CF matrix ====================
  if (verbose) message("\n[1/6] Reading CF matrix...")
  raw_cf <- read.table(cffile, header = TRUE, sep = "\t", row.names = 1)
  if (verbose) message("  - Raw CF matrix: ", nrow(raw_cf), " cells, ", ncol(raw_cf), " mutations")
  
  # ==================== 2. Read and align input table ====================
  if (verbose) message("\n[2/6] Reading input table...")
  raw_input <- read.table(inputfile, header = TRUE, sep = "\t", row.names = 1)
  if (verbose) message("  - Input table: ", nrow(raw_input), " cells, ", ncol(raw_input), " mutations")
  
  # ==================== 3. Sort mutation matrix ====================
  if (verbose) message("\n[3/6] Sorting mutation matrix...")
  sorted_cf <- sort_mutation_matrix(raw_cf)
  input_table <- raw_input[match(rownames(sorted_cf), rownames(raw_input)),
                           match(colnames(sorted_cf), colnames(raw_input))]
  if (verbose) message("  - Sorted: ", nrow(sorted_cf), " cells, ", ncol(sorted_cf), " mutations")
  
  # ==================== 4. Find zero barcodes and build tree ====================
  if (verbose) message("\n[4/6] Building phylogenetic tree...")
  zero_barcodes <- find_zero_barcode(sorted_cf)
  
  # Process cffile (remove zero barcodes) - generate in a temp location
  # Since prepare_data doesn't have an outputpath, use a temporary directory
  temp_dir <- tempdir()
  cffile_no0 <- remove_zero_barcode(cffile, output_dir = temp_dir)
  tree_data <- converTree::cf2treedata(cffile_no0)
  tree <- dat2tree(tree_data)
  
  # Read no-zero matrix
  cf_no0 <- read.table(cffile_no0, header = TRUE, sep = "\t")
  if (verbose) message("  - Tree: ", ape::Ntip(tree), " tips, ", ape::Nnode(tree), " nodes")
  
  # Clean up temp file (optional)
  unlink(cffile_no0)
  
  # ==================== 5. Process annotations ====================
  if (verbose) message("\n[5/6] Processing annotations...")
  annotations <- NULL
  if (!is.null(annotation_file) && file.exists(annotation_file)) {
    df_all_info <- read.table(annotation_file, header = TRUE, sep = "\t")
    annotations <- process_annotations_with_zero(df_all_info, zero_barcodes)
    if (verbose) message("  - Loaded annotations: ", paste(names(annotations), collapse = ", "))
  } else {
    if (verbose) message("  - No annotation file provided")
  }
  
  # ==================== 6. Read flipping data ====================
  if (verbose) message("\n[6/6] Reading flipping data...")
  flipping_data <- read.table(features_file, header = TRUE, sep = "\t")
  total_flipping <- read.table(total_flipping_file, header = TRUE, sep = "\t")
  colnames(total_flipping) <- c("total_flipping_1_to_0", "total_flipping_0_to_1", 
                                 "total_flipping_NA_to_0", "total_flipping_NA_to_1")
  
  # Read clone file if exists (for heatmap)
  clone_data <- NULL
  if (file.exists(clone_file)) {
    clone_data <- read.csv(clone_file, header = TRUE, stringsAsFactors = FALSE)
    if (verbose) message("  - Loaded clone data: ", nrow(clone_data), " rows")
  }
  
  # ==================== Create PhyloData object ====================
  data <- PhyloData(
    raw_cf = raw_cf,
    input_table = input_table,
    sorted_cf = sorted_cf,
    tree = tree,
    tree_data = tree_data,
    tree_no0 = cf_no0,
    zero_barcodes = zero_barcodes,
    annotations = annotations,
    clone_data = clone_data,
    flipping_data = flipping_data,
    total_flipping = total_flipping
  )
  
  if (verbose) {
    message("\n", paste(rep("=", 60), collapse = ""))
    message("Data preparation complete!")
    message(paste(rep("=", 60), collapse = ""))
    print(data)
  }
  
  invisible(data)
}