#' Plot heatmap from PhyloSOLID output
#'
#' This function generates a heatmap visualization from PhyloSOLID output files.
#'
#' @param inputpath Path to the PhyloSOLID phylo directory (optional if data provided)
#' @param outputpath Path to output directory for saving plots
#' @param data A PhyloData object (optional). If provided, uses this instead of reading files.
#' @param tree A phylo object (optional). If provided, uses this tree instead of building one.
#' @param ordered_metadata_file Path to ordered_metadata_for_heatmap.txt file.
#'   If NULL, will look for it in the outputpath.
#' @param colors Character vector of colors for the heatmap (0 and 1 values).
#'   Default: c("#D4E8F0", "#7D2224")
#' @param na_color Color for NA/missing values. Default: "white"
#' @param show_colnames Whether to show column names (mutation IDs). Default: TRUE
#' @param show_rownames Whether to show row names (cell IDs). Default: FALSE
#' @param cluster_cols Whether to cluster columns (mutations). Default: FALSE
#' @param border_color Border color for heatmap cells. Default: NA
#' @param verbose Print progress messages. Default: TRUE
#'
#' @return Invisibly returns the heatmap matrix and the tree object
#' @export
plot_heatmap <- function(
    inputpath = NULL,
    outputpath,
    data = NULL,
    tree = NULL,
    ordered_metadata_file = NULL,
    colors = c("#D4E8F0", "#7D2224"),
    na_color = "white",
    show_colnames = TRUE,
    show_rownames = FALSE,
    cluster_cols = FALSE,
    border_color = NA,
    verbose = TRUE
) {
  
  start_time <- Sys.time()
  
  if (verbose) {
    message("\n", paste(rep("=", 60), collapse = ""))
    message("PhyloSOLIDvis - Heatmap Visualization")
    message(paste(rep("=", 60), collapse = ""))
    message("Start time: ", start_time)
  }
  
  if (!dir.exists(outputpath)) {
    dir.create(outputpath, recursive = TRUE)
    if (verbose) message("Created output directory: ", outputpath)
  }
  
  if (is.null(ordered_metadata_file)) {
    ordered_metadata_file <- file.path(outputpath, "ordered_metadata_for_heatmap.txt")
  }
  
  use_prepared_data <- !is.null(data) && is_PhyloData(data)
  
  if (use_prepared_data) {
    if (verbose) message("Using provided PhyloData object")
    cfmatrix_file <- NULL
    clone_data <- data$clone_data
    input_table <- data$input_table
  } else if (!is.null(inputpath)) {
    if (verbose) message("Reading from files: ", inputpath)
    cfmatrix_file <- file.path(inputpath, "final_cleaned_M_full_basedPivots.filtered_sites_inferred.CFMatrix")
    clone_file <- file.path(inputpath, "df_barcode_clones_from_phylo_tree.csv")
    input_file <- file.path(inputpath, "final_cleaned_I_full_withNA3_for_circosPlot.txt")
    
    required_files <- c(cfmatrix_file, clone_file, input_file)
    missing_files <- required_files[!file.exists(required_files)]
    if (length(missing_files) > 0) {
      stop("Required files not found:\n", paste("  -", missing_files, collapse = "\n"))
    }
    
    clone_data <- read.csv(clone_file, header = TRUE, stringsAsFactors = FALSE)
    input_table <- read.csv(input_file, sep = "\t", row.names = 1)
  } else {
    stop("Either inputpath or data must be provided")
  }
  
  if (verbose) message("\n[1/5] Getting phylogenetic tree...")
  
  if (!is.null(tree) && inherits(tree, "phylo")) {
    phylo_tree <- tree
    if (verbose) message("  - Using provided tree (", ape::Ntip(phylo_tree), " tips, ", ape::Nnode(phylo_tree), " nodes)")
  } else if (use_prepared_data && !is.null(data$tree)) {
    phylo_tree <- data$tree
    if (verbose) message("  - Using tree from PhyloData (", ape::Ntip(phylo_tree), " tips, ", ape::Nnode(phylo_tree), " nodes)")
  } else if (!is.null(cfmatrix_file) && file.exists(cfmatrix_file)) {
    if (verbose) message("  - Building tree from CF matrix...")
    treedata <- converTree::cf2treedata(cfmatrix_file)
    nt2 <- treedata
    phylo_tree <- treedata %>% ape::as.phylo()
    
    index_label <- function(node) {
      nt2[nt2$node == node, "label"]
    }
    phylo_tree[["tip.label"]] <- Map(index_label, phylo_tree[["tip.label"]]) %>% unlist()
    phylo_tree$node.label <- treedata$label[match(phylo_tree$node.label, treedata$node)]
    
    if (is.null(phylo_tree$edge.length) || length(phylo_tree$edge.length) == 0) {
      phylo_tree$edge.length <- rep(1, nrow(phylo_tree$edge))
    }
    
    tree_rds <- file.path(outputpath, "phylo_tree.rds")
    saveRDS(phylo_tree, tree_rds)
    if (verbose) message("  - Tree saved to: ", tree_rds)
    if (verbose) message("  - Tree has ", ape::Ntip(phylo_tree), " tips and ", ape::Nnode(phylo_tree), " nodes")
  } else {
    stop("No tree available. Please provide a tree, PhyloData object, or inputpath.")
  }
  
  if (verbose) message("\n[2/5] Reading clone color data...")
  
  if (!is.null(clone_data)) {
    our_color <- clone_data
    if (!"our_color" %in% colnames(our_color)) {
      our_color$our_color <- our_color$color
      our_color$color <- NULL
    }
    if (verbose) message("  - Loaded ", nrow(our_color), " cells with clone colors")
  } else {
    our_color <- data.frame(label = phylo_tree$tip.label, our_color = "Unknown")
    if (verbose) message("  - No clone color data available, using default")
  }
  
  if (verbose) message("\n[3/5] Getting tip order from tree...")
  
  tip_data <- our_color
  rownames(tip_data) <- tip_data$label
  tip_data$our_color <- as.character(tip_data$our_color)
  tip_data$our_color[tip_data$our_color == "" | is.na(tip_data$our_color)] <- "Unknown"
  tip_data$our_color <- factor(tip_data$our_color)
  
  p <- ggtree(phylo_tree, branch.length = "none") %<+% tip_data +
    geom_tree(color = "gray50", size = 0.5, lineend = "round") +
    theme_tree2() +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.margin = margin(0.5, 0.5, 0.5, 1, "cm")
    ) +
    coord_cartesian(clip = "off")
  
  tip_order <- get_taxa_name(p)
  if (verbose) message("  - Tip order obtained (", length(tip_order), " tips)")
  
  if (verbose) message("\n[4/5] Preparing heatmap data...")
  
  if (file.exists(ordered_metadata_file)) {
    matrix_circos <- read.csv(ordered_metadata_file, sep = "\t", row.names = 1)
    if (verbose) message("  - Loaded ordered metadata (", nrow(matrix_circos), " rows, ", ncol(matrix_circos), " columns)")
    matrix_circos_no_mutid <- matrix_circos[rownames(matrix_circos) != "mutid_in_heatmap", , drop = FALSE]
    cols_to_use <- colnames(matrix_circos_no_mutid)
  } else {
    if (verbose) message("  - ordered_metadata_file not found: ", ordered_metadata_file)
    if (verbose) message("  - Using input matrix directly")
    matrix_circos_no_mutid <- NULL
    cols_to_use <- colnames(input_table)
  }
  
  if (!is.null(input_table)) {
    if (!is.null(matrix_circos_no_mutid)) {
      I_for_our_tree_ordered <- input_table[tip_order, cols_to_use]
    } else {
      I_for_our_tree_ordered <- input_table[tip_order, ]
    }
    
    plot_matrix_our <- as.matrix(I_for_our_tree_ordered)
    plot_matrix_our_numeric <- apply(plot_matrix_our, 2, as.numeric)
    rownames(plot_matrix_our_numeric) <- rownames(plot_matrix_our)
    plot_matrix_our_numeric[plot_matrix_our_numeric == 3] <- NA
    
    if (verbose) message("  - Heatmap matrix: ", nrow(plot_matrix_our_numeric), " cells, ", ncol(plot_matrix_our_numeric), " mutations")
  } else {
    stop("No input data available for heatmap")
  }
  
  if (verbose) message("\n[5/5] Drawing heatmap...")
  
  pdf_file <- file.path(outputpath, "heatmap_and_histograms_for_our_tree.pdf")
  pdf(pdf_file, width = 10, height = 8)
  
  pheatmap::pheatmap(
    plot_matrix_our_numeric,
    color = colors,
    na_col = na_color,
    cluster_rows = FALSE,
    cluster_cols = cluster_cols,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    legend = FALSE,
    border_color = border_color
  )
  
  dev.off()
  if (verbose) message("  - Heatmap saved to: ", pdf_file)
  
  end_time <- Sys.time()
  if (verbose) {
    message("\n", paste(rep("=", 60), collapse = ""))
    message("COMPLETED successfully!")
    message("Time elapsed: ", round(as.numeric(difftime(end_time, start_time, units = "secs")), 2), " seconds")
    message("Output directory: ", outputpath)
    message("  - Heatmap PDF: ", pdf_file)
    message(paste(rep("=", 60), collapse = ""), "\n")
  }
  
  invisible(list(
    tree = phylo_tree,
    matrix = plot_matrix_our_numeric,
    tip_order = tip_order,
    output_files = pdf_file
  ))
}