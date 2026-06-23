#' Plot heatmap from PhyloSOLID output
#'
#' This function generates a heatmap visualization from PhyloSOLID output files.
#' It uses the ordered metadata, CF matrix, and clone information to create
#' a publication-ready heatmap with proper cell ordering.
#'
#' @param inputpath Path to the PhyloSOLID phylo directory
#' @param outputpath Path to output directory for saving plots
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
#'
#' @examples
#' \dontrun{
#' plot_heatmap(
#'   inputpath = "results/Org4S15D63/03_tree_building/mutation_integrator/phylo/",
#'   outputpath = "figures/heatmap/"
#' )
#' }
plot_heatmap <- function(
    inputpath,
    outputpath,
    ordered_metadata_file = NULL,
    colors = c("#D4E8F0", "#7D2224"),
    na_color = "white",
    show_colnames = TRUE,
    show_rownames = FALSE,
    cluster_cols = FALSE,
    border_color = NA,
    verbose = TRUE
) {
  
  # Start timing
  start_time <- Sys.time()
  
  if (verbose) {
    message("\n", paste(rep("=", 60), collapse = ""))
    message("PhyloSOLIDvis - Heatmap Visualization")
    message(paste(rep("=", 60), collapse = ""))
    message("Start time: ", start_time)
  }
  
  # Create output directory
  if (!dir.exists(outputpath)) {
    dir.create(outputpath, recursive = TRUE)
    if (verbose) message("Created output directory: ", outputpath)
  }
  
  # Set file paths
  cfmatrix_file <- file.path(inputpath, "final_cleaned_M_full_basedPivots.filtered_sites_inferred.CFMatrix")
  clone_file <- file.path(inputpath, "df_barcode_clones_from_phylo_tree.csv")
  input_file <- file.path(inputpath, "final_cleaned_I_full_withNA3_for_circosPlot.txt")
  
  # Validate input files
  required_files <- c(cfmatrix_file, clone_file, input_file)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    stop("Required files not found:\n", paste("  -", missing_files, collapse = "\n"))
  }
  
  # Find ordered_metadata_file if not provided
  if (is.null(ordered_metadata_file)) {
    ordered_metadata_file <- file.path(outputpath, "ordered_metadata_for_heatmap.txt")
  }
  
  # ==================== 1. 从 cfmatrix 转换成 phylo 树 ====================
  if (verbose) message("\n[1/5] Building phylogenetic tree from CF matrix...")
  
  treedata <- converTree::cf2treedata(cfmatrix_file)
  nt2 <- treedata
  phylo_tree <- treedata %>% ape::as.phylo()
  
  # 更新 tip label
  index_label <- function(node) {
    nt2[nt2$node == node, "label"]
  }
  phylo_tree[["tip.label"]] <- Map(index_label, phylo_tree[["tip.label"]]) %>% unlist()
  phylo_tree$node.label <- treedata$label[match(phylo_tree$node.label, treedata$node)]
  
  # 确保有 branch length
  if (is.null(phylo_tree$edge.length) || length(phylo_tree$edge.length) == 0) {
    phylo_tree$edge.length <- rep(1, nrow(phylo_tree$edge))
  }
  
  tree_rds <- file.path(outputpath, "phylo_tree.rds")
  saveRDS(phylo_tree, tree_rds)
  if (verbose) message("  - Tree saved to: ", tree_rds)
  if (verbose) message("  - Tree has ", ape::Ntip(phylo_tree), " tips and ", ape::Nnode(phylo_tree), " nodes")
  
  # ==================== 2. 读取颜色数据 ====================
  if (verbose) message("\n[2/5] Reading clone color data...")
  
  our_color <- read.csv(clone_file, header = TRUE, stringsAsFactors = FALSE)
  our_color$our_color <- our_color$color
  our_color$color <- NULL
  if (verbose) message("  - Loaded ", nrow(our_color), " cells with clone colors")
  
  # ==================== 3. 准备数据并获取 tip 顺序 ====================
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
  
  # ==================== 4. 准备热图数据 ====================
  if (verbose) message("\n[4/5] Preparing heatmap data...")
  
  if (file.exists(ordered_metadata_file)) {
    matrix_circos <- read.csv(ordered_metadata_file, sep = "\t", row.names = 1)
    if (verbose) message("  - Loaded ordered metadata (", nrow(matrix_circos), " rows, ", ncol(matrix_circos), " columns)")
    matrix_circos_no_mutid <- matrix_circos[rownames(matrix_circos) != "mutid_in_heatmap", , drop = FALSE]
  } else {
    if (verbose) message("  - ordered_metadata_file not found: ", ordered_metadata_file)
    if (verbose) message("  - Proceeding without ordered metadata")
    matrix_circos_no_mutid <- NULL
  }
  
  I_for_our_tree <- read.csv(input_file, sep = "\t", row.names = 1)
  if (verbose) message("  - Loaded input matrix (", nrow(I_for_our_tree), " cells, ", ncol(I_for_our_tree), " mutations)")
  
  if (!is.null(matrix_circos_no_mutid)) {
    I_for_our_tree_ordered <- I_for_our_tree[tip_order, colnames(matrix_circos_no_mutid)]
  } else {
    I_for_our_tree_ordered <- I_for_our_tree[tip_order, ]
  }
  
  plot_matrix_our <- as.matrix(I_for_our_tree_ordered)
  plot_matrix_our_numeric <- apply(plot_matrix_our, 2, as.numeric)
  rownames(plot_matrix_our_numeric) <- rownames(plot_matrix_our)
  plot_matrix_our_numeric[plot_matrix_our_numeric == 3] <- NA
  
  if (verbose) message("  - Heatmap matrix: ", nrow(plot_matrix_our_numeric), " cells, ", ncol(plot_matrix_our_numeric), " mutations")
  
  # ==================== 5. 绘制热图 ====================
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
  
  # png_file <- file.path(outputpath, "heatmap_preview.png")
  # png(png_file, width = 800, height = 600)
  
  # pheatmap::pheatmap(
  #   plot_matrix_our_numeric,
  #   color = colors,
  #   na_col = na_color,
  #   cluster_rows = FALSE,
  #   cluster_cols = cluster_cols,
  #   show_rownames = FALSE,
  #   show_colnames = show_colnames,
  #   legend = FALSE,
  #   border_color = border_color
  # )
  
  # dev.off()
  # if (verbose) message("  - Heatmap preview saved to: ", png_file)
  
  end_time <- Sys.time()
  if (verbose) {
    message("\n", paste(rep("=", 60), collapse = ""))
    message("COMPLETED successfully!")
    message("Time elapsed: ", round(as.numeric(difftime(end_time, start_time, units = "secs")), 2), " seconds")
    message("Output directory: ", outputpath)
    message("  - Tree: ", tree_rds)
    message("  - Heatmap PDF: ", pdf_file)
    # message("  - Heatmap PNG: ", png_file)
    message(paste(rep("=", 60), collapse = ""), "\n")
  }
  
  invisible(list(
    tree = phylo_tree,
    matrix = plot_matrix_our_numeric,
    tip_order = tip_order,
    # output_files = c(tree_rds, pdf_file, png_file)
    output_files = c(tree_rds, pdf_file)
  ))
}