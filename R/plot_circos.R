#' Plot circular phylogenetic tree
#'
#' @param inputpath Path to input data directory
#' @param outputpath Path to output directory
#' @param annotation_file Path to annotation file
#' @param target_mut Target mutation to highlight (default: "no")
#' @param selected_mutlist Selected mutations (default: "all")
#' @param manual_fp_file Manual false positive file (default: "no")
#' @param tip_label_offset Offset for tip labels (default: 6)
#' @param tip_label_size Size of tip labels (default: 2.5)
#' @param tip_point_size Size of tip points (default: 0.5)
#' @param heatmap_width Width of heatmap (default: 0.3)
#' @param heatmap_circos_offset Offset for heatmap (default: 0.05)
#' @param flipping_point_size Size of flipping points (default: 1.3)
#' @param plot_height Plot height (default: 12)
#' @param plot_width Plot width (default: 18)
#' @param verbose Print progress (default: TRUE)
#'
#' @return Invisible list with plot results
#' @export
plot_circos <- function(
    inputpath,
    outputpath,
    annotation_file,
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
    verbose = TRUE
) {
  
  # Start timing
  start_time <- Sys.time()
  
  if (verbose) {
    message("\n", paste(rep("=", 60), collapse = ""))
    message("PhyloSOLIDvis - Circular Tree Visualization")
    message(paste(rep("=", 60), collapse = ""))
    message("Start time: ", start_time)
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
  all_mutid <- colnames(cf_table)
  
  # Read input table
  raw_input_table <- read.table(inputfile, header = TRUE, sep = "\t", row.names = 1)
  input_table <- raw_input_table[match(rownames(cf_table), rownames(raw_input_table)), 
                                  match(colnames(cf_table), colnames(raw_input_table))]
  
  # Remove all 0 spots
  outputree_dat_zero_barcode <- find_zero_barcode(cf_table)
  
  # Process cffile
  cffile_no0 <- remove_zero_barcode(cffile)
  cf_tableno0 <- read.table(cffile_no0, header = TRUE, sep = "\t")
  if (verbose) message("  - Spots drawn: ", dim(cf_tableno0)[1])
  
  # Build tree
  if (verbose) message("\n[3/6] Building phylogenetic tree...")
  tree_dat <- converTree::cf2treedata(cffile_no0)
  tree <- dat2tree(tree_dat)
  if (verbose) message("  - Tree: ", ape::Ntip(tree), " tips, ", ape::Nnode(tree), " nodes")
  
  # Process annotations
  if (verbose) message("\n[4/6] Processing annotations...")
  df_all_info <- read.table(annotation_file, header = TRUE, sep = "\t")
  
  # Process annotation data
  anno_color <- c()
  
  # Cell cluster
  if ("cluster_info" %in% colnames(df_all_info)) {
    df_cluster <- df_all_info[, c("barcode", "cluster_info")]
    df_cluster$cluster_info <- trimws(as.character(df_cluster$cluster_info))
    df_cluster <- df_cluster[!df_cluster$barcode %in% outputree_dat_zero_barcode, ]
    colnames(df_cluster) <- c("label", "cluster")
    cluster_levels <- sort(unique(df_cluster$cluster))
    df_cluster$cluster <- factor(df_cluster$cluster, levels = cluster_levels)
    set.seed(42)
    df_cluster_color <- grDevices::rgb(runif(length(cluster_levels)), runif(length(cluster_levels)), runif(length(cluster_levels)))
    names(df_cluster_color) <- cluster_levels
    anno_color <- c(anno_color, df_cluster_color)
  }
  
  # Cell type
  if ("cell_type" %in% colnames(df_all_info)) {
    df_celltype <- df_all_info[, c("barcode", "cell_type")]
    df_celltype$cell_type <- trimws(as.character(df_celltype$cell_type))
    df_celltype <- df_celltype[!df_celltype$barcode %in% outputree_dat_zero_barcode, ]
    colnames(df_celltype) <- c("label", "celltype")
    celltype_levels <- sort(unique(df_celltype$celltype))
    df_celltype$celltype <- factor(df_celltype$celltype, levels = celltype_levels)
    df_celltype_color <- paletteer::paletteer_d("ggthemes::stata_s1rcolor")[seq_along(celltype_levels)]
    names(df_celltype_color) <- celltype_levels
    anno_color <- c(anno_color, df_celltype_color)
  }
  
  # Sample
  if ("sample" %in% colnames(df_all_info)) {
    df_sample <- df_all_info[, c("barcode", "sample")]
    df_sample$sample <- trimws(as.character(df_sample$sample))
    df_sample <- df_sample[!df_sample$barcode %in% outputree_dat_zero_barcode, ]
    colnames(df_sample) <- c("label", "sample")
    df_sample_color <- c("#ebb16e", "#b2aad3")
    names(df_sample_color) <- c('P6_ST_vis_rep2', 'P6_ST_vis_rep1')
    anno_color <- c(anno_color, df_sample_color)
  }
  
  # Tumor score
  tumor_breaks <- NULL
  if ("tumor_score" %in% colnames(df_all_info)) {
    df_tumor <- df_all_info[, c("barcode", "tumor_score")]
    df_tumor$tumor_score <- as.numeric(as.character(df_tumor$tumor_score))
    df_tumor <- df_tumor[!df_tumor$barcode %in% outputree_dat_zero_barcode, ]
    min_value <- min(df_tumor$tumor_score[df_tumor$tumor_score != 0])
    max_value <- max(df_tumor$tumor_score)
    colnames(df_tumor) <- c("label", "tumor_score")
    df_tumor_color <- c("#00a3c4", "#ffebf6", "#ff64be")
    tumor_breaks <- c(0, min_value, max_value)
    anno_color <- c(anno_color, df_tumor_color)
  }
  
  # B cell proportion
  bcell_breaks <- NULL
  if ("B_cell_prop" %in% colnames(df_all_info)) {
    df_Bcell <- df_all_info[, c("barcode", "B_cell_prop")]
    df_Bcell$B_cell_prop <- as.numeric(as.character(df_Bcell$B_cell_prop))
    df_Bcell <- df_Bcell[!df_Bcell$barcode %in% outputree_dat_zero_barcode, ]
    min_value <- min(df_Bcell$B_cell_prop[df_Bcell$B_cell_prop != 0])
    max_value <- max(df_Bcell$B_cell_prop)
    median_value <- median(df_Bcell$B_cell_prop)
    colnames(df_Bcell) <- c("label", "B_cell_prop")
    df_Bcell_color <- c("#F4D166", "#DBE8B4", "#24693D")
    bcell_breaks <- c(min_value, median_value, max_value)
    anno_color <- c(anno_color, df_Bcell_color)
  }
  
  # Sort cells and identify subclones
  if (verbose) message("\n[5/6] Identifying subclones...")
  cf_tableno0_for_sorting_mat <- cf_tableno0[, 2:ncol(cf_tableno0), drop = FALSE]
  rownames(cf_tableno0_for_sorting_mat) <- cf_tableno0[, 1]
  mat <- as.matrix(cf_tableno0_for_sorting_mat)
  
  cell_clone_structure <- recursive_clone_partition(rownames(mat), mat, by_row = TRUE)
  sorted_cells <- extract_mutation_order(cell_clone_structure)
  sorted_mut_indices <- extract_mutation_order(recursive_clone_partition(1:ncol(mat), mat, by_row = FALSE))
  sorted_mutations <- colnames(mat)[sorted_mut_indices]
  
  sorted_mat <- mat[sorted_cells, sorted_mutations, drop = FALSE]
  write.table(sorted_mat, file = file.path(outputpath, "sorted_cf_matrix.txt"), sep = "\t", quote = FALSE)
  
  cf_tableno0_sorted <- data.frame(cellIDxmutID = rownames(sorted_mat), sorted_mat, stringsAsFactors = FALSE, check.names = FALSE)
  cf_tableno0_sorted[, -1] <- lapply(cf_tableno0_sorted[, -1], as.numeric)
  rownames(cf_tableno0_sorted) <- rownames(cf_tableno0)
  
  subclones <- subclone_finder(cf_tableno0_sorted)
  subclones$clone <- unlist(subclones$clone)
  
  # Prepare subclone info
  subclones_transformed <- tidyr::separate_rows(subclones, cells, sep = ",")
  
  if (exists("df_cluster")) {
    subclones_transformed <- subclones_transformed %>%
      left_join(df_cluster, by = c("cells" = "label"))
  } else if (exists("df_celltype")) {
    subclones_transformed <- subclones_transformed %>%
      left_join(df_celltype, by = c("cells" = "label"))
  }
  
  # Prepare leaf data
  leaf_df <- tree_dat[tree_dat$label %in% tree$tip.label, ]
  leaf_df <- leaf_df %>% left_join(subclones_transformed, by = c("label" = "cells"))
  leaf_df$clone[is.na(leaf_df$clone)] <- "root_clone"
  
  if (target_mut != "no") {
    target_clone <- cf_tableno0_sorted[cf_tableno0_sorted[[target_mut]] == 1, "cellIDxmutID"]
    leaf_df <- leaf_df %>% 
      mutate(target_clone = ifelse(label %in% target_clone, "clone_under_target_mut", "others"))
  }
  
  # Order cells by clone
  rownames(cf_tableno0_sorted) <- cf_tableno0_sorted$cellIDxmutID
  mutation_mat <- cf_tableno0_sorted[, -1, drop = FALSE]
  leaf_df$mutation_count <- rowSums(mutation_mat[leaf_df$label, , drop = FALSE])
  leaf_df <- leaf_df %>%
    arrange(subclone_size, clone_order) %>%
    group_by(clone_order) %>%
    arrange(desc(mutation_count), .by_group = TRUE) %>%
    ungroup()
  
  clone_order_tree <- ape::rotateConstr(tree, leaf_df$label)
  
  # Get mutation info
  mutant_cell_size <- data.frame(
    mutant = colnames(cf_tableno0)[-1],
    mutant_cell_size = colSums(cf_tableno0[, -1])
  )
  
  # Create merged_df and leaf_df_cell_size
  merged_df <- leaf_df %>%
    left_join(tree_dat, by = c("parent" = "node"))
  
  merged_df2 <- tidyr::separate_rows(merged_df, label.y, sep = "\\|")
  
  leaf_df_cell_size <- merged_df2 %>%
    left_join(mutant_cell_size, by = c("label.y" = "mutant"))
  
  # Prepare flipping data
  input_table_long <- input_table %>% tibble::as_tibble(rownames = "barcode") %>%
    tidyr::pivot_longer(-barcode, names_to = "mutation", values_to = "input")
  
  cf_table_long <- cf_table %>% tibble::as_tibble(rownames = "barcode") %>%
    tidyr::pivot_longer(-barcode, names_to = "mutation", values_to = "cf")
  
  all_long <- merge(input_table_long, cf_table_long)
  all_long$flipped_type <- paste(all_long$input, all_long$cf, sep = ">")
  
  all_wide <- all_long %>%
    select(-input, -cf) %>%
    tidyr::pivot_wider(values_from = flipped_type, names_from = mutation, values_fill = "not_flipped") %>%
    tibble::column_to_rownames(., "barcode")
  
  all_wide <- as.data.frame(all_wide)
  all_mat <- as.matrix(all_wide)
  
  # Color definitions
  color_define <- structure(c("#2c81be", "#e17259", "#d9d9d9"),
                           names = c("absence", "presence", "missing"))
  
  unique_types <- unique(all_long$flipped_type)
  no_flipped_index <- strsplit(unique_types, ">") %>%
    lapply(function(x) var(x) == 0) %>%
    unlist()
  
  colors_before <- dplyr::case_when(
    startsWith(unique_types, "0") ~ color_define["absence"],
    startsWith(unique_types, "1") ~ color_define["presence"],
    startsWith(unique_types, "3") ~ color_define["missing"]
  )
  names(colors_before) <- unique_types
  
  colors_after <- dplyr::case_when(
    endsWith(unique_types, "0") ~ color_define["absence"],
    endsWith(unique_types, "1") ~ color_define["presence"],
    endsWith(unique_types, "3") ~ color_define["missing"]
  )
  names(colors_after) <- unique_types
  
  # Prepare sub data
  reordered_heatmap <- cf_tableno0[match(leaf_df$label, cf_tableno0$cellIDxmutID), colnames(cf_table)]
  
  sub <- all_mat[!rownames(all_mat) %in% outputree_dat_zero_barcode, ]
  colnames_order <- colnames(reordered_heatmap)[colnames(reordered_heatmap) %in% colnames(sub)]
  sub <- sub[clone_order_tree$tip.label, colnames_order]
  sub <- as.data.frame(sub)
  sub$cellIDxmutID <- rownames(sub)
  
  sub_long <- sub %>%
    pivot_longer(-cellIDxmutID, names_to = "mutants", values_to = "mut_count")
  
  # Clone colors
  clone_order <- unique(subclones[order(subclones$clone_order, decreasing = TRUE), ]$clone)
  roots_index <- grep("root", clone_order)
  colfunc <- grDevices::colorRampPalette(c("#7b767a", "#a49b95", "#dad7d8"))
  roots_color <- colfunc(length(roots_index))
  
  clone_index <- !grepl("root", clone_order)
  clone_num <- sum(clone_index)
  clone_color_o <- grDevices::hcl(h = seq(15, 375, length = clone_num + 1), c = 100, l = 65)[1:clone_num]
  clone_color <- vector(length = length(clone_order))
  clone_color[roots_index] <- roots_color
  clone_color[clone_index] <- clone_color_o
  names(clone_color) <- clone_order
  
  # Subclone colors
  subclones_transformed2 <- tidyr::separate_rows(subclones, subclone, sep = ",")
  temp_cell_size <- unique(stats::na.omit(leaf_df_cell_size[, c("subclone_num", "mutant_cell_size", "label.y")]))
  subclones_transformed3 <- subclones_transformed2 %>%
    left_join(temp_cell_size, by = c("subclone_num" = "subclone_num"))
  
  subclones_transformed4 <- subclones_transformed3 %>%
    arrange(desc(mutant_cell_size), clone_order)
  
  clone_subclone_names <- unique(subclones_transformed4[, c("clone", "subclone_num", "label.y", "clone_order")])
  clone_subclone_names <- clone_subclone_names %>% 
    left_join(subclones[, c("subclone_num", "subclone")], by = c("subclone_num" = "subclone_num"))
  clone_subclone_names <- clone_subclone_names[order(clone_subclone_names$clone_order), ]
  clone_subclone_names$color <- clone_color[clone_subclone_names$clone]
  
  subclone_color <- clone_color[clone_subclone_names$clone]
  names(subclone_color) <- clone_subclone_names$subclone
  
  subclone_mut_color <- clone_color[clone_subclone_names$clone]
  names(subclone_mut_color) <- clone_subclone_names$label.y
  
  # Mutation colors
  if (selected_mutlist == "all" | is.na(selected_mutlist)) {
    pivot_muts <- colnames(cf_table)
  } else {
    pivot_muts <- unlist(strsplit(selected_mutlist, ",\\s*"))
  }
  color_muts <- pals::viridis(length(pivot_muts))
  names(color_muts) <- pivot_muts
  
  sub_long$mutants_symbol <- paste("M", as.numeric(factor(sub_long$mutants, levels = pivot_muts)), sep = "")
  for_sort_symbol <- unique(sub_long[, c("mutants", "mutants_symbol")])
  pivot_muts_symbol <- for_sort_symbol$mutants_symbol
  pivot_muts_sub <- for_sort_symbol$mutants
  
  # Handle manual fp flipping
  if (manual_fp_file == "no") {
    new_colors_before <- colors_before
    new_colors_after <- colors_after
  } else {
    df_manual_fp <- read.table(manual_fp_file, sep = '\t', header = TRUE)
    sub_long$mut_count <- ifelse(paste(sub_long$cellIDxmutID, sub_long$mutants) %in% 
                                 paste(df_manual_fp$cellID, df_manual_fp$mutID),
                                 paste0("manual1", sub(".*(>.)", "\\1", sub_long$mut_count)), 
                                 sub_long$mut_count)
    new_colors_before <- c(colors_before, "manual1>0" = "#FFFF00", "manual1>1" = "#FFFF00")
    new_colors_after <- c(colors_after, "manual1>0" = "#2c81be", "manual1>1" = "#e17259")
    unique_mut_counts <- names(new_colors_before)
    no_flipped_index <- sapply(unique_mut_counts, function(x) {
      if ("manual" %in% x) {
        FALSE
      } else {
        x_split <- strsplit(x, ">")[[1]]
        length(x_split) == 2 && x_split[1] == x_split[2]
      }
    })
  }
  
  labels_dict <- setNames(format_flipping_label(names(new_colors_before)), names(new_colors_before))
  
  # Create title theme
  title_theme <- ggplot2::element_text(size = 20, colour = "#3c3c3c", angle = 0)
  
  # Save tree data
  saveRDS(clone_order_tree, file.path(outputpath, "CNVtree_data_clone_order_tree.rds"))
  
  # Plot tree
  num_cells <- length(clone_order_tree$tip.label)
  
  # Trim node labels
  clone_order_tree_trim <- clone_order_tree
  max_mutnum_on_one_node <- max(sapply(clone_order_tree$node.label, function(x) {
    if (!is.na(x)) return(length(strsplit(x, "\\|")[[1]]))
    else return(0)
  }))
  max_mutnum_on_one_node <- pmin(max_mutnum_on_one_node, 5)
  clone_order_tree_trim$node.label <- sapply(clone_order_tree_trim$node.label, 
    function(x) trim_label(x, target_mut = target_mut, default_n = max_mutnum_on_one_node))
  
  # Build tree plot
  if (num_cells >= 72) {
    if (target_mut == "no") {
      No_target_clone_mut_color <- rep("white", length(clone_order_tree_trim$node.label))
      names(No_target_clone_mut_color) <- sapply(clone_order_tree_trim$node.label, function(x) gsub("\\|", "\n", x))
      
      p1_tree <- ggtree(clone_order_tree_trim, size = 0.1, branch.length = "none", ladderize = FALSE, layout = "radial", color = "black") +
        geom_tippoint(color = "black", size = tip_point_size) +
        xlim(-1, NA) +
        geom_label2(aes(subset = !isTip, label = gsub("\\|", "\\\n", label), fill = label),
                    color = "#3c3c3c", show.legend = FALSE, alpha = 0.6) +
        geom_tiplab(align = FALSE, offset = tip_label_offset, linesize = 0, size = tip_label_size, show.legend = FALSE) +
        scale_fill_manual(values = No_target_clone_mut_color)
    } else {
      target_clone_mut_color <- c("red", rep("white", length(pivot_muts) - 1))
      names(target_clone_mut_color) <- c(target_mut, setdiff(pivot_muts, target_mut))
      
      p1_tree <- ggtree(clone_order_tree, size = 0.1, branch.length = "none", ladderize = FALSE, layout = "radial", aes(color = target_clone)) %<+% leaf_df[, c("label", "target_clone")] +
        geom_tippoint(aes(color = target_clone), size = tip_point_size) +
        xlim(-1, NA) +
        scale_color_manual(values = c("clone_under_target_mut" = "red", "others" = "black"), guide = "none") +
        geom_label2(aes(subset = !isTip, label = gsub("\\|", "\\\n", label), fill = label),
                    color = "#3c3c3c", show.legend = FALSE, alpha = 0.6) +
        geom_tiplab(align = FALSE, offset = tip_label_offset, linesize = 0, size = tip_label_size, show.legend = FALSE) +
        scale_fill_manual(values = target_clone_mut_color, guide = "none")
    }
  } else {
    angle_per_cell <- 5
    total_angle <- angle_per_cell * num_cells
    if (total_angle > 360) total_angle <- 360
    open_angle <- 360 - total_angle
    
    if (target_mut == "no") {
      No_target_clone_mut_color <- rep("white", length(clone_order_tree_trim$node.label))
      names(No_target_clone_mut_color) <- sapply(clone_order_tree_trim$node.label, function(x) gsub("\\|", "\n", x))
      
      p1_tree <- ggtree(clone_order_tree_trim, size = 0.3, branch.length = "none", ladderize = FALSE, layout = "fan", color = "black", open.angle = open_angle) +
        geom_tippoint(color = "black", size = tip_point_size) +
        xlim(-1, NA) +
        geom_label2(aes(subset = !isTip, label = gsub("\\|", "\\\n", label), fill = label),
                    color = "#3c3c3c", show.legend = FALSE, alpha = 0.6) +
        geom_tiplab(align = FALSE, offset = tip_label_offset, linesize = 0, size = tip_label_size, show.legend = FALSE) +
        scale_fill_manual(values = No_target_clone_mut_color)
    } else {
      target_clone_mut_color <- c("red", rep("white", length(pivot_muts) - 1))
      names(target_clone_mut_color) <- c(target_mut, setdiff(pivot_muts, target_mut))
      
      p1_tree <- ggtree(clone_order_tree, size = 0.1, branch.length = "none", ladderize = FALSE, layout = "radial", aes(color = target_clone)) %<+% leaf_df[, c("label", "target_clone")] +
        geom_tippoint(aes(color = target_clone), size = tip_point_size) +
        xlim(-1, NA) +
        scale_color_manual(values = c("clone_under_target_mut" = "red", "others" = "black"), guide = "none") +
        geom_label2(aes(subset = !isTip, label = gsub("\\|", "\\\n", label), fill = label),
                    color = "#3c3c3c", show.legend = FALSE, alpha = 0.6) +
        geom_tiplab(align = FALSE, offset = tip_label_offset, linesize = 0, size = tip_label_size, show.legend = FALSE) +
        scale_fill_manual(values = target_clone_mut_color, guide = "none")
    }
  }
  
  # Add cell annotation layers
  p2_anno <- p1_tree
  
  if (exists("df_celltype")) {
    p2_anno <- p2_anno + new_scale_fill() +
      geom_fruit(data = df_celltype, geom = geom_tile, pwidth = 0.2, width = 0.2,
                 mapping = aes(y = label, x = "celltype", fill = celltype),
                 offset = 0.03, color = NA) +
      scale_fill_manual(values = df_celltype_color,
                        guide = guide_legend(title = "Cell Type", title.theme = title_theme,
                                             ncol = 3, title.position = "top",
                                             override.aes = list(shape = 18, size = 5), order = 1)) +
      theme(legend.direction = "horizontal", legend.box = "vertical", legend.box.just = "left",
            legend.title = element_text(size = 16), legend.text = element_text(size = 14))
  }
  
  if (exists("df_cluster")) {
    p2_anno <- p2_anno + new_scale_fill() +
      geom_fruit(data = df_cluster, geom = geom_tile, pwidth = 0.2, width = 0.2,
                 mapping = aes(y = label, x = "cluster", fill = cluster),
                 offset = 0.03, color = NA) +
      scale_fill_manual(values = df_cluster_color,
                        guide = guide_legend(title = "Cluster", title.theme = title_theme,
                                             ncol = 6, title.position = "top",
                                             override.aes = list(shape = 18, size = 5), order = 2)) +
      theme(legend.direction = "horizontal", legend.box = "vertical", legend.box.just = "left",
            legend.title = element_text(size = 16), legend.text = element_text(size = 14))
  }
  
  if (exists("df_sample")) {
    p2_anno <- p2_anno + new_scale_fill() +
      geom_fruit(data = df_sample, geom = geom_tile, pwidth = 0.2, width = 0.2,
                 mapping = aes(y = label, x = "sample", fill = sample),
                 offset = 0.03, color = NA) +
      scale_fill_manual(values = df_sample_color,
                        guide = guide_legend(title = "Sample", title.theme = title_theme,
                                             ncol = 4, title.position = "top",
                                             override.aes = list(shape = 18, size = 5), order = 3)) +
      theme(legend.direction = "horizontal", legend.box = "vertical", legend.box.just = "left",
            legend.title = element_text(size = 16), legend.text = element_text(size = 14))
  }
  
  if (exists("df_tumor")) {
    p2_anno <- p2_anno + new_scale_fill() +
      geom_fruit(data = df_tumor, geom = geom_tile, pwidth = 0.2, width = 0.2,
                 mapping = aes(y = label, x = "tumor_score", fill = tumor_score),
                 offset = 0.03, color = NA) +
      scale_fill_gradientn(colors = df_tumor_color, values = scales::rescale(tumor_breaks),
                           guide = guide_colorbar(title = "Tumor Score", title.theme = title_theme,
                                                  title.position = "top", barwidth = 10, barheight = 0.5, order = 4)) +
      theme(legend.direction = "horizontal", legend.box = "vertical", legend.box.just = "left",
            legend.title = element_text(size = 16), legend.text = element_text(size = 14))
  }
  
  if (exists("df_Bcell")) {
    p2_anno <- p2_anno + new_scale_fill() +
      geom_fruit(data = df_Bcell, geom = geom_tile, pwidth = 0.2, width = 0.2,
                 mapping = aes(y = label, x = "B_cell_prop", fill = B_cell_prop),
                 offset = 0.03, color = NA) +
      scale_fill_gradientn(colors = df_Bcell_color, values = scales::rescale(bcell_breaks),
                           guide = guide_colorbar(title = "B Cell Proportion", title.theme = title_theme,
                                                  title.position = "top", barwidth = 10, barheight = 0.5, order = 5)) +
      theme(legend.direction = "horizontal", legend.box = "vertical", legend.box.just = "left",
            legend.title = element_text(size = 16), legend.text = element_text(size = 14))
  }
  
  # Add heatmap
  guide_filp_point <- guide_legend(title = "Genotyper (raw > inferred)", title.position = "top",
                                   override.aes = list(size = 10), ncol = 1, title.theme = title_theme)
  
  if (target_mut == "no") {
    p3_heatmap <- p2_anno + new_scale_fill() +
      geom_fruit_list(
        geom_fruit(data = sub_long, geom = geom_tile, pwidth = heatmap_width,
                   mapping = aes(y = cellIDxmutID, x = factor(mutants_symbol, levels = pivot_muts_symbol),
                                 fill = as.character(mut_count)),
                   offset = heatmap_circos_offset, color = "white", size = 0.2,
                   axis.params = list(axis = "x", line.size = NA, line.color = NA, text.size = 0, vjust = 0.5)),
        scale_fill_manual(values = new_colors_before, breaks = names(new_colors_before)[!no_flipped_index],
                          labels = labels_dict[names(new_colors_before)[!no_flipped_index]], guide = guide_filp_point),
        new_scale_fill(),
        geom_fruit(data = sub_long, geom = geom_point, pwidth = heatmap_width,
                   mapping = aes(y = cellIDxmutID, x = factor(mutants_symbol, levels = pivot_muts_symbol),
                                 fill = as.character(mut_count)),
                   offset = heatmap_circos_offset, size = flipping_point_size, shape = 21, color = "transparent"),
        scale_fill_manual(values = new_colors_after, breaks = names(new_colors_after)[!no_flipped_index],
                          labels = labels_dict[names(new_colors_after)[!no_flipped_index]], guide = guide_filp_point)
      ) +
      theme(legend.text = element_text(size = 20), axis.title.x = element_blank(), axis.text.x = element_blank())
  } else {
    highlight_mutid <- sub_long %>%
      filter(mutants == gsub(":", ".", target_mut)) %>%
      pull(mutants_symbol) %>% unique() %>% as.character()
    
    p3_heatmap <- p2_anno + new_scale_fill() +
      geom_fruit_list(
        geom_fruit(data = sub_long, geom = geom_tile, pwidth = heatmap_width,
                   mapping = aes(y = cellIDxmutID, x = factor(mutants_symbol, levels = pivot_muts_symbol),
                                 fill = as.character(mut_count)),
                   offset = heatmap_circos_offset,
                   color = ifelse(factor(sub_long$mutants_symbol, levels = pivot_muts_symbol) == highlight_mutid, "black", "white"),
                   size = ifelse(factor(sub_long$mutants_symbol, levels = pivot_muts_symbol) == highlight_mutid, 0.2, 0.2),
                   axis.params = list(axis = "x", line.size = NA, line.color = NA, text.size = 0, vjust = 0.5)),
        scale_fill_manual(values = new_colors_before, breaks = names(new_colors_before)[!no_flipped_index],
                          labels = labels_dict[names(new_colors_before)[!no_flipped_index]], guide = guide_filp_point),
        new_scale_fill(),
        geom_fruit(data = sub_long, geom = geom_point, pwidth = heatmap_width,
                   mapping = aes(y = cellIDxmutID, x = factor(mutants_symbol, levels = pivot_muts_symbol),
                                 fill = as.character(mut_count)),
                   offset = heatmap_circos_offset, size = flipping_point_size, shape = 21, color = "transparent"),
        scale_fill_manual(values = new_colors_after, breaks = names(new_colors_after)[!no_flipped_index],
                          labels = labels_dict[names(new_colors_after)[!no_flipped_index]], guide = guide_filp_point)
      ) +
      theme(legend.text = element_text(size = 20), axis.title.x = element_blank(), axis.text.x = element_blank())
  }
  
  # Get ordered metadata for heatmap
  circostree <- p1_tree$data
  tip_data <- circostree[circostree$isTip, ]
  tip_data$angle <- atan2(tip_data$y, tip_data$x)
  sorted_tips <- tip_data[order(tip_data$angle), ]
  tip_label_order <- sorted_tips$label
  
  cols_to_sort <- for_sort_symbol$mutants
  
  metadata_for_heatmap_temp <- cf_tableno0
  sorting_metadata_for_heatmap <- metadata_for_heatmap_temp[, c("cellIDxmutID", cols_to_sort)]
  sorting_metadata_for_heatmap$cellIDxmutID <- factor(sorting_metadata_for_heatmap$cellIDxmutID, levels = tip_label_order)
  ordered_metadata_for_heatmap <- sorting_metadata_for_heatmap[order(sorting_metadata_for_heatmap$cellIDxmutID), ]
  ordered_metadata_for_heatmap[] <- lapply(ordered_metadata_for_heatmap, as.character)
  new_row <- c("mutid_in_heatmap", for_sort_symbol$mutants_symbol)
  ordered_metadata_for_heatmap_final <- rbind(new_row, ordered_metadata_for_heatmap)
  write.table(ordered_metadata_for_heatmap_final, file.path(outputpath, "ordered_metadata_for_heatmap.txt"), 
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  # Generate mutation id corresponding to the ordered tip.lables
  df_muts_corresponding_to_ordered_tiplabels <- data.frame(
    cellIDxmutID = ordered_metadata_for_heatmap_final$cellIDxmutID[2:nrow(ordered_metadata_for_heatmap_final)], 
    mutation_in_heatmap = character(nrow(ordered_metadata_for_heatmap_final) - 1)
  )
  
  for (i in 2:nrow(ordered_metadata_for_heatmap_final)) {
    row_data <- ordered_metadata_for_heatmap_final[i, cols_to_sort]
    ones_indices <- which(row_data == 1)
    if (length(ones_indices) > 0) {
      max_index <- max(ones_indices)
      mutation_least_cells <- cols_to_sort[max_index]
      df_muts_corresponding_to_ordered_tiplabels$mutation_in_heatmap[i - 1] <- mutation_least_cells
    } else {
      df_muts_corresponding_to_ordered_tiplabels$mutation_in_heatmap[i - 1] <- NA
    }
  }
  
  write.table(df_muts_corresponding_to_ordered_tiplabels, 
              file.path(outputpath, "df_muts_corresponding_to_ordered_tiplabels_by_anticlockwise.txt"), 
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  # Final plot
  main_plot <- p3_heatmap + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  pdf_lastfix <- paste0(format(Sys.time(), "%m%d_%H%M%S"), "_", substr(uuid::UUIDgenerate(), 1, 8))
  
  # Save main plot
  if (target_mut == "no") {
    svg_filename <- file.path(outputpath, paste0("No_target.circle_tree_output_as_point.", pdf_lastfix, ".svg"))
    pdf_filename <- file.path(outputpath, paste0("No_target.circle_tree_output_as_point.", pdf_lastfix, ".pdf"))
  } else {
    svg_filename <- file.path(outputpath, paste0("target_", target_mut, ".circle_tree_output_as_point.", pdf_lastfix, ".svg"))
    pdf_filename <- file.path(outputpath, paste0("target_", target_mut, ".circle_tree_output_as_point.", pdf_lastfix, ".pdf"))
  }
  
  ggplot2::ggsave(filename = svg_filename, plot = main_plot, 
                  width = plot_width * 2, height = plot_height * 2)
  ggplot2::ggsave(filename = pdf_filename, plot = main_plot, 
                  width = plot_width * 2, height = plot_height * 2)
  
  # Save circos annotation legend
  p_legend_circos <- cowplot::plot_grid(cowplot::get_legend(main_plot))
  legend_circos_svg <- file.path(outputpath, "legend_components.circos_annotation.svg")
  legend_circos_pdf <- file.path(outputpath, "legend_components.circos_annotation.pdf")
  ggplot2::ggsave(filename = legend_circos_svg, plot = p_legend_circos, 
                  width = plot_width, height = plot_height)
  ggplot2::ggsave(filename = legend_circos_pdf, plot = p_legend_circos, 
                  width = plot_width, height = plot_height)
  
  # Save total flipping count legend
  total_flipping <- read.table(total_flipping_file, header = TRUE, sep = "\t")
  colnames(total_flipping) <- c("total_flipping_0_to_1", "total_flipping_1_to_0", "total_flipping_NA_to_0", "total_flipping_NA_to_1")
  
  data_total_flipping <- tibble::tibble(
    type = c("False negative", "False positive", "Fill NA with 0", "Fill NA with 1"),
    category = c("0>1", "1>0", "NA>0", "NA>1"),
    count = c(
      total_flipping$total_flipping_0_to_1 %||% 0,
      total_flipping$total_flipping_1_to_0 %||% 0,
      total_flipping$total_flipping_NA_to_0 %||% 0,
      total_flipping$total_flipping_NA_to_1 %||% 0
    )
  )
  data_total_flipping$labels <- paste(data_total_flipping$category, ":", data_total_flipping$count)
  
  p_total_flipping <- ggplot(data_total_flipping, aes(x = labels, y = category, color = labels)) + 
    geom_text(aes(label = labels)) +
    scale_color_manual(values = rep("#3c3c3c", length(unique(data_total_flipping$labels))),
                       breaks = unique(data_total_flipping$labels),
                       guide = guide_legend(title = "Flipping count of all mutations",
                                            title.position = "top",
                                            override.aes = list(size = 8, label = data_total_flipping$type),
                                            ncol = 1,
                                            title.theme = title_theme)) +
    theme(legend.position = "right", legend.text = element_text(size = 20), legend.title = element_text())
  
  legend_total_flipping <- cowplot::get_legend(p_total_flipping)
  p_legend_total_flipping <- cowplot::plot_grid(legend_total_flipping)
  
  legend_flipping_svg <- file.path(outputpath, "legend_components.total_flipping_count.svg")
  legend_flipping_pdf <- file.path(outputpath, "legend_components.total_flipping_count.pdf")
  ggplot2::ggsave(filename = legend_flipping_svg, plot = p_legend_total_flipping, 
                  width = plot_width, height = plot_height)
  ggplot2::ggsave(filename = legend_flipping_pdf, plot = p_legend_total_flipping, 
                  width = plot_width, height = plot_height)
  
  # End timing
  end_time <- Sys.time()
  if (verbose) {
    message("\n", paste(rep("=", 60), collapse = ""))
    message("COMPLETED successfully!")
    message("Time elapsed: ", round(as.numeric(difftime(end_time, start_time, units = "secs")), 2), " seconds")
    message("Output directory: ", outputpath)
    message(paste(rep("=", 60), collapse = ""), "\n")
  }
  
  invisible(list(
    main_plot = main_plot,
    svg_file = svg_filename,
    pdf_file = pdf_filename
  ))
}