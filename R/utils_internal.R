#' Remove zero barcodes from CF matrix
#' @keywords internal
remove_zero_barcode <- function(cffile, output_dir = NULL) {
  # Determine output file path
  if (!is.null(output_dir)) {
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    base_name <- basename(cffile)
    cffile_no0 <- file.path(output_dir, paste0(base_name, ".no0"))
  } else {
    # Original behavior: save in same directory
    cffile_no0 <- paste0(cffile, ".no0")
  }
  
  lines <- readLines(cffile)
  filtered_lines <- suppressWarnings(
    lines[c(TRUE, sapply(lines[-1], function(line) {
      sum(as.numeric(strsplit(line, "\t")[[1]][-1])) != 0
    }))]
  )
  writeLines(filtered_lines, cffile_no0)
  return(cffile_no0)
}

#' Process annotations and generate colors (with zero barcodes filtering)
#' @keywords internal
process_annotations_with_zero <- function(df_all_info, zero_barcodes) {
  result <- list()
  
  if ("cluster_info" %in% colnames(df_all_info)) {
    temp_cluster <- df_all_info[, c("barcode", "cluster_info")]
    temp_cluster$cluster_info <- trimws(as.character(temp_cluster$cluster_info))
    temp_cluster <- temp_cluster[!temp_cluster$barcode %in% zero_barcodes, ]
    temp_cluster <- temp_cluster[!is.na(temp_cluster$cluster_info) & temp_cluster$cluster_info != "", ]
    
    if (nrow(temp_cluster) > 0) {
      df_cluster <- data.frame(label = temp_cluster$barcode, cluster = temp_cluster$cluster_info, stringsAsFactors = FALSE)
      cluster_levels <- sort(unique(df_cluster$cluster))
      df_cluster$cluster <- factor(df_cluster$cluster, levels = cluster_levels)
      set.seed(42)
      df_cluster_color <- grDevices::rgb(
        runif(length(cluster_levels)), 
        runif(length(cluster_levels)), 
        runif(length(cluster_levels))
      )
      names(df_cluster_color) <- cluster_levels
      result$cluster <- list(data = df_cluster, colors = df_cluster_color)
    }
  }
  
  if ("cell_type" %in% colnames(df_all_info)) {
    temp_celltype <- df_all_info[, c("barcode", "cell_type")]
    temp_celltype$cell_type <- trimws(as.character(temp_celltype$cell_type))
    temp_celltype <- temp_celltype[!temp_celltype$barcode %in% zero_barcodes, ]
    temp_celltype <- temp_celltype[!is.na(temp_celltype$cell_type) & temp_celltype$cell_type != "", ]
    
    if (nrow(temp_celltype) > 0) {
      df_celltype <- data.frame(label = temp_celltype$barcode, celltype = temp_celltype$cell_type, stringsAsFactors = FALSE)
      celltype_levels <- sort(unique(df_celltype$celltype))
      df_celltype$celltype <- factor(df_celltype$celltype, levels = celltype_levels)
      set.seed(142)
      df_celltype_color <- grDevices::rgb(
        runif(length(celltype_levels)), 
        runif(length(celltype_levels)), 
        runif(length(celltype_levels))
      )
      names(df_celltype_color) <- celltype_levels
      result$celltype <- list(data = df_celltype, colors = df_celltype_color)
    }
  }
  
  if ("sample" %in% colnames(df_all_info)) {
    temp_sample <- df_all_info[, c("barcode", "sample")]
    temp_sample$sample <- trimws(as.character(temp_sample$sample))
    temp_sample <- temp_sample[!temp_sample$barcode %in% zero_barcodes, ]
    temp_sample <- temp_sample[!is.na(temp_sample$sample) & temp_sample$sample != "", ]
    
    if (nrow(temp_sample) > 0) {
      df_sample <- data.frame(label = temp_sample$barcode, sample = temp_sample$sample, stringsAsFactors = FALSE)
      sample_levels <- unique(df_sample$sample)
      df_sample_color <- c("#ebb16e", "#b2aad3")[seq_along(sample_levels)]
      names(df_sample_color) <- sample_levels
      result$sample <- list(data = df_sample, colors = df_sample_color)
    }
  }
  
  if ("tumor_score" %in% colnames(df_all_info)) {
    temp_tumor <- df_all_info[, c("barcode", "tumor_score")]
    temp_tumor$tumor_score <- as.numeric(as.character(temp_tumor$tumor_score))
    temp_tumor <- temp_tumor[!temp_tumor$barcode %in% zero_barcodes, ]
    temp_tumor <- temp_tumor[!is.na(temp_tumor$tumor_score), ]
    
    if (nrow(temp_tumor) > 0) {
      df_tumor <- data.frame(label = temp_tumor$barcode, tumor_score = temp_tumor$tumor_score, stringsAsFactors = FALSE)
      min_value <- min(df_tumor$tumor_score[df_tumor$tumor_score != 0], na.rm = TRUE)
      max_value <- max(df_tumor$tumor_score, na.rm = TRUE)
      result$tumor <- list(
        data = df_tumor,
        colors = c("#00a3c4", "#ffebf6", "#ff64be"),
        breaks = c(0, min_value, max_value)
      )
    }
  }
  
  if ("B_cell_prop" %in% colnames(df_all_info)) {
    temp_bcell <- df_all_info[, c("barcode", "B_cell_prop")]
    temp_bcell$B_cell_prop <- as.numeric(as.character(temp_bcell$B_cell_prop))
    temp_bcell <- temp_bcell[!temp_bcell$barcode %in% zero_barcodes, ]
    temp_bcell <- temp_bcell[!is.na(temp_bcell$B_cell_prop), ]
    
    if (nrow(temp_bcell) > 0) {
      df_Bcell <- data.frame(label = temp_bcell$barcode, B_cell_prop = temp_bcell$B_cell_prop, stringsAsFactors = FALSE)
      min_value <- min(df_Bcell$B_cell_prop[df_Bcell$B_cell_prop != 0], na.rm = TRUE)
      max_value <- max(df_Bcell$B_cell_prop, na.rm = TRUE)
      median_value <- median(df_Bcell$B_cell_prop, na.rm = TRUE)
      result$Bcell <- list(
        data = df_Bcell,
        colors = c("#F4D166", "#DBE8B4", "#24693D"),
        breaks = c(min_value, median_value, max_value)
      )
    }
  }
  
  return(result)
}

#' Add target clone information
#' @keywords internal
add_target_clone <- function(leaf_df, cf_tableno0_sorted, target_mut) {
  target_clone <- cf_tableno0_sorted[cf_tableno0_sorted[[target_mut]] == 1, "cellIDxmutID"]
  leaf_df <- leaf_df %>% 
    dplyr::mutate(
      target_clone = ifelse(label %in% target_clone, "clone_under_target_mut", "others")
    )
  return(leaf_df)
}

#' Sort leaf data frame
#' @keywords internal
sort_leaf_df <- function(leaf_df, cf_tableno0_sorted) {
  mutation_mat <- cf_tableno0_sorted[, -1, drop = FALSE]
  leaf_df$mutation_count <- rowSums(mutation_mat[leaf_df$label, , drop = FALSE])
  
  leaf_df <- leaf_df %>%
    dplyr::arrange(subclone_size, clone_order) %>%
    dplyr::group_by(clone_order) %>%
    dplyr::arrange(dplyr::desc(mutation_count), .by_group = TRUE) %>%
    dplyr::ungroup()
  
  return(leaf_df)
}

#' Recursive clone partition
#' @keywords internal
recursive_clone_partition <- function(available_items, mat, parent_members = NULL, by_row = TRUE) {
  if (length(available_items) == 0) return(NULL)
  
  if (!is.null(parent_members)) {
    if (by_row) {
      mat_subset <- mat[available_items, parent_members, drop = FALSE]
      row_sums <- rowSums(mat_subset)
      if (max(row_sums) == 0) return(NULL)
      leader_idx <- which.max(row_sums)
      leader <- available_items[leader_idx]
      members <- parent_members[which(mat_subset[leader_idx, ] == 1)]
    } else {
      mat_subset <- mat[parent_members, available_items, drop = FALSE]
      col_sums <- colSums(mat_subset)
      if (max(col_sums) == 0) return(NULL)
      leader_idx <- which.max(col_sums)
      leader <- available_items[leader_idx]
      members <- parent_members[which(mat_subset[, leader_idx] == 1)]
    }
  } else {
    if (by_row) {
      item_sums <- rowSums(mat[available_items, , drop = FALSE])
      leader_idx <- which.max(item_sums)
      leader <- available_items[leader_idx]
      members <- which(mat[leader, ] == 1)
    } else {
      item_sums <- colSums(mat[, available_items, drop = FALSE])
      leader_idx <- which.max(item_sums)
      leader <- available_items[leader_idx]
      members <- which(mat[, leader] == 1)
    }
  }
  
  if (by_row) {
    clone_items <- available_items[rowSums(mat[available_items, members, drop = FALSE]) > 0]
    clone_items <- clone_items[order(rowSums(mat[clone_items, members, drop = FALSE]), decreasing = TRUE)]
  } else {
    clone_items <- available_items[colSums(mat[members, available_items, drop = FALSE]) > 0]
    clone_items <- clone_items[order(colSums(mat[members, clone_items, drop = FALSE]), decreasing = TRUE)]
  }
  
  current_clone <- list(leader = leader, clone_items = clone_items, members = members)
  
  unused_items <- setdiff(available_items, clone_items)
  if (length(unused_items) > 0) {
    parallel_clone <- recursive_clone_partition(unused_items, mat, by_row = by_row)
    if (!is.null(parallel_clone)) {
      current_clone$parallel <- parallel_clone
    }
  }
  
  return(current_clone)
}

#' Extract ordering from clone structure
#' @keywords internal
extract_mutation_order <- function(clone_structure) {
  ordered_items <- c()
  
  process_clone <- function(clone) {
    if (is.null(clone)) return()
    new_items <- setdiff(clone$clone_items, ordered_items)
    if (length(new_items) > 0) {
      ordered_items <<- c(ordered_items, new_items)
    }
    if (!is.null(clone$parallel)) {
      process_clone(clone$parallel)
    }
  }
  
  process_clone(clone_structure)
  return(ordered_items)
}

#' Format flipping label for legend
#' @keywords internal
format_flipping_label <- function(x) {
  x <- gsub("0", "non-mutant", x)
  x <- gsub("1", "mutant", x)
  x <- gsub("3", "missing", x)
  x <- gsub(">", " -> ", x)
  return(x)
}

#' Insert newlines into long text
#' @keywords internal
insert_newline <- function(text, n = 180) {
  split_text <- strsplit(text, NULL)[[1]]
  split_text <- split(split_text, ceiling(seq_along(split_text) / n))
  result <- sapply(split_text, paste, collapse = "")
  final_result <- paste(result, collapse = "\n")
  return(final_result)
}

#' Save outputs
#' @keywords internal
save_outputs <- function(plot_result, outputpath, target_mut, plot_width, plot_height) {
  pdf_lastfix <- paste0(format(Sys.time(), "%m%d_%H%M%S"), "_", substr(uuid::UUIDgenerate(), 1, 8))
  
  if (target_mut == "no") {
    svg_filename <- file.path(outputpath, paste0("No_target.circle_tree_output_as_point.", pdf_lastfix, ".svg"))
  } else {
    svg_filename <- file.path(outputpath, paste0("target_", target_mut, ".circle_tree_output_as_point.", pdf_lastfix, ".svg"))
  }
  
  pdf_filename <- gsub("svg$", "pdf", svg_filename)
  
  ggplot2::ggsave(filename = svg_filename, plot = plot_result$main_plot, 
         width = plot_width * 2, height = plot_height * 2)
  ggplot2::ggsave(filename = pdf_filename, plot = plot_result$main_plot, 
         width = plot_width * 2, height = plot_height * 2)
  
  return(c(svg_filename, pdf_filename))
}

#' Find zero barcodes (rows with all zeros)
#' @keywords internal
find_zero_barcode <- function(table) {
  zero_barcode <- rownames(table)[rowSums(table) == 0]
  return(zero_barcode)
}

#' Find zero mutations (columns with all zeros)
#' @keywords internal
find_zero_mutation <- function(table) {
  zero_mutation <- colnames(table)[colSums(table) == 0]
  return(zero_mutation)
}

#' Trim node label for display
#' @keywords internal
trim_label <- function(label_str, target_mut, default_n = 2) {
  if (is.na(label_str) || label_str == "") return("")
  
  parts <- unlist(strsplit(label_str, "\\|"))
  
  if (target_mut %in% parts) {
    parts_trim <- target_mut
  } else {
    parts_trim <- head(parts, default_n)
  }
  
  paste(parts_trim, collapse = "\n")
}

#' Null coalescing operator
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Process annotations for plot (backward compatible)
#' @keywords internal
process_annotations <- function(df_all_info) {
  # This is kept for backward compatibility
  process_annotations_with_zero(df_all_info, c())
}