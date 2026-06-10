#' Sort mutation matrix using recursive clone partitioning
#'
#' Implements a recursive algorithm to sort mutations and cells based on
#' clonal relationships, optimizing the visualization layout.
#'
#' @param df Input data frame with mutations as columns and cells as rows
#' @return Sorted data frame
#' @export
#'
#' @examples
#' sorted_mat <- sort_mutation_matrix(raw_cf_table)
sort_mutation_matrix <- function(df) {
  mat <- as.matrix(df)
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)
  
  find_leader_clone <- function(available_items, by_row = TRUE) {
    if (length(available_items) == 0) return(NULL)
    
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
    
    if (length(members) == 0) return(NULL)
    return(list(leader = leader, members = members))
  }
  
  recursive_clone_partition <- function(available_items, parent_members = NULL, by_row = TRUE) {
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
      leader_clone <- find_leader_clone(available_items, by_row)
      if (is.null(leader_clone)) return(NULL)
      leader <- leader_clone$leader
      members <- leader_clone$members
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
      parallel_clone <- recursive_clone_partition(unused_items, by_row = by_row)
      if (!is.null(parallel_clone)) {
        current_clone$parallel <- parallel_clone
      }
    }
    
    return(current_clone)
  }
  
  extract_ordering <- function(clone_structure) {
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
  
  # Main sorting
  cell_clone_structure <- recursive_clone_partition(rownames(mat), by_row = TRUE)
  sorted_rows <- extract_ordering(cell_clone_structure)
  
  mut_clone_structure <- recursive_clone_partition(1:n_cols, by_row = FALSE)
  sorted_cols <- extract_ordering(mut_clone_structure)
  
  # Handle remaining items
  remaining_rows <- setdiff(rownames(mat), sorted_rows)
  if (length(remaining_rows) > 0) {
    remaining_rows <- remaining_rows[order(-rowSums(mat[remaining_rows, , drop = FALSE]))]
    sorted_rows <- c(sorted_rows, remaining_rows)
  }
  
  remaining_cols <- setdiff(1:n_cols, sorted_cols)
  if (length(remaining_cols) > 0) {
    remaining_cols <- remaining_cols[order(-colSums(mat[, remaining_cols, drop = FALSE]))]
    sorted_cols <- c(sorted_cols, remaining_cols)
  }
  
  return(df[sorted_rows, colnames(df)[sorted_cols], drop = FALSE])
}