#' Identify subclones and clones from conflict-free matrix
#'
#' Analyzes a conflict-free mutation matrix to identify subclones and their 
#' hierarchical relationships using a recursive superset algorithm.
#'
#' @param cf_tableno0 Data frame of conflict-free matrix (first column is cell IDs)
#' @return Data frame with subclone information including clone assignment,
#'         clone size, and clone order
#' @export
#'
#' @examples
#' subclones <- subclone_finder(cf_matrix)
subclone_finder <- function(cf_tableno0) {
  # Convert to long format
  cf_tableno0_long <- cf_tableno0 %>%
    tidyr::gather(key = "mutation", value = "presence", -cellIDxmutID) %>%
    dplyr::filter(presence == 1) %>%
    dplyr::arrange(cellIDxmutID, mutation)
  
  # Reorder mutations by frequency
  muts_factor <- factor(cf_tableno0_long$mutation)
  level_counts <- table(muts_factor)
  ordered_levels <- names(sort(level_counts, decreasing = TRUE))
  muts_factor_reordered <- factor(muts_factor, levels = ordered_levels)
  cf_tableno0_long$mutation_num <- as.character(as.numeric(muts_factor_reordered))
  
  # Identify subclones
  subclones <- cf_tableno0_long %>%
    dplyr::group_by(cellIDxmutID) %>%
    dplyr::summarize(
      subclone = paste(mutation, collapse = ","),
      subclone_num = paste(mutation_num, collapse = ","),
      .groups = "drop"
    ) %>%
    dplyr::group_by(subclone, subclone_num) %>%
    dplyr::summarize(
      cells = paste(cellIDxmutID, collapse = ","),
      .groups = "drop"
    )
  
  subclones$subclone_size <- stringr::str_count(subclones$subclone, ",") + 1
  
  # Convert subclone_num to integer lists
  subclones$subclone_int <- lapply(subclones$subclone_num, function(x) {
    sort(as.integer(unlist(strsplit(x, ","))))
  })
  
  # Superset check function
  is_superset <- function(subclone_set, other_set) {
    all(other_set %in% subclone_set)
  }
  
  # Merge subsets recursively
  merge_subsets <- function(clone_mapping) {
    to_remove <- rep(FALSE, length(clone_mapping))
    for (i in seq_along(clone_mapping)) {
      for (j in seq_along(clone_mapping)) {
        if (i != j && !to_remove[j] && is_superset(clone_mapping[[i]], clone_mapping[[j]])) {
          clone_mapping[[i]] <- unique(c(clone_mapping[[i]], clone_mapping[[j]]))
          to_remove[j] <- TRUE
        }
      }
    }
    return(clone_mapping[!to_remove])
  }
  
  # Create clone mapping
  clone_mapping <- setNames(subclones$subclone_int, paste0("clone", seq_along(subclones$subclone_int)))
  clone_mapping <- merge_subsets(clone_mapping)
  
  # Map to subclones
  subclones$clone_list <- sapply(subclones$subclone_int, function(x) {
    matching_clone <- which(sapply(clone_mapping, function(y) is_superset(y, x)))
    if (length(matching_clone) == 0) return(NA)
    return(names(clone_mapping)[matching_clone])
  })
  
  # Rename clones
  subclones$count <- as.character(subclones$clone_list) %>% stringr::str_count(",") + 1
  subclones$clone <- NA
  subclones$clone[subclones$count == 1] <- unlist(subclones$clone_list[subclones$count == 1])
  
  contatins_all <- sapply(subclones$clone_list, function(x) {
    return(all(names(clone_mapping) %in% x))
  })
  subclones$clone[subclones$count != 1 & subclones$count == max(subclones$count) & contatins_all] <- "root_clone"
  
  subroot_who <- sapply(subclones$clone_list[subclones$count != 1 & subclones$count < max(subclones$count)], 
                         function(x) paste0("subroot:(", paste(x, collapse = ","), ")"))
  subclones$clone[subclones$count != 1 & subclones$count < max(subclones$count)] <- subroot_who
  
  # Order subclones
  subroot_list <- subclones$clone_list
  names(subroot_list) <- rownames(subclones)
  
  contain_matrix <- matrix(FALSE, nrow = length(subroot_list), ncol = length(subroot_list))
  for (i in 1:length(subroot_list)) {
    for (j in 1:length(subroot_list)) {
      if (i != j && all(subroot_list[[i]] %in% subroot_list[[j]])) {
        contain_matrix[i, j] <- TRUE
      }
    }
  }
  
  lengths <- sapply(subroot_list, length)
  order_index <- order(-rowSums(contain_matrix), -lengths)
  sorted_list <- subroot_list[order_index]
  subclones <- subclones[names(sorted_list), ]
  subclones$clone_order <- 1:nrow(subclones)
  
  # Re-number clones
  clones <- subclones$clone
  numbers <- gsub(".*?(\\d+)$", "\\1", clones)
  numbers <- as.numeric(numbers)
  unique_numbers <- unique(na.omit(numbers))
  new_numbers <- setNames(seq_along(unique_numbers), unique_numbers)
  
  match <- regexpr("\\d+$", clones, perl = TRUE)
  last_number <- regmatches(clones, match)
  subclones$clone[which(attr(match, "match.length") != -1)] <- paste0("clone", new_numbers[last_number])
  
  test <- subclones$clone[subclones$count != 1 & subclones$count < max(subclones$count)]
  if (length(test) > 0) {
    replacement_function <- function(x) paste0("clone", new_numbers[x])
    test2 <- gsubfn::gsubfn("clone(\\d+)", replacement_function, test)
    subclones$clone[subclones$count != 1 & subclones$count < max(subclones$count)] <- test2
  }
  
  # Final summary
  subclones <- subclones %>%
    dplyr::group_by(clone) %>%
    dplyr::reframe(clone_size = sum(subclone_size), dplyr::across(dplyr::everything()))
  
  return(subclones)
}