#' Convert data frame to phylogenetic tree with branch lengths
#'
#' Converts a tree data frame to a phylo object where branch lengths represent mutation counts.
#'
#' @param dat Tree data frame with node, parent, and label columns
#' @return Phylogenetic tree object (phylo class)
#' @export
#'
#' @examples
#' tree <- dat2tree(tree_data_frame)
dat2tree <- function(dat) {
  count_muts <- function(dat) {
    tmp <- stringr::str_count(dat$label, "\\|")
    return(tmp + 1)
  }
  
  mutation_count <- count_muts(dat)
  dat$branch.length <- mutation_count
  
  tree_ <- dat %>% ape::as.phylo()
  
  index_label <- function(node) {
    dat[dat$node == node, ]$label
  }
  
  index_length <- function() {
    root_node <- dat[dat$node == dat$parent, ]$parent
    dat3 <- dat[dat$node != root_node, ]
    dat3$branch.length
  }
  
  tree_[["edge.length"]] <- index_length()
  tree_[["tip.label"]] <- sapply(tree_[["tip.label"]], index_label) %>% unlist()
  tree_[["node.label"]] <- sapply(tree_[["node.label"]], index_label) %>% unlist()
  
  return(tree_)
}