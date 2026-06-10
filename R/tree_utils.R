#' Convert CF matrix to tree data structure
#'
#' Converts a conflict-free matrix to a tree data structure compatible with ggtree.
#'
#' @param cffile Path to CF matrix file
#' @return Tree data object compatible with ggtree
#' @export
#'
#' @examples
#' tree_data <- cf2treedata("path/to/CFMatrix.txt")
cf2treedata <- function(cffile) {
  tab <- read.table(cffile, header = TRUE, sep = "\t", row.names = 1)
  tab[tab > 0] <- 1
  
  dist_mat <- dist(tab, method = "binary")
  tree <- hclust(dist_mat, method = "ward.D2")
  phylo_tree <- ape::as.phylo(tree)
  tree_data <- treeio::as.treedata(phylo_tree)
  
  # Add mutation labels to nodes
  for (i in 1:nrow(tree_data@data)) {
    if (!tree_data@data$isTip[i]) {
      tips <- treeio::offspring(tree_data, i, tiponly = TRUE)
      if (length(tips) > 0) {
        node_muts <- colnames(tab)[colSums(tab[tips, , drop = FALSE]) > 0]
        tree_data@data$label[i] <- paste(node_muts, collapse = "|")
      }
    }
  }
  
  return(tree_data)
}

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
    tmp <- stringr::str_count(dat$label, fixed("|"))
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