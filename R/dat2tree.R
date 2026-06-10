#' Convert data frame to phylogenetic tree with branch lengths
#'
#' Converts a tree data frame to a phylo object where branch lengths represent mutation counts.
#'
#' @param dat Tree data frame from converTree or phylo object
#' @return Phylogenetic tree object (phylo class)
#' @export
#'
#' @examples
#' tree <- dat2tree(tree_data_frame)
dat2tree <- function(dat) {
  # If dat is already a phylo object, return it directly
  if (inherits(dat, "phylo")) {
    return(dat)
  }
  
  # Check if dat has the expected columns for conversion
  if (!"parent" %in% colnames(dat) || !"node" %in% colnames(dat)) {
    # Try direct conversion to phylo
    tree_ <- tryCatch({
      ape::as.phylo(dat)
    }, error = function(e) {
      stop("Cannot convert to phylo object. Expected a phylo object or a data frame with 'parent' and 'node' columns.")
    })
    return(tree_)
  }
  
  # Original conversion logic for data frame with parent/node columns
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
    root_node <- dat[dat$node == dat$parent, ]$parent[1]
    dat3 <- dat[dat$node != root_node, ]
    dat3$branch.length
  }
  
  tree_[["edge.length"]] <- index_length()
  tree_[["tip.label"]] <- sapply(tree_[["tip.label"]], index_label) %>% unlist()
  tree_[["node.label"]] <- sapply(tree_[["node.label"]], index_label) %>% unlist()
  
  return(tree_)
}