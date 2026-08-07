#' PhyloData Class for Storing Phylogenetic Data
#'
#' An S3 class that stores all phylogenetic data after one-time processing.
#' This avoids repeated file I/O and computation across different visualization functions.
#'
#' @keywords internal
PhyloData <- function(
    raw_cf = NULL,
    input_table = NULL,
    sorted_cf = NULL,
    tree = NULL,
    tree_data = NULL,
    tree_no0 = NULL,
    zero_barcodes = NULL,
    annotations = NULL,
    clone_data = NULL,
    flipping_data = NULL,
    total_flipping = NULL
) {
  structure(
    list(
      raw_cf = raw_cf,
      input_table = input_table,
      sorted_cf = sorted_cf,
      tree = tree,
      tree_data = tree_data,
      tree_no0 = tree_no0,
      zero_barcodes = zero_barcodes,
      annotations = annotations,
      clone_data = clone_data,
      flipping_data = flipping_data,
      total_flipping = total_flipping
    ),
    class = "PhyloData"
  )
}

#' Check if object is PhyloData
#' @export
is_PhyloData <- function(x) {
  inherits(x, "PhyloData")
}

#' Print method for PhyloData
#' @export
print.PhyloData <- function(x, ...) {
  cat("PhyloData object\n")
  cat("----------------\n")
  cat("Raw CF matrix:", if (!is.null(x$raw_cf)) paste0(nrow(x$raw_cf), "x", ncol(x$raw_cf)) else "NULL", "\n")
  cat("Sorted CF matrix:", if (!is.null(x$sorted_cf)) paste0(nrow(x$sorted_cf), "x", ncol(x$sorted_cf)) else "NULL", "\n")
  cat("Tree:", if (!is.null(x$tree)) paste0(ape::Ntip(x$tree), " tips, ", ape::Nnode(x$tree), " nodes") else "NULL", "\n")
  cat("Tree data:", if (!is.null(x$tree_data)) paste0(nrow(x$tree_data), " rows") else "NULL", "\n")
  cat("Tree no0:", if (!is.null(x$tree_no0)) paste0(nrow(x$tree_no0), " rows") else "NULL", "\n")
  cat("Zero barcodes:", if (!is.null(x$zero_barcodes)) length(x$zero_barcodes) else "NULL", "\n")
  cat("Annotations:", if (!is.null(x$annotations)) paste(names(x$annotations), collapse = ", ") else "None", "\n")
  cat("Clone data:", if (!is.null(x$clone_data)) paste0(nrow(x$clone_data), " rows") else "NULL", "\n")
  cat("Flipping data:", if (!is.null(x$flipping_data)) paste0(nrow(x$flipping_data), " rows") else "NULL", "\n")
  invisible(x)
}