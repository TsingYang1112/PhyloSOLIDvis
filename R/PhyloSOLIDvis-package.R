#' PhyloSOLIDvis: Visualization Suite for PhyloSOLID Phylogenetic Trees
#'
#' @description
#' Generate publication-ready circular phylogenetic tree visualizations with 
#' integrated mutation heatmaps, cell annotations, and genotype flipping markers.
#'
#' @docType package
#' @name PhyloSOLIDvis
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import ggtree
#' @import ggtreeExtra
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import converTree
#' @import svglite
#' @importFrom ape as.phylo rotateConstr Ntip Nnode
#' @importFrom cowplot get_legend plot_grid
#' @importFrom ggnewscale new_scale_fill
#' @importFrom grDevices rgb colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales rescale
#' @importFrom stringr str_c str_count str_split str_detect
#' @importFrom tibble column_to_rownames as_tibble
#' @importFrom utils packageVersion read.table write.table
#' @importFrom stats na.omit
#' @importFrom pheatmap pheatmap
#' @importFrom ggtree ggtree geom_tree theme_tree2 get_taxa_name
#' @importFrom ggplot2 theme element_blank margin coord_cartesian
#' @export run_all
#' @export prepare_data
#' @export plot_circos
#' @export plot_heatmap
#' @export order_matrix
#' @export sort_mutation_matrix
#' @export subclone_finder
#' @export dat2tree
#' @export is_PhyloData
#' @export print.PhyloData

## usethis namespace: end
NULL

#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}