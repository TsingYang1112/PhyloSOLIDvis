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
#' @importFrom ape as.phylo rotateConstr
#' @importFrom cowplot get_legend plot_grid
#' @importFrom ggnewscale new_scale_fill
#' @importFrom stringr str_c str_count str_split
#' @importFrom utils packageVersion read.table write.table
## usethis namespace: end
NULL