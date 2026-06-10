#' PhyloSOLIDvis: Visualization Suite for PhyloSOLID Phylogenetic Trees
#'
#' @description
#' Generate publication-ready circular phylogenetic tree visualizations with 
#' integrated mutation heatmaps, cell annotations, and genotype flipping markers.
#'
#' @section Main Features:
#' \itemize{
#'   \item \strong{Circos Plot}: Circular tree visualization with mutation heatmaps
#'   \item \strong{Clone Partitioning}: Recursive algorithm for optimal ordering
#'   \item \strong{Annotation Integration}: Cell types, clusters, samples, scores
#'   \item \strong{Genotype Flipping}: Visualize input vs inferred discrepancies
#'   \item \strong{Adaptive Parameters}: Auto-estimation based on data dimensions
#' }
#'
#' @section Acknowledgments:
#' The author thanks Yonghe Xia for initial conceptualization of the visualization 
#' approach and for providing the converTree package, which is used for tree 
#' construction from conflict-free matrices.
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
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom cowplot get_legend plot_grid
#' @importFrom ggforce geom_arc_bar
#' @importFrom ggnewscale new_scale_fill
#' @importFrom ggtext element_markdown
#' @importFrom grid gpar unit
#' @importFrom gridExtra grid.arrange
#' @importFrom gsubfn gsubfn
#' @importFrom paletteer paletteer_d
#' @importFrom pals viridis
#' @importFrom patchwork plot_layout
#' @importFrom pheatmap pheatmap
#' @importFrom phylogram as.dendrogram
#' @importFrom Polychrome palette36.colors
#' @importFrom rlang `%||%`
#' @importFrom stringr str_c str_count str_detect str_replace_all str_split
#' @importFrom tibble as_tibble
#' @importFrom treeio as.treedata offspring
#' @importFrom utils packageVersion read.table write.table
## usethis namespace: end
NULL

#' Check and warn about converTree installation
#' 
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  if (!requireNamespace("converTree", quietly = TRUE)) {
    packageStartupMessage(
      "Note: converTree package is required for tree construction.\n",
      "It will be automatically installed from GitHub when needed.\n",
      "Or manually install with: devtools::install_github('xiayh17/converTree')"
    )
  }
}