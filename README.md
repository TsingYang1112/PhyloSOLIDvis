# PhyloSOLIDvis

**Visualization Suite for PhyloSOLID Phylogenetic Trees**

## Overview

PhyloSOLIDvis generates publication-ready circular (circos-style) phylogenetic tree visualizations with:

- ✨ **Circos-style circular plots** with integrated mutation heatmaps
- 🔬 **Cell annotation layers** (types, clusters, samples, tumor scores)
- 🧬 **Genotype flipping visualization** (input vs inferred states)
- 🎯 **Target mutation highlighting**
- 📊 **Adaptive parameter estimation**

## Installation

### Prerequisites

Make sure you have R (>= 4.0) installed.

### One-line installation (recommended)

```r
install.packages("remotes")
remotes::install_github("TsingYang1112/PhyloSOLIDvis")
```

### Alternative with devtools

```r
devtools::install_github("TsingYang1112/PhyloSOLIDvis")
```

### Manual dependency installation (if needed)

```r
# Install converTree from GitHub (required)
remotes::install_github("xiayh17/converTree")

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ggtree", "ggtreeExtra", "treeio", "ComplexHeatmap"))
```

## Quick Start

```r
library(PhyloSOLIDvis)

result <- plot_circos(
  inputpath = "path/to/PhyloSOLID/output/",
  outputpath = "path/to/figures/",
  annotation_file = "path/to/annotations.txt"
)
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `inputpath` | (required) | Path to directory containing PhyloSOLID output files |
| `outputpath` | (required) | Path to output directory for saving plots |
| `annotation_file` | (required) | Path to cell annotation file (TSV format) |
| `target_mut` | "no" | Target mutation ID to highlight |
| `tip_label_offset` | 10 | Offset distance for tip labels |
| `tip_label_size` | 2.5 | Font size for tip labels |
| `tip_point_size` | 0.5 | Size of tip points |
| `heatmap_width` | 0.3 | Width of heatmap track |
| `heatmap_circos_offset` | 0.04 | Offset for heatmap from tree |
| `flipping_point_size` | 0.2 | Size of flipping marker points |
| `plot_height` | 12 | Output plot height (inches) |
| `plot_width` | 18 | Output plot width (inches) |
| `verbose` | TRUE | Print progress messages |

## Required Input Files

| File | Description |
|------|-------------|
| `final_cleaned_I_full_withNA3_for_circosPlot.txt` | Input genotype matrix |
| `final_cleaned_M_full_basedPivots.filtered_sites_inferred.CFMatrix` | Conflict-free matrix |
| `df_flipping_count_for_each_mut.txt` | Per-mutation flipping statistics |
| `df_total_flipping_count.txt` | Total flipping statistics |

### Annotation File Format

The annotation file (TSV format) should contain a `barcode` column and can include any of the following optional columns:

| Column | Description |
|--------|-------------|
| `cluster_info` | Cell cluster assignments |
| `cell_type` | Cell type annotations |
| `sample` | Sample identifiers |
| `tumor_score` | Numeric tumor scores |
| `B_cell_prop` | B cell proportions |

## Output Files

| File | Description |
|------|-------------|
| `target_*.svg/pdf` | Main circular plot |
| `legend_components.circos_annotation.svg/pdf` | Legend |
| `sorted_cf_matrix.txt` | Sorted mutation matrix |
| `clone_and_subclone_table.txt` | Clone assignments |
| `ordered_metadata_for_heatmap.txt` | Heatmap ordering data |
| `df_muts_corresponding_to_ordered_tiplabels_by_anticlockwise.txt` | Mutation-cell mapping |

## Example

```r
# Basic usage
library(PhyloSOLIDvis)

result <- plot_circos(
  inputpath = "results/phylosolid/",
  outputpath = "figures/circos/",
  annotation_file = "data/annotations.txt",
  verbose = TRUE
)

# With target mutation highlighting
result <- plot_circos(
  inputpath = "results/phylosolid/",
  outputpath = "figures/circos/",
  annotation_file = "data/annotations.txt",
  target_mut = "chr11_65426524_T_C",
  tip_label_offset = 12,
  tip_label_size = 3
)

# Check output files
print(result$output_files)
```

## Troubleshooting

### Network issues during installation

If you encounter network errors, try using a CRAN mirror:

```r
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
remotes::install_github("TsingYang1112/PhyloSOLIDvis")
```

### Missing dependencies

If dependencies fail to install, install them manually:

```r
# CRAN packages
install.packages(c("ape", "circlize", "cowplot", "dplyr", "ggplot2"))

# Bioconductor packages
BiocManager::install(c("ggtree", "ggtreeExtra", "treeio", "ComplexHeatmap"))

# GitHub packages
remotes::install_github("xiayh17/converTree")
```

## Acknowledgments

The author thanks **Yonghe Xia** for initial conceptualization of the visualization approach and for providing the [converTree](https://github.com/xiayh17/converTree) package, which is used for tree construction from conflict-free matrices.

## Author

**Qing Yang**  
*Westlake University, School of Life Sciences*  
Email: qing.yang@westlake.edu.cn  
ORCID: 0009-0005-7366-6879  
GitHub: [@TsingYang1112](https://github.com/TsingYang1112)

## License

MIT © Westlake University

## Citation

If you use PhyloSOLIDvis in your research, please cite:

```bibtex
@software{Yang_PhyloSOLIDvis_2024,
  author = {Qing Yang},
  title = {PhyloSOLIDvis: Visualization Suite for PhyloSOLID Phylogenetic Trees},
  year = {2024},
  publisher = {GitHub},
  url = {https://github.com/TsingYang1112/PhyloSOLIDvis}
}
```
