# PhyloSOLIDvis

**Visualization Suite for PhyloSOLID Phylogenetic Trees**

## Overview

PhyloSOLIDvis generates publication-ready circular (circos-style) phylogenetic tree visualizations with:

- ✨ **Circos-style circular plots** with integrated mutation heatmaps
- 🔬 **Cell annotation layers** (types, clusters, samples, tumor scores, B cell proportions)
- 🧬 **Genotype flipping visualization** (input vs inferred states)
- 🎯 **Target mutation highlighting**
- 📊 **Adaptive parameter estimation**

## Installation

### One-line installation (recommended)

```r
remotes::install_github("TsingYang1112/PhyloSOLIDvis", dependencies = TRUE)
```

### Alternative with devtools

```r
devtools::install_github("TsingYang1112/PhyloSOLIDvis", dependencies = TRUE)
```

## Quick Start

**Note:** `ggplot2` must be loaded before `PhyloSOLIDvis` to avoid namespace conflicts.

```r
library(ggplot2)
library(PhyloSOLIDvis)

result <- plot_circos(
  inputpath = "path/to/workdir/sampleid/03_tree_building/mutation_integrator/phylo/",
  outputpath = "path/to/circos_plot/",
  annotation_file = "path/to/annotations.txt"
)
```

### Locating the correct input directory

After running PhyloSOLID with `--workdir ./results` and `--sample SAMPLE_ID`, the required files are located in:

```
workdir/
└── SAMPLE_ID/
    └── 03_tree_building/
        └── mutation_integrator/
            └── phylo/
                ├── final_cleaned_I_full_withNA3_for_circosPlot.txt
                ├── final_cleaned_M_full_basedPivots.filtered_sites_inferred.CFMatrix
                ├── df_flipping_count_for_each_mut.txt
                └── df_total_flipping_count.txt
```

**Set `inputpath` to this `phylo/` directory.**

### Using built-in demo data

```r
# Locate demo data
demo_path <- system.file("examples/input", package = "PhyloSOLIDvis")
list.files(demo_path)

# Run with demo data (requires annotation file)
result <- plot_circos(
  inputpath = demo_path,
  outputpath = "output/",
  annotation_file = "path/to/annotation.txt"
)
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `inputpath` | (required) | Path to `phylo/` directory containing PhyloSOLID output files |
| `outputpath` | (required) | Path to output directory for saving plots |
| `annotation_file` | (required) | Path to cell annotation file (TSV format) |
| `target_mut` | "no" | Target mutation ID to highlight |
| `selected_mutlist` | "all" | Comma-separated list of mutations to include |
| `manual_fp_file` | "no" | Path to manual false positive annotation file |
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

All required files are located in the `phylo/` directory after running PhyloSOLID:

| File | Description |
|------|-------------|
| `final_cleaned_I_full_withNA3_for_circosPlot.txt` | Input genotype matrix |
| `final_cleaned_M_full_basedPivots.filtered_sites_inferred.CFMatrix` | Conflict-free matrix |
| `df_flipping_count_for_each_mut.txt` | Per-mutation flipping statistics |
| `df_total_flipping_count.txt` | Total flipping statistics |

### Annotation File Format

The annotation file (TSV format) should contain a `barcode` column and can include:

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
| `No_target.circle_tree_output_as_point.*.svg/pdf` | Main circular plot |
| `legend_components.circos_annotation.svg/pdf` | Circos annotation legend |
| `legend_components.total_flipping_count.svg/pdf` | Flipping count legend |
| `sorted_cf_matrix.txt` | Sorted mutation matrix |
| `ordered_metadata_for_heatmap.txt` | Heatmap ordering data |
| `CNVtree_data_clone_order_tree.rds` | Tree structure data |
| `df_muts_corresponding_to_ordered_tiplabels_by_anticlockwise.txt` | Mutation-cell mapping |

## Example

```r
library(ggplot2)
library(PhyloSOLIDvis)

# Basic usage with PhyloSOLID output
result <- plot_circos(
  inputpath = "results/Org4S15D63/03_tree_building/mutation_integrator/phylo/",
  outputpath = "figures/circos/",
  annotation_file = "data/annotations.txt",
  verbose = TRUE
)

# With target mutation highlighting
result <- plot_circos(
  inputpath = "results/Org4S15D63/03_tree_building/mutation_integrator/phylo/",
  outputpath = "figures/circos/",
  annotation_file = "data/annotations.txt",
  target_mut = "chr11_65426524_T_C",
  tip_label_offset = 8,
  tip_label_size = 3
)
```

## Troubleshooting

### ggplot2 must be loaded first

Always load `ggplot2` before `PhyloSOLIDvis`:

```r
library(ggplot2)
library(PhyloSOLIDvis)
```

### Network issues during installation

Try using a CRAN mirror:

```r
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
remotes::install_github("TsingYang1112/PhyloSOLIDvis", dependencies = TRUE)
```

### Missing dependencies

Install dependencies manually:

```r
# CRAN packages
install.packages(c("ape", "circlize", "cowplot", "dplyr", "ggplot2", "tidyr"))

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
Email: yangqing@westlake.edu.cn  
ORCID: 0009-0005-7366-6879  
GitHub: [@TsingYang1112](https://github.com/TsingYang1112)

## License

MIT © Westlake University

## Citation

If you use PhyloSOLID in your research, please cite:

1. Yang, Q. et al. PhyloSOLID: Robust phylogeny reconstruction from single-cell data despite inherent error and sparsity. (2026) doi:10.64898/2026.02.04.703905.
