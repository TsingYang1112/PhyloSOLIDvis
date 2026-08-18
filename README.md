# PhyloSOLIDvis

**Visualization Suite for PhyloSOLID Phylogenetic Trees**

## Overview

PhyloSOLIDvis generates publication-ready circular (circos-style) phylogenetic tree visualizations with:

- ✨ **Circos-style circular plots** with integrated mutation heatmaps
- 🔬 **Cell annotation layers** (types, clusters, samples, tumor scores, B cell proportions)
- 🧬 **Genotype flipping visualization** (input vs inferred states)
- 🎯 **Target mutation highlighting**
- 📊 **Heatmap visualization** for mutation profiles
- 🚀 **Optimized performance** with shared data workflow (2.2x faster)

## Example Output

Below is an example circos plot generated from PhyloSOLID output (351 cells, 17 mutations):

![Circos Plot](inst/demo_circos.png)

## Installation

### Recommended: Conda environment (avoids compilation issues)

For users who encounter compilation errors (especially with `igraph`), we recommend using conda:

```bash
# Create environment with specific R version
conda create -n phylosolid_env r-base=4.5.3 -c conda-forge
conda activate phylosolid_env

# Install core dependencies via conda (pre-compiled)
mamba install -c conda-forge \
  r-ape r-circlize r-ggplot2 r-dplyr r-tidyr r-cowplot \
  r-igraph r-readr r-vroom r-v8

# Install Bioconductor packages
conda install -c bioconda \
  r-biocmanager r-ggtree r-treeio r-complexheatmap

# Start R and install PhyloSOLIDvis
R
```

Then in R:

```r
# Install converTree dependency
remotes::install_github("xiayh17/converTree")

# Install PhyloSOLIDvis
remotes::install_github("TsingYang1112/PhyloSOLIDvis", dependencies = TRUE)
```

### Standard R installation

```r
# CRAN packages
install.packages(c("ape", "circlize", "cowplot", "dplyr", "ggplot2", "tidyr", 
                   "readr", "vroom", "ggridges", "ggtext", "glue"))

# Bioconductor packages
BiocManager::install(c("ggtree", "ggtreeExtra", "treeio", "ComplexHeatmap"))

# GitHub packages
remotes::install_github("xiayh17/converTree")
remotes::install_github("TsingYang1112/PhyloSOLIDvis", dependencies = TRUE)
```

### One-line installation (if dependencies are already satisfied)

```r
remotes::install_github("TsingYang1112/PhyloSOLIDvis", dependencies = TRUE)
```

## Quick Start

**Note:** `ggplot2` must be loaded before `PhyloSOLIDvis` to avoid namespace conflicts.

```r
library(ggplot2)
library(PhyloSOLIDvis)

# Set your paths
# inputpath: Point to the phylo/ directory from your PhyloSOLID run
inputpath <- "path/to/workdir/sampleid/03_tree_building/mutation_integrator/phylo/"
outputpath <- "path/to/your/output/directory/"
```

### Option 1: One-click complete pipeline (simplest, recommended)

```r
results <- run_all(
  inputpath = inputpath,
  outputpath = outputpath,
  annotation_file = "path/to/annotations.txt",  # optional
  verbose = TRUE
)
```

### Option 2: Step-by-step with shared data (more flexible, 2.2x faster)

```r
# Step 1: Prepare data once (reads files, sorts matrix, builds tree)
data <- prepare_data(
  inputpath = inputpath,
  annotation_file = "path/to/annotations.txt",  # optional
  verbose = TRUE
)

# Step 2: Generate circos plot
circos <- plot_circos(
  data = data,
  outputpath = outputpath,
  tip_label_offset = 10,        # Adjust label distance from tree
  flipping_point_size = 0.2,    # Adjust size of genotype flipping markers
  heatmap_width = 0.3,          # Adjust heatmap track width
  verbose = TRUE
)

# Step 3: Generate heatmap (reuses tree from circos)
heatmap <- plot_heatmap(
  data = data,
  tree = circos$clone_order_tree,  # Reuse tree!
  outputpath = outputpath,
  verbose = TRUE
)
```

**Choose Option 1** if you just want the default plots with minimal code.  
**Choose Option 2** if you want to customize parameters or reuse data for multiple plots.

### Performance Comparison

| Method | File I/O | Tree Building | Speedup |
|--------|----------|---------------|---------|
| Traditional (separate steps) | Multiple times | Multiple times | 1x |
| **Option 1 (run_all)** | **Once** | **Once** | **2.2x** |
| **Option 2 (shared data)** | **Once** | **Once** | **2.2x** |

## Key Parameters for Fine-Tuning

The three most important parameters for adjusting your circos plot:

| Parameter | Default | Description | When to Adjust |
|-----------|---------|-------------|----------------|
| **`tip_label_offset`** | 10 | Distance of tip labels from the tree | Labels overlapping → Increase; Too far → Decrease |
| **`flipping_point_size`** | 0.2 | Size of genotype flipping markers | Too small → Increase; Too large → Decrease |
| **`heatmap_width`** | 0.3 | Width of the heatmap track | Too narrow → Increase; Too wide → Decrease |

### Quick Adjustment Guide

```r
# Default look
circos_default <- plot_circos(
  data = data,
  outputpath = outputpath,
  tip_label_offset = 10,
  flipping_point_size = 0.2,
  heatmap_width = 0.3
)

# Cleaner look with labels closer and markers larger
circos_clean <- plot_circos(
  data = data,
  outputpath = outputpath,
  tip_label_offset = 6,          # Bring labels closer
  flipping_point_size = 1.3,     # Make markers more visible
  heatmap_width = 0.4            # Widen heatmap
)

# For dense trees with many cells
circos_dense <- plot_circos(
  data = data,
  outputpath = outputpath,
  tip_label_offset = 12,         # Push labels out to avoid overlap
  flipping_point_size = 0.1,     # Smaller markers to reduce clutter
  heatmap_width = 0.25           # Narrower heatmap for more tree space
)
```

## Locating the Correct Input Directory

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
                ├── df_total_flipping_count.txt
                └── df_barcode_clones_from_phylo_tree.csv
```

**Set `inputpath` to this `phylo/` directory.**

## Functions

| Function | Description | Recommended For |
|----------|-------------|-----------------|
| `run_all()` | Complete pipeline (data prep + circos + heatmap) | **Most users** |
| `prepare_data()` | Prepare shared data object (read once) | Efficient step-by-step workflow |
| `plot_circos()` | Generate circular phylogenetic tree plot | Custom circos plots |
| `plot_heatmap()` | Generate heatmap visualization | Custom heatmaps |

## Parameters Reference

### prepare_data Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `inputpath` | (required) | Path to `phylo/` directory |
| `annotation_file` | NULL | Path to annotation file (TSV, optional) |
| `verbose` | TRUE | Print progress messages |

**Returns:** A `PhyloData` object containing all prepared data

### plot_circos Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `inputpath` | NULL | Path to `phylo/` directory (optional if `data` provided) |
| `outputpath` | (required) | Path to output directory |
| `data` | NULL | PhyloData object from `prepare_data()` |
| `annotation_file` | NULL | Path to annotation file |
| `target_mut` | "no" | Target mutation to highlight |
| `selected_mutlist` | "all" | Comma-separated list of mutations |
| `manual_fp_file` | "no" | Manual false positive file |
| **`tip_label_offset`** | 10 | **Distance of tip labels** |
| `tip_label_size` | 2.5 | Font size for tip labels |
| `tip_point_size` | 0.5 | Size of tip points |
| **`heatmap_width`** | 0.3 | **Width of heatmap track** |
| `heatmap_circos_offset` | 0.04 | Offset for heatmap from tree |
| **`flipping_point_size`** | 0.2 | **Size of flipping markers** |
| `plot_height` | 12 | Output plot height (inches) |
| `plot_width` | 18 | Output plot width (inches) |
| `verbose` | TRUE | Print progress messages |

**Returns:** A list containing `main_plot`, `svg_file`, `pdf_file`, `clone_order_tree`, `leaf_df`, `subclones`

### plot_heatmap Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `inputpath` | NULL | Path to `phylo/` directory (optional if `data` provided) |
| `outputpath` | (required) | Path to output directory |
| `data` | NULL | PhyloData object from `prepare_data()` |
| `tree` | NULL | phylo object (reuse from circos for best performance) |
| `ordered_metadata_file` | NULL | Path to ordered_metadata_for_heatmap.txt |
| `colors` | c("#D4E8F0", "#7D2224") | Colors for 0 and 1 values |
| `na_color` | "white" | Color for NA/missing values |
| `show_colnames` | TRUE | Show column names (mutation IDs) |
| `show_rownames` | FALSE | Show row names (cell IDs) |
| `cluster_cols` | FALSE | Cluster columns (mutations) |
| `border_color` | NA | Border color for heatmap cells |
| `verbose` | TRUE | Print progress messages |

**Returns:** A list containing `tree`, `matrix`, `tip_order`, `output_files`

### run_all Parameters

`run_all()` accepts all parameters from `plot_circos()` and `plot_heatmap()` plus:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `run_adjusted` | FALSE | Run adjusted plot with different parameters |
| `adjusted_tip_label_offset` | 6 | Tip label offset for adjusted plot |
| `adjusted_flipping_point_size` | 1.3 | Flipping point size for adjusted plot |
| `heatmap_colors` | c("#D4E8F0", "#7D2224") | Colors for heatmap |
| `heatmap_na_color` | "white" | Color for NA values in heatmap |

## Required Input Files

All required files are located in the `phylo/` directory after running PhyloSOLID:

| File | Description |
|------|-------------|
| `final_cleaned_I_full_withNA3_for_circosPlot.txt` | Input genotype matrix |
| `final_cleaned_M_full_basedPivots.filtered_sites_inferred.CFMatrix` | Conflict-free matrix |
| `df_flipping_count_for_each_mut.txt` | Per-mutation flipping statistics |
| `df_total_flipping_count.txt` | Total flipping statistics |
| `df_barcode_clones_from_phylo_tree.csv` | Clone assignment and color information |

### Annotation File Format (Optional)

The annotation file (TSV format) should contain a `barcode` column and can include:

| Column | Description |
|--------|-------------|
| `cluster_info` | Cell cluster assignments |
| `cell_type` | Cell type annotations |
| `sample` | Sample identifiers |
| `tumor_score` | Numeric tumor scores |
| `B_cell_prop` | B cell proportions |

**Note:** `annotation_file` is optional. If not provided, the plot will be generated without annotation layers (tree + heatmap only).

## Output Files

After running, your output directory contains:

| File | Description |
|------|-------------|
| `No_target.circle_tree_output_as_point.*.svg/pdf` | Main circular plot |
| `heatmap_and_histograms_for_our_tree.pdf` | Heatmap |
| `legend_components.circos_annotation.svg/pdf` | Circos annotation legend |
| `legend_components.total_flipping_count.svg/pdf` | Flipping count legend |
| `sorted_cf_matrix.txt` | Sorted mutation matrix |
| `ordered_metadata_for_heatmap.txt` | Heatmap ordering data |
| `CNVtree_data_clone_order_tree.rds` | Tree structure data |
| `df_muts_corresponding_to_ordered_tiplabels_by_anticlockwise.txt` | Mutation-cell mapping |
| `phylo_tree.rds` | Phylogenetic tree object |

## Interactive Online Inspection Platform

For interactive post-analysis and quality control, upload your PhyloSOLIDvis outputs to:

**Web Server**: [https://phylosolid.westlake.edu.cn](https://phylosolid.westlake.edu.cn)

> **🔐 Access:** Register with your **institutional email address** (e.g., `@*.edu`, `@*.ac.cn`, `@*.edu.cn`). Registration is free for academic users.

**Features:**
- Interactive Circos exploration with zoom and cell-level query
- Tree topology evaluation and quality control
- Genotype flip statistics for assessing tree reliability
- Optional upload of IGV snapshots and ANNOVAR annotation results

### Required Upload Files

| Upload Field | Required File | Description |
|:-------------|:--------------|:------------|
| **Binary matrix for heatmap layer** | `ordered_metadata_for_heatmap.txt` | Ordered cell IDs with metadata |
| **Ordered cell IDs with top mutations** | `df_muts_corresponding_to_ordered_tiplabels_by_anticlockwise.txt` | Mutation order |
| **Circos plot (SVG)** | `No_target.circle_tree_output_as_point.*.svg` | Main circular tree plot |
| **Legend figure** | `legend_components.circos_annotation.svg` | Annotation legend |
| **Global flip statistics** | `df_total_flipping_count.txt` | Total flipping statistics |
| **Per-mutation flip statistics** | `df_flipping_count_for_each_mut.txt` | Per-mutation statistics |

### Where to Find These Files

```
your_output_directory/
├── ordered_metadata_for_heatmap.txt
├── df_muts_corresponding_to_ordered_tiplabels_by_anticlockwise.txt
├── No_target.circle_tree_output_as_point.*.svg
├── legend_components.circos_annotation.svg
├── df_total_flipping_count.txt          (from phylo/ directory)
└── df_flipping_count_for_each_mut.txt   (from phylo/ directory)
```

### Deep Inspection Features

At the bottom of the webpage, you can upload supplementary files:

- **ANNOVAR variant annotation results**: Gene context, genomic regions, functional predictions
- **IGV snapshot figures** (PNG): Raw sequencing reads at each mutation site
  - Each PNG named with exact mutation ID (e.g., `chr1_39034563_T_A.png`)

## Run Demo

After installation, test with built-in demo data:

```r
library(ggplot2)
library(PhyloSOLIDvis)

demo_path <- system.file("examples/input", package = "PhyloSOLIDvis")
output_path <- "demo_output/"
annotation_file <- system.file("examples/annotation.txt", package = "PhyloSOLIDvis")

# Method 1: One-click
results <- run_all(
  inputpath = demo_path,
  outputpath = output_path,
  annotation_file = annotation_file,
  verbose = TRUE
)

# Method 2: Step-by-step
data <- prepare_data(demo_path, annotation_file, verbose = TRUE)
circos <- plot_circos(
  data = data,
  outputpath = output_path,
  tip_label_offset = 6,
  flipping_point_size = 1.3,
  verbose = TRUE
)
heatmap <- plot_heatmap(
  data = data,
  tree = circos$clone_order_tree,
  outputpath = output_path,
  verbose = TRUE
)

# View generated files
list.files(output_path)
```

Or run the demo script:

```bash
wget https://raw.githubusercontent.com/TsingYang1112/PhyloSOLIDvis/main/run_demo.R
Rscript run_demo.R
```

## Troubleshooting

### ggplot2 must be loaded first

```r
library(ggplot2)  # Always load first
library(PhyloSOLIDvis)
```

### Common mistake: Passing data object to run_all()

```r
# ❌ Wrong: run_all() expects a path string
data <- prepare_data(inputpath)
results <- run_all(inputpath = data, outputpath = outputpath)  # Error!

# ✅ Correct: run_all() expects a path
results <- run_all(inputpath = inputpath, outputpath = outputpath)

# ✅ Or use data object with plot functions directly
circos <- plot_circos(data = data, outputpath = outputpath)
heatmap <- plot_heatmap(data = data, tree = circos$clone_order_tree, outputpath = outputpath)
```

### Network issues during installation

```r
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
remotes::install_github("TsingYang1112/PhyloSOLIDvis", dependencies = TRUE)
```

### Missing dependencies

```r
# CRAN packages
install.packages(c("ape", "circlize", "cowplot", "dplyr", "ggplot2", "tidyr"))

# Bioconductor packages
BiocManager::install(c("ggtree", "ggtreeExtra", "treeio", "ComplexHeatmap"))

# GitHub packages
remotes::install_github("xiayh17/converTree")
```

### R version mismatch issues

If you encounter `undefined symbol: R_getRegisteredNamespace` errors, this indicates packages were compiled under a different R version. Reinstall the problematic packages from source:

```r
install.packages(c("readr", "vroom"), type = "source")
remotes::install_github("TsingYang1112/PhyloSOLIDvis")
```

## Version History

### v1.1.0 (2026-08-07)
- 🚀 **Major performance optimization**: Added `prepare_data()` and `PhyloData` class
- ⚡ **2.2x faster**: Data read once, tree built once across all functions
- 🔄 **Backward compatible**: All existing functions still work
- 📦 **Shared data workflow**: Efficient step-by-step analysis
- 🐛 Fixed swapped legend labels (False positive/False negative)
- 📝 Improved documentation with key parameter guidance

### v1.0.0
- Initial release
- Core functions: `run_all()`, `plot_circos()`, `plot_heatmap()`, `order_matrix()`
- Support for annotation layers and target mutation highlighting

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
