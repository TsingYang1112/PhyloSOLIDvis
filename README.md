
# PhyloSOLIDvis

[![R-CMD-check](https://github.com/yourusername/PhyloSOLIDvis/workflows/R-CMD-check/badge.svg)](https://github.com/yourusername/PhyloSOLIDvis/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Visualization Suite for PhyloSOLID Phylogenetic Trees**

## Overview

PhyloSOLIDvis provides comprehensive visualization tools for phylogenetic trees, featuring:

- ✨ **Circos-style circular plots** with integrated mutation heatmaps
- 🔬 **Cell annotation layers** (types, clusters, samples, scores)
- 🧬 **Genotype flipping visualization** (input vs inferred)
- 🎯 **Target mutation highlighting**
- 📊 **Adaptive parameter estimation**

## Installation

### One-line installation (all dependencies auto-installed)

```r
# From GitHub
devtools::install_github("yourusername/PhyloSOLIDvis")
```

### Or with dependency script

```bash
# Download and run
wget https://raw.githubusercontent.com/yourusername/PhyloSOLIDvis/main/install_deps.R
Rscript install_deps.R
```

## Quick Start

```r
library(PhyloSOLIDvis)

result <- plot_circos(
  inputpath = "path/to/PhyloSOLID/output/",
  outputpath = "path/to/figures/",
  annotation_file = "path/to/annotations.txt",
  target_mut = "chr11_65426524_T_C"
)
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| tip_label_offset | 10 | Offset distance for tip labels |
| tip_label_size | 2.5 | Font size for tip labels |
| tip_point_size | 0.5 | Size of tip points |
| heatmap_width | 0.3 | Width of heatmap track |
| heatmap_circos_offset | 0.04 | Offset for heatmap |
| flipping_point_size | 0.2 | Size of flipping markers |
| plot_height | 12 | Output plot height (inches) |
| plot_width | 18 | Output plot width (inches) |

## Required Input Files

- `final_cleaned_I_full_withNA3_for_circosPlot.txt`
- `final_cleaned_M_full_basedPivots.filtered_sites_inferred.CFMatrix`
- `df_flipping_count_for_each_mut.txt`
- `df_total_flipping_count.txt`

## Output Files

- `target_*.svg/pdf` - Main circular plot
- `legend_components.circos_annotation.svg/pdf` - Legend
- `sorted_cf_matrix.txt` - Sorted mutation matrix
- `clone_and_subclone_table.txt` - Clone assignments

## Acknowledgments

The author thanks **Yonghe Xia** for initial conceptualization of the visualization approach and for providing the [converTree](https://github.com/xiayh17/converTree) package, which is used for tree construction from conflict-free matrices.

## Author

**Qing Yang** - *Westlake University, School of Life Sciences*  
Email: qing.yang@westlake.edu.cn  
ORCID: 0009-0005-7366-6879

## License

MIT © Westlake University

## Citation

```bibtex
@software{Yang_PhyloSOLIDvis_2024,
  author = {Qing Yang},
  title = {PhyloSOLIDvis: Visualization Suite for PhyloSOLID Phylogenetic Trees},
  year = {2024},
  publisher = {GitHub},
  url = {https://github.com/yourusername/PhyloSOLIDvis}
}
```
