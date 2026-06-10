# Example usage of PhyloSOLIDvis package

library(PhyloSOLIDvis)

# Basic usage
result <- plot_circos(
  inputpath = "path/to/PhyloSOLID/output/",
  outputpath = "path/to/figures/",
  annotation_file = "path/to/annotations.txt",
  verbose = TRUE
)

# With target mutation highlighting
result <- plot_circos(
  inputpath = "path/to/PhyloSOLID/output/",
  outputpath = "path/to/figures/",
  annotation_file = "path/to/annotations.txt",
  target_mut = "chr11_65426524_T_C",
  tip_label_offset = 10,
  tip_label_size = 2.5,
  verbose = TRUE
)

# Check results
str(result)