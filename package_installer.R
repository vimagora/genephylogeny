# Install CRAN libraries

install.packages("tidyverse")
install.packages("ape")

# Install Bioconductor libraries

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
BiocManager::install("ggtree")
BiocManager::install("treeio")
