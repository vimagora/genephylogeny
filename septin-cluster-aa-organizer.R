rm(list = ls())
library(tidyverse)
library(ggtree)
library(treeio)
library(ape)
library(phytools)
library(Biostrings)

seq_file <- "OG0000195_clean2_Myo2.fa"
tree_file <- "OG0000195_clean2_Myo2_mafft_80trim.treefile"
data_dir <- "local_data"
phyl_dir <- file.path(data_dir, 'phylogeny_analysis')
tree_dir <- file.path(phyl_dir, 'iqtree_files')
seq_dir <- file.path(phyl_dir, 'seq_files')

tree_path <- file.path(tree_dir, tree_file)
seq_path <- file.path(seq_dir, seq_file)

aa_seqs <- readAAStringSet(seq_path)
tree <- read.tree(tree_path)

paralog_seq <-list(set1 = c("Rhimi59_2-1895217",
                            "Rhimi59_2-1893231",
                            "Rhimi59_2-1823729",
                            "Rhimi59_2-1923850",
                            "Rhimi59_2-1876489"),
                   set2 = c("Phybl2-189122",
                            "Phybl2-154374",
                            "Phybl2-136727",
                            "Phybl2-186967"),
                   set3 = c("Phybl2-152587",
                            "Phybl2-162338"),
                   set4 = c("Phybl2-185673",
                            "Phybl2-156601"))


# Extract the sequences
paralog_seqs <- lapply(paralog_seq, function(x) aa_seqs[x])

for (x in names(paralog_seqs)) {
  sequences <- paralog_seqs[[x]]
  writeXStringSet(sequences,file.path(phyl_dir, 'paralog_seqs', paste0(x,"_clean2_Myo2.fa")))
}

# Manually align and curate for paralogs
remove_seqs <- c("Schpo1-829",
                 "Malgl1-2908",
                 "Malgl1-1035",
                 "Schpo1-3667",
                 "Rhimi59_2-1919851",
                 "Rhimi59_2-1616090",
                 "Rhimi59_2-1893231",
                 "Rhimi59_2-1823729",
                 "Rhimi59_2-1923850",
                 "Rhimi59_2-1876489",
                 "Phybl2-189122",
                 "Phybl2-154374",
                 "Phybl2-136727",
                 "Phybl2-152587",
                 "Phybl2-156601",
                 "Mellp2_3-87088",
                 "Yarli1-66281",
                 "Schpo1-4123")

# Remove the selected paralogs
clean3_aa_seq <- aa_seqs[!(names(aa_seqs) %in% remove_seqs)]

writeXStringSet(clean3_aa_seq,file.path(phyl_dir, 'separated_seqs', "OG0000195_clean2_Myo3.fa"))

# Now to separate all of them
tree <- root(tree, outgroup = "Sacce1-Myo2", resolve.root = TRUE)

# According to https://alexknyshov.github.io/R/page3.html modify root branches proportionally
tree$edge.length[which(!(tree$edge[,1] %in% tree$edge[,2]))] <- sum(tree$edge.length[which(!(tree$edge[,1] %in% tree$edge[,2]))])/2

ggtree(tree) +
  geom_tiplab() +
  geom_nodelab(aes(label = node), size = 2, nudge_x = 0.01, nudge_y = 0.5)

groups <- c("CDC3","CDC10", "CDC11", "CDC12", "SPR3", "SPR28", "SHS1","BHR1plus")

nodes <- c(788, 938, 781, 1096, 1250, 1367, 1356, 1263)

names(nodes) <- groups

# Get the descendants of the target node
desc_nodes <- lapply(nodes, function(x) phytools::getDescendants(tree,x))

# The tips are the nodes with numbers less than or equal to the number of tips (Ntip)
tip_nodes <- lapply(desc_nodes, function(x) x[x <= Ntip(tree)])

# Get the labels for these tip nodes
tip_labels <- lapply(tip_nodes, function(x) tree$tip.label[x])

# Extract the sequences
septin_seqs <- lapply(tip_labels, function(x) aa_seqs[x])

for (group in groups) {
  sequences <- septin_seqs[[group]]
  clean_seqs <- sequences[!(names(sequences) %in% remove_seqs)]
  writeXStringSet(clean_seqs,file.path(phyl_dir, 'separated_seqs', paste0(group,"_clean2_Myo2.fa")))
}







