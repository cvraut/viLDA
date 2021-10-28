# note it takes 8GB RAM and 30 min to install celda
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install("celda")

library(celda)
set.seed("19890418")

simsce <- simulateCells("celda_CG", S = 5, K = 5, L = 10, G = 200, CRange = c(30, 50))

library(SingleCellExperiment)
dim(simsce$counts)
# 207 cells in 5 groups with 200 genes
# generative model is the true LDA model



