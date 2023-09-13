options(timeout = 10000000000000000000000000)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GSAR")
BiocManager::install("GSEABenchmarkeR")
BiocManager::install("KEGGgraph")
BiocManager::install("CHRONOS")
BiocManager::install("graphite")
BiocManager::install("limma")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("annotate")

require(dplyr)
require(ggplot2)
require(purrr)
require(pcaPP) # GS dependency
require(igraph) # GS dependency
require(limma) # DEA tool / GS dependency
require(LEANR) # LEAN package
require(doMC) # LEAN dependency
require(GSAR) # p53 data
# For KEGG network
require(org.Hs.eg.db)
require(KEGGgraph)
require(XML)
require(CHRONOS)
require(graphite)