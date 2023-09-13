options(timeout = 10000000000000000000000000)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("curatedBladderData")
BiocManager::install("curatedOvarianData")
BiocManager::install("GSAR")
BiocManager::install("KEGGgraph")
BiocManager::install("CHRONOS")
BiocManager::install("graphite")
BiocManager::install("DOSE")
BiocManager::install("limma")
BiocManager::install("org.Hs.eg.db")

require(mND)
require(dplyr)
require(ggplot2)
require(purrr)
require(pcaPP) # GS dependency
require(igraph) # GS dependency
require(limma) # GS dependency
require(curatedBladderData) # GS data
require(curatedOvarianData) # GS data
require(GSAR) # p53 data
require(DOSE) # Enrichment tool
# For KEGG network
require(org.Hs.eg.db)
require(KEGGgraph)
require(XML)
require(CHRONOS)
require(graphite)

