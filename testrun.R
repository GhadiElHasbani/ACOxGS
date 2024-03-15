
#######################
#     Dependencies
#######################

library(dplyr)
library(ggplot2)
library(purrr)
require(pcaPP) # GS dependency
library(igraph) # GS dependency
library(limma) # DEA tool / GS dependency
library(GSAR) # p53 data
library(GSEABenchmarkeR) # GSEABenchmark data
library(org.Hs.eg.db)
library(annotate) #map NCBI Gene IDs to gene symbols
library(LEANR) # LEAN package
library(doMC) # LEAN dependency
# For KEGG network
library(BiocManager)
library(org.Hs.eg.db)
library(KEGGgraph)
library(XML)
library(graphite)
library(CHRONOS)

setwd("C:/Users/Ghadi/Desktop/GS_optimization")
setwd("~/Desktop/GS_optimization/")
setwd("~/Desktop/GeneSurrounder_optimization/GS_optimization")

#setwd("....")

source("GeneSurrounder/GeneSurrounder.R")
source("ant_colony.R")

#######################
#   Getting data 
#######################

###### GSEA data

  # TP53 
data("p53DataSet")

df <- p53DataSet
class_labels <- c(rep("WT", 17), rep("MUT", 33))


###### GSEAbenchmarkR data

geo2kegg <- loadEData("geo2kegg", preproc = TRUE)
#geo2kegg <- maPreproc(geo2kegg)

ids <- c()
sizes <- c()
for(id in names(geo2kegg)){
  ids<- c(ids, id)
  idx <- which(names(geo2kegg) == id)
  eset <- geo2kegg[[idx]]
  sizes <- c(sizes, dim(eset)[2])
}
bm_ds_ids <- data.frame(ids, sizes) %>% dplyr::arrange(desc(sizes))

  # GSE20291: two groups of age and gender matched groups of 14 Parkinson and 19 Control subjects
idx <- which(names(geo2kegg) == "GSE20291")

  # GSE5281_VCX: 19 Alzheimer's and 12 normal-aged brain samples
idx <- which(names(geo2kegg) == "GSE5281_VCX")

  # GSE8762: 12 moderate stage Huntington's patients (8 female and 4 male) and 10 age-matched controls (5 female and 5 male)
idx <- which(names(geo2kegg) == "GSE8762")

  # GSE5281_HIP: Alzheimer's
idx <- which(names(geo2kegg) == "GSE5281_HIP")

  # GSE5281_EC: Alzheimer's
idx <- which(names(geo2kegg) == "GSE5281_EC")

  # General
eset <- geo2kegg[[idx]]
df <- eset@assays@data@listData$exprs
class_labels <- as.character(eset@colData@listData$GROUP)
gene_symbols <- unlist(lookUp(eset@NAMES, 'org.Hs.eg', 'SYMBOL'))
sum(rownames(df) != names(gene_symbols))
rownames(df)<- gene_symbols

  ### for cross-study concordance
idx1 <- which(names(geo2kegg) == "GSE5281_VCX")
eset1 <- geo2kegg[[idx1]]
idx2 <- which(names(geo2kegg) == "GSE5281_HIP")
eset2 <- geo2kegg[[idx2]]
idx3 <- which(names(geo2kegg) == "GSE5281_EC")
eset3 <- geo2kegg[[idx3]]

df1 <- eset1@assays@data@listData$exprs
df2 <- eset2@assays@data@listData$exprs
df3 <- eset3@assays@data@listData$exprs
common_genes <- rownames(df1)[rownames(df1) %in% rownames(df2)][rownames(df1)[rownames(df1) %in% rownames(df2)] %in% rownames(df3)]

eset <- geo2kegg[[idx]]
df <- eset@assays@data@listData$exprs[common_genes,]
class_labels <- as.character(eset@colData@listData$GROUP)
gene_symbols <- unlist(lookUp(eset@NAMES, 'org.Hs.eg', 'SYMBOL'))
sum(rownames(df) != names(gene_symbols))
rownames(df)<- gene_symbols

##########################################################
#     Building KEGG Network (DO NOT RUN) 
##########################################################

EG_IDs <- mget(rownames(df), revmap(org.Hs.egSYMBOL),ifnotfound=NA)
KEGG_IDs <- mget(as.character(EG_IDs), org.Hs.egPATH,ifnotfound=NA)

x <- org.Hs.egPATH
# Get the entrez gene identifiers that are mapped to a KEGG pathway ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
pathway_ids <- unique(unlist(xx))[-c(which(unique(unlist(xx)) %in% c("00072", "00300", "00472", "04320", 
                                                                     "00460", "00471"#, "01100"
                                                                     )))] 
                                                                      #Synthesis and degradation of ketone bodies pathway (hsa00072) not found for homo sapiens
                                                                      #Lysine Biosynthesis (hsa00300) not found for homo sapiens
                                                                      #D-Arginine and D-ornithine metabolism (hsa00472) not found for homo sapiens
                                                                      #Dorso-ventral axis formation (hsa04320) not found for homo sapiens
                                                                      #Cyanoamino acid metabolism (hsa00460) not found for homo sapiens
                                                                      #D-Glutamine and D-glutamate metabolism (hsa00471) not found for homo sapiens
graphs <- c()
options(timeout=300)
for(i in 1:length(pathway_ids)){
  pathway_id <- pathway_ids[i]
  tmp <- paste("hsa", pathway_id, ".xml", sep = "")
  
  if(!file.exists(tmp)){
    kgmls <- retrieveKGML(pathway_id, destfile = tmp, organism = "hsa")
    #sfiles <- system.file(paste("hsa", pathway_ids[1], ".xml", sep = ""), package = "KEGGgraph")
    #sfiles <- parseKGML(sfiles)
  }
  graphs <- c(graphs, parseKGML2Graph(tmp))
}
#paths <- downloadPathways(org='hsa', pathways=pathway_ids[1])
#graphs <- createPathwayGraphs(org='hsa', pathways=paths)
KEGG_graph <- mergeGraphs(as.list(graphs), edgemode = "directed")
KEGG_A <- as(KEGG_graph, "matrix")

x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
rownames(KEGG_A) <- xx[translateKEGGID2GeneID(rownames(KEGG_A), "hsa")]
colnames(KEGG_A) <- xx[translateKEGGID2GeneID(colnames(KEGG_A), "hsa")]

sum(rowSums(KEGG_A) == 0) # -> 2683 (p53), 2681 (GSE20291), 2681 (GSE5281_VCX), 2681 (GSE8762)
#saveRDS(KEGG_A, "Networks/KEGG_A_.....rds")

#Make graph undirected
for(i in 1:nrow(KEGG_A)){
  for(j in 1:ncol(KEGG_A)){
    if(i > j){
      if(KEGG_A[i, j] != KEGG_A[j, i]){
        KEGG_A[i, j] = max(c(KEGG_A[i, j], KEGG_A[j, i]))
        KEGG_A[j, i] = max(c(KEGG_A[i, j], KEGG_A[j, i]))
      }
    }
  }
}
isSymmetric(KEGG_A) # -> TRUE
sum(rowSums(KEGG_A) == 0) # -> 1832 (p53), 1832 (GSE20291), 1832 (GSE5281_VCX), 1832 (GSE8762)
#saveRDS(KEGG_A, "Networks/KEGG_A_undir_.....rds")

##########################
# Reading RDS file to CSV
##########################

am_rds <- readRDS("Networks/KEGG_A_undir_......rds")
write.csv(am_rds, "Networks/KEGG_A_undir_......csv", row.names=TRUE)

am_rds <- readRDS("NetworksKEGG_A_undir_p53.rds")
write.csv(am_rds, "Networks/KEGG_A_undir_p53.csv", row.names=TRUE)

am_rds <- readRDS("Networks/KEGG_A_undir_GE20291.rds")
write.csv(am_rds, "Networks/KEGG_A_undir_GE20291.csv", row.names=TRUE)

#############################################
# perform networkX LCC extraction --> python
#############################################

#############################################
# Importing LCC CSV & checking AM correctness
#############################################

#lcc<- as.matrix(read.csv('./Networks/KEGG_A_undir_....._LCC.csv', row.names = 1))
lcc<- as.matrix(read.csv('./Networks/KEGG_A_undir_GE20291_LCC.csv', row.names = 1))
#View(lcc)
ncol(lcc)
nrow(lcc)
sum(rownames(lcc) != colnames(lcc))
sum(!rownames(lcc) %in% colnames(lcc))
sum(!colnames(lcc) %in% rownames(lcc))

rownames(lcc)[!rownames(lcc) %in% colnames(lcc)] <- sub('-', '.', rownames(lcc)[!rownames(lcc) %in% colnames(lcc)])
sum(rownames(lcc) != colnames(lcc))

#########################
#     Parallelized ACO
#########################

library(parallel)
source("ant_colony.R")
source("OSCC.R")

max_cores <- detectCores()

seed <- 7
dataset <- "20291"
KEGG_A_LCC <- lcc                                     
aco_modules <- optimize_modules(df, KEGG_A_LCC, n_ants = 40, n_iter = 5, starting_capacity = 1.0, 
                                n_cores = min(5, max_cores - 2), alpha = 0.6, beta = 1.2, 
                                seed = seed, n_resamples = 100, class_labels = class_labels, penalize_dist = TRUE)

saveRDS(aco_modules, paste('Objects/Run ', seed, paste("/aco_gamma_trial", dataset, "40", seed, sep="_"), '.rds', sep = ''))
saveRDS(aco_modules, 'aco_gamma_trial_20291_40_7.rds')
#gamma: n_iter = 5, starting_capacity = 1.0, n_cores = min(5, max_cores - 2), alpha = 0.6, beta = 1.2, n_resamples = 100, penalize_dist = TRUE

#########################
#     GS
#########################

source("GeneSurrounder/GeneSurrounder.R")
source("GeneSurrounder/run_geneSurrounder.R")
seed <- 127
gs_modules <- run_geneSurrounder(adj.matrix = KEGG_A_LCC,
                                 gene_scores_per_sample = df,
                                 class_labels = class_labels,
                                 nresamples = 100,
                                 seed = seed,
                                 num.Sphere.resamples = 100,
                                 gene.id = names(aco_modules),
                                 decay_only = FALSE,
                                 #file_name = "gs_results_{.....}_{seed - 120}.csv", #<- fill before running
                                 cores = 1 # Set to 1 for Windows
)
saveRDS(gs_modules, paste('Objects/Run ', seed - 120, paste('/gs_results',dataset, seed - 120, sep="_"), '.rds', sep = ''))
saveRDS(gs_modules, 'gs_results_20291_7.rds')
#########################
#     Module Graphs
#########################

library(visNetwork)

#random module generator
plot_rand_module <- function(aco_modules, seed = 126, gene = "none"){
  set.seed(seed)
  if(gene == "none"){
    rand_gene <- sample(names(aco_modules), 1)
  }else{
    rand_gene <- gene
  }
  p_decay_min <- min(aco_modules[[rand_gene]]$stats[,'p_fishers'])
  print(paste("Module center: ", rand_gene))
  print(paste("Module p-Fisher: ", p_decay_min))
  
  set.seed(seed)
  low_id <- sample(which(aco_modules[[rand_gene]]$stats[,'p_fishers'] == p_decay_min), 1)
  module <- aco_modules[[rand_gene]]$modules$modules[[low_id]]
  A_subset <- KEGG_A_LCC[rownames(KEGG_A_LCC) %in% module, colnames(KEGG_A_LCC) %in% module]
  
  if(length(module) > 1){
    mod_graph <- graph_from_adjacency_matrix(A_subset) %>% 
      set_vertex_attr("class", value = append(rep('neighbor', length(module) - 1), 'center', after = which(rownames(A_subset) == rand_gene) - 1))
    vis_mod <- toVisNetworkData(mod_graph)
    visNetwork(nodes = vis_mod$nodes, edges = vis_mod$edges) %>%
      visIgraphLayout(layout = "layout_nicely") %>%
      visOptions(nodesIdSelection = TRUE, selectedBy = "class") %>%
      visNodes(color = vis_mod$nodes$class)
  } else{
    mod_graph <- graph_from_adjacency_matrix(A_subset) %>% 
      set_vertex_attr("class", value = c('center'))
    vis_mod <- toVisNetworkData(mod_graph)
    visNetwork(nodes = vis_mod$nodes, edges = vis_mod$edges) %>%
      visOptions(nodesIdSelection = TRUE, selectedBy = "class") %>%
      visNodes(color = vis_mod$nodes$class)
  }
  
  
}

plot_rand_module(aco_modules, seed = 113, gene = "PIK3CG")

#Genes where module is gene itself only:
#ALG14
#BRIP1
#CDC7

#CXCL14

#########################
#     LEAN
#########################
######### NO SEED OPTION --> not reproducible

runLEAN <- function(gene_scores_per_sample,
                    adj_mat,
                    n_cores = 1,
                    seed,
                    n_resamples = 1000,
                    class_labels){
  require(limma)
  require(LEANR)
  require(doMC)
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  print("Preprocessing...")
  gene_scores_per_sample <- gene_scores_per_sample[rownames(gene_scores_per_sample) %in% rownames(adj_mat),]
  adj_mat <- adj_mat[rownames(adj_mat) %in% rownames(gene_scores_per_sample), ]
  adj_mat <- adj_mat[ ,colnames(adj_mat) %in% rownames(gene_scores_per_sample)]
  int_network <- graph_from_adjacency_matrix(adj_mat, mode = "undirected")
  
  # Calc gene level statistics & a null set of gene level stats (shuffle phenotype labels)
  calcGeneTStats <- function(expr,
                             classLabels){
    # Calc gene level statistics
    # This code is being written using CurOvGradeKEGGnets[[2]]
    #
    # Args:
    #    expr: is a matrix of genes by samples
    #    classLabels: is a vector of class labels (e.g. high vs low)
    # Returns:
    #    observedStats: a vector of observed moderated t-statistics
    
    # cf. d715_timecourse_contrasts, network_review_GSEAhat ?
    # Should I save the fit so I have the gene p values etc...?
    
    require(limma)
    print("Calculating observed gene stats using limma...")
    desMat =  model.matrix(~factor(classLabels))
    # treat is a limma function
    # Given a microarray linear model fit, compute moderated t-statistic, etc
    fit =  treat(lmFit(expr,desMat))
    observedStats = as.data.frame(fit$p.value[,2])
    
    res <- observedStats %>% pull()
    names(res) <- rownames(observedStats)
    return(res)
    
  }

  gene_stats <- calcGeneTStats(gene_scores_per_sample, class_labels)
  print("Running LEAN method...")
  
  return(LEANR::run.lean(gene_stats, int_network, verbose = TRUE, ranked = FALSE, n_reps = n_resamples, ncores = n_cores))
}

KEGG_A_LCC <- lcc    
seed <- 3
lean_modules <- runLEAN(df, KEGG_A_LCC, n_cores = NULL, seed = seed, 
                        n_resamples = 100, class_labels = class_labels)

#saveRDS(lean_modules, 'lean_modules_.....rds')
saveRDS(lean_modules, paste('Objects/Run ', seed, '/lean_modules_', seed, '.rds', sep = ''))

