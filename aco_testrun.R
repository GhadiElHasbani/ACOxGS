
#######################
#     Dependencies
#######################

#options(timeout = 10000000)
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("curatedBladderData")
#BiocManager::install("curatedOvarianData")
#BiocManager::install("GSAR")
#BiocManager::install("GSEABenchmarkeR")
#BiocManager::install("KEGGgraph")
#BiocManager::install("CHRONOS")
#BiocManager::install("graphite")
#BiocManager::install("DOSE")
#BiocManager::install("limma")

#library(devtools)
#install_github("GSEA-MSigDB/GSEA_R")

source(system.file('extdata', 'Run.GSEA.R', package = 'GSEA'))
#library(mND)
library(dplyr)
library(ggplot2)
library(purrr)
library(preprocessCore) # preprocessing (QN)
require(pcaPP) # GS dependency
library(igraph) # GS dependency
library(limma) # DE tool / GS dependency
library(curatedBladderData) # GS data
library(curatedOvarianData) # GS data
library(GSAR) # p53 data
library(GSEABenchmarkeR) # GSEABenchmark data
library(org.Hs.eg.db)
library(annotate) #map NCBI Gene IDs to gene symbols
library(LEANR) # LEAN package
library(doMC) # LEAN dependency
library(DOSE) # Enrichment tool
library(mlr3) # Importing benchmark data from OSF
library(readr) # Importing benchmark data from OSF
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
#setwd("C:/Users/Frederick/Desktop/GitHub_Repos/GS_optimization")
source("old_genesurrounder/GeneSurrounder.R")
source("ant_colony.R")

#######################
#   ACO on mND data 
#######################

data(A)
int_network <- graph_from_adjacency_matrix(A, mode = "undirected")
distance <- calcAllPairsDistances(int_network, directionPaths = "all", networkName = "int_network")
#norm_distance <- normalize_adj_mat(distance)
normalized_distance <- t(as.data.frame(t(distance)) %>% mutate_all(~ ./max(.)))

data(X0)
gene_scores <- data.frame(X0) %>% select(L2)

aco_modules_single <- optimize_module(gene_id = rownames(distance)[1],
                                      gene_scores = gene_scores,
                                      adj_mat = A,
                                      distance = distance,
                                      normalized_distance = normalized_distance,
                                      alpha = 0.6,
                                      beta = 1.2,
                                      n_iter = 3,
                                      n_ants = 10,
                                      seed = 3)


#######################
#   Getting data 
#######################

## Bladder

#data(package="curatedBladderData")
data(GSE13507_eset)

eset <- GSE13507_eset
df <- exprs(eset)
class_labels <- phenoData(eset)@data$summarystage

data(GSE31684_eset)

eset <- GSE31684_eset
df <- exprs(eset)
class_labels <- phenoData(eset)@data$summarystage

data(GSE32894_eset)

eset <- GSE32894_eset
df <- exprs(eset)
class_labels <- phenoData(eset)@data$summarystage

## Ovarian

  ## version 1
data(GSE14764_eset)
data(GSE17260_eset)
data(GSE9891_eset)
df1 <- exprs(GSE14764_eset)
df2 <- exprs(GSE17260_eset)
df3 <- exprs(GSE9891_eset)
common_genes <- rownames(df1)[rownames(df1) %in% rownames(df2)][rownames(df1)[rownames(df1) %in% rownames(df2)] %in% rownames(df3)]

eset <- GSE14764_eset
df <- exprs(eset)[common_genes,]
class_labels <- phenoData(eset)@data$summarygrade

eset <- GSE17260_eset
df <- exprs(eset)[common_genes,]
class_labels <- phenoData(eset)@data$summarygrade

eset <- GSE9891_eset
df <- exprs(eset)
class_labels <- phenoData(eset)@data$summarygrade
na_filter <- !is.na(class_labels)
df <- df[common_genes,na_filter]
class_labels <- class_labels[na_filter]

  ## version 2
data(GSE14764_eset)

eset <- GSE14764_eset
df <- exprs(eset)
class_labels <- phenoData(eset)@data$summarygrade

data(GSE17260_eset)

eset <- GSE17260_eset
df <- exprs(eset)
class_labels <- phenoData(eset)@data$summarygrade

data(GSE9891_eset)

eset <- GSE9891_eset
df <- exprs(eset)
class_labels <- phenoData(eset)@data$summarygrade
na_filter <- !is.na(class_labels)
df <- df[,na_filter]14
class_labels <- class_labels[na_filter]


###### stage I lung adenocarcinoma (GSE20189): 73 adenocarcinoma cases, 80 controls

# Retrieve file that contains URL for each dataset
download.file("https://osf.io/rcxmf/download", destfile = "URLs.txt")

# Extract URLs for datasets
urls <- read.table("URLs.txt", sep = "\t", header = F, stringsAsFactors = F, row.names = 1)

# Loop through URLs, download and extract the data to local directories
options(timeout = 10000000)
for (i in c(6, 16, 23)) {
  dataSetName <- rownames(urls)[i]
  url <- urls[i,1]
  
  destFilePath <- paste0(dataSetName, ".tar.gz")
  download.file(url, destfile = destFilePath, mode = "wb")
  
  destDirPath <- paste0(getwd(), "/", dataSetName)
  untar(destFilePath, exdir=destDirPath)
  file.remove(destFilePath)
}

file.remove("URLs.txt")

parseData <- function(geoID) {
  # Get the available analysis files, which contain class and covariate information
  analysisFiles <- list.files(path=paste0(getwd(), "/", geoID, "/Analysis"), pattern = ("txt$"), full.names = T)
  
  # Retrieve class and covariate information for one of the analysis files (as a demo)
  metaData <- read_tsv(analysisFiles[1], show_col_types = FALSE)
  
  # Retrieve the expression data
  expressionData <- suppressWarnings(read_tsv(paste0(getwd(), "/", geoID, "/", geoID, "_Expression.txt.gz"), show_col_types = FALSE))
  colnames(expressionData)[1] <- "SampleID"
  
  # Merge the class and covariate data with the expression data
  mergedData <- inner_join(metaData, expressionData, by="SampleID")
  
  # Remove the SampleID column, which is no longer needed
  mergedData <- mergedData[-which(colnames(mergedData) == "SampleID")]
  
  # Convert the tibble to a simple data frame because that's what mlr likes
  mergedData <- as.data.frame(mergedData)
  
  # Convert the Class column to a factor
  mergedData$Class <- factor(mergedData$Class)
  
  return(mergedData)
}

data <- list()

# Iterate through datasets and classification algorithms
for (geoID in rownames(urls)[c(6, 16, 23)]) {
  Sys.setenv("VROOM_CONNECTION_SIZE" = "1000000") # Increase connection buffer size
  data[[geoID]] <- parseData(geoID)
}

print("Lung adenocarcinoma")
summary(data[[1]][1:5])
dim(data[[1]])
print("Autism")
summary(data[[2]][1:5])
dim(data[[2]])
print("Oral squamous cell carcinoma")
summary(data[[3]][1:5])
dim(data[[3]])

### Lung adenocarcinoma
df <- as.data.frame(t(data[[1]] %>% dplyr::select(-c(Class, Smoking_Status.Former, Smoking_Status.Never))))
first <- function(x){x[[1]]}
annots <- mapIds(org.Hs.eg.db, keys = rownames(df), column = "SYMBOL", keytype = "ENSEMBL", multiVals = first)
df <- df[!is.na(annots),]
rownames(df) <- annots[!is.na(annots)]
class_labels <- vapply(data[[1]]$Class, FUN  = function(x){if(x == "Case"){1}else{0}}, FUN.VALUE = 0)

### Autism
df <- as.data.frame(t(data[[2]] %>% dplyr::select(-c(Class, Subject_Age, `Father__Age-Years`, `Mother__Age-Years`))))
first <- function(x){x[[1]]}
annots <- mapIds(org.Hs.eg.db, keys = rownames(df), column = "SYMBOL", keytype = "ENSEMBL", multiVals = first)
#annots[!is.na(annots)][annots[!is.na(annots)] == "ASPRV1"]
#annots[!is.na(annots)][annots[!is.na(annots)] == "FAM163A"]
annots["ENSG00000179818"] <- NA
annots["ENSG00000243062"] <- NA
df <- df[!is.na(annots),]
rownames(df) <- annots[!is.na(annots)]
class_labels <- vapply(data[[2]]$Class, FUN  = function(x){if(x == "autism"){1}else{0}}, FUN.VALUE = 0)

#Golightly, N. P., Bell, A., Bischoff, A. I., Hollingsworth, P. D., & Piccolo, S. R. (2018). Curated compendium of human transcriptional biomarker data. Scientific data, 5(1), 1-8.
#https://osf.io/ssk3t/
#adapted from https://osf.io/4n62k

###### GSEA data

preprocess_GSEA <- function(file, log = FALSE, from_file = TRUE){
  if(from_file){
    raw_data <- read.delim(file=file, skip=2)
  } else {
    raw_data <- file
  }
  print(head(raw_data))
  
  not_mapped <- c()
  for(i in raw_data$NAME){
    not_mapped <- c(not_mapped, as.logical(length(grep("_at", i))))
  }
  print(paste("# Probes not mapped: ", sum(not_mapped), "/", length(not_mapped), sep = ""))
  raw_data <- raw_data[!not_mapped, ]
  
  raw_named_data <- data.frame(t(raw_data %>% dplyr::select(-c(DESCRIPTION, NAME))))
  colnames(raw_named_data) <- raw_data$NAME
  class_labels <- vapply(rownames(raw_named_data), FUN = function(x) {strsplit(x, '_')[[1]][1]}, FUN.VALUE = "", USE.NAMES = FALSE)
  
  duplicated_columns <- duplicated(as.list(raw_named_data))
  
  print(paste("# Duplicate probe intensities: ", sum(duplicated_columns), "/", length(duplicated_columns), sep = ""))
  
  print(paste('Missing values: ', sum(is.na(raw_named_data))))
  
  qn_data_ALL <- data.frame(normalize.quantiles(as.matrix(raw_named_data[class_labels == unique(class_labels)[1],]), keep.names = TRUE))
  qn_data_AML <- data.frame(normalize.quantiles(as.matrix(raw_named_data[class_labels == unique(class_labels)[2],]), keep.names = TRUE))
  qn_data <- rbind(qn_data_ALL, qn_data_AML)
  
  if(log){
    scaled_qn_data <- apply(1 + qn_data, MARGIN = 2, FUN = log2)
  } else {
    #cbrt_qn_data <- apply(qn_data, MARGIN = 2, FUN = unction(x) {sign(x) * abs(x)^(1/3)})
    scaled_qn_data <- apply(qn_data, MARGIN = 2, FUN = scale)
    rownames(scaled_qn_data) <- rownames(qn_data)
  }
  
  return(list(data = t(scaled_qn_data), labels = class_labels))
}

  # TP53 
data("p53DataSet")

df <- p53DataSet
class_labels <- c(rep("WT", 17), rep("MUT", 33))

  # Leukemia
#idx <- which(names(geo2kegg) == "GSE9476")
#eset <- geo2kegg[[idx]]
#df <- eset@assays@data@listData$exprs
#class_labels <- eset@colData@listData$GROUP
res <- preprocess_GSEA("Leukemia_collapsed_symbols.gct.txt")
df <- res$data
class_labels <- res$labels

  # Lung Cancer
res <- preprocess_GSEA("Lung_Mich_collapsed_symbols_common_Mich_Bost.gct.txt")
df <- res$data
class_labels <- res$labels

###### GSEAbenchmarkR data
normalize_df <- function(raw_named_data, class_labels, log2 = TRUE, log = TRUE, qn = TRUE, tqn = TRUE, center = FALSE, mean = FALSE, sqrt = FALSE, cbrt = FALSE){
  duplicated_columns <- duplicated(as.list(raw_named_data))
  
  print(paste("# Duplicate probe intensities: ", sum(duplicated_columns), "/", length(duplicated_columns), sep = ""))
  
  print(paste('Missing values: ', sum(is.na(raw_named_data))))
  
  if(qn){
    if(tqn){
      qn_data_1 <- data.frame(normalize.quantiles(as.matrix(raw_named_data[class_labels == unique(class_labels)[1],]), keep.names = TRUE))
      qn_data_2 <- data.frame(normalize.quantiles(as.matrix(raw_named_data[class_labels == unique(class_labels)[2],]), keep.names = TRUE))
      qn_data <- rbind(qn_data_1, qn_data_2)
    } else {
      qn_data <- data.frame(normalize.quantiles(as.matrix(raw_named_data, keep.names = TRUE)))
    }
  } else {
    qn_data <- raw_named_data
  }
  
  if(log){
    if(log2){
      scaled_qn_data <- apply(1 + qn_data, MARGIN = 2, FUN = log2)
    }else{
      #scaled_qn_data <- apply(qn_data, MARGIN = 2, FUN = log10)
      scaled_qn_data <- log10(1 + qn_data)
    }
  } else {
    if(mean){
      scaled_qn_data <- apply(qn_data, MARGIN = 2, FUN = function(x){(x - median(x))/(sd(x))})
    }else{
      if(sqrt){
        scaled_qn_data <- apply(qn_data, MARGIN = 2, FUN = function(x) {abs(x)^(1/2)})
      }else{
        if(cbrt){
          scaled_qn_data <- apply(qn_data, MARGIN = 2, FUN = function(x) {sign(x) * abs(x)^(1/3)})
        }else{
          scaled_qn_data <- apply(qn_data, MARGIN = 2, FUN = scale, center = center)
          rownames(scaled_qn_data) <- rownames(qn_data)
        }
      }
    }
  }
  return(data.frame(t(scaled_qn_data)))
}

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

  # GSE19188: 91 tumor- and 65 adjacent normal lung tissue samples.
idx <- which(names(geo2kegg) == "GSE19188")

  # GSE42057: 42 control subjects and 94 subjects with varying severity of COPD had PBMC gene expression profiles generated.
idx <- which(names(geo2kegg) == "GSE42057")

  # GSE18842: 91 samples studied, 46 non-small cell lung cancer and 45 controls.
idx <- which(names(geo2kegg) == "GSE18842")

  # GSE11906: 42 healthy non-smokers, 49 healthy smokers, 11 symptomatic smokers, 22 smokers with lone emphysema with normal spirometry, and 20 smokers with COPD.
idx <- which(names(geo2kegg) == "GSE11906")

  # GSE15471: 35 pancreatic cancer patients and 35 controls
idx <- which(names(geo2kegg) == "GSE15471")

  # GSE9476: normal hematopoietic cells from 38 healthy donors and leukemic blasts from 26 AML patients.
idx <- which(names(geo2kegg) == "GSE9476")

  # GSE23878: 35 colorectal cancer samples versus 24 normal samples. (19-19 in this package)
idx <- which(names(geo2kegg) == "GSE23878")

  # GSE20291: two groups of age and gender matched groups of 14 Parkinson and 19 Control subjects
idx <- which(names(geo2kegg) == "GSE20291")

  # GSE32676: 42 (25 in this package) human PDAC tumors and 7 non-malignant pancreas samples 
idx <- which(names(geo2kegg) == "GSE32676")

  # GSE5281_VCX: 19 Alzheimer's and 12 normal-aged brain samples
idx <- which(names(geo2kegg) == "GSE5281_VCX")

  # GSE8762: 12 moderate stage Huntington's patients (8 female and 4 male) and 10 age-matched controls (5 female and 5 male)
idx <- which(names(geo2kegg) == "GSE8762")

  # GSE16515: 36 Pancreatic tumor samples and 16 normal samples
idx <- which(names(geo2kegg) == "GSE16515")

  # GSE38666_epithelia: ovarian cancer
idx <- which(names(geo2kegg) == "GSE38666_epithelia")

  # GSE1145: heart failure
idx <- which(names(geo2kegg) == "GSE1145")

  # GSE30153: 17 patients with quiescent lupus versus 9 controls
idx <- which(names(geo2kegg) == "GSE30153")

  # GSE19420: Type 2 diabetes
idx <- which(names(geo2kegg) == "GSE19420")

  # GSE5281_HIP: Alzheimer's
idx <- which(names(geo2kegg) == "GSE5281_HIP")

  # GSE4183: Colon cancer
idx <- which(names(geo2kegg) == "GSE4183")

  # GSE4107: Colorectal cancer
idx <- which(names(geo2kegg) == "GSE4107")

  # GSE8762: Lymphocyte samples from 12 moderate stage Huntington's patients (8 female and 4 male) and 10 age-matched controls (5 female and 5 male).
idx <- which(names(geo2kegg) == "GSE8762")

  # GSE5281_EC: Alzheimer's
idx <- which(names(geo2kegg) == "GSE5281_EC")

  # GSE1297: hippocampal gene expression of 9 control and 22 Incipient Alzheimer's subjects of varying severity
idx <- which(names(geo2kegg) == "GSE1297") 

  # GSE20153: 8 Parkinson lymphobast cell lines versus 8 control cell lines
idx <- which(names(geo2kegg) == "GSE20153")

  # GSE20164: Substantia nigra samples from 6 Parkinson's and 5 control subjects were obtained
idx <- which(names(geo2kegg) == "GSE20164")

  # GSE16759: parietal lobe tissue from 4 AD patients and 4 age-matched controls
idx <- which(names(geo2kegg) == "GSE16759")

  # General
eset <- geo2kegg[[idx]]
df <- eset@assays@data@listData$exprs
class_labels <- as.character(eset@colData@listData$GROUP)
gene_symbols <- unlist(lookUp(eset@NAMES, 'org.Hs.eg', 'SYMBOL'))
sum(rownames(df) != names(gene_symbols))
rownames(df)<- gene_symbols
#ggplot2::ggplot(data.frame(t(df)), ggplot2::aes(x = !!as.name(rownames(df)[150]))) + ggplot2::geom_density()
#df <- normalize_df(data.frame(t(df)), class_labels, qn = FALSE, tqn = TRUE, log = TRUE, log2 = TRUE, center = FALSE, mean = FALSE, sqrt = FALSE, cbrt = FALSE)
#ggplot(data.frame(t(df)), aes(x = !!as.name(rownames(df)[150]))) + geom_density()

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

sum(rowSums(KEGG_A) == 0) # -> 2657, 2658, 2683 (p53), 2681 (ovarian 2,3), 2760 (leukemia), 2760 (GSE18842), 2681 (GSE15471), 2681 (GSE9476), 2681 (GSE42057), 2681 (GSE20291), 2681 (GSE5281_VCX), 2681 (GSE20189), 2681 (GSE25507), 2681 (GSE8762)
#saveRDS(KEGG_A, "KEGG_A_GE13507.rds")
#saveRDS(KEGG_A, "KEGG_A_leukemia.rds")
#saveRDS(KEGG_A, "KEGG_A_GE18842.rds")
#saveRDS(KEGG_A, "KEGG_A_GE15471.rds")
#saveRDS(KEGG_A, "KEGG_A_GE9476.rds")
#saveRDS(KEGG_A, "KEGG_A_GE42057.rds")
#saveRDS(KEGG_A, "KEGG_A_GE20291.rds")
#saveRDS(KEGG_A, "KEGG_A_GE5281VCX.rds")
#saveRDS(KEGG_A, "KEGG_A_GE20189.rds")
#saveRDS(KEGG_A, "KEGG_A_GE25507.rds")
saveRDS(KEGG_A, "KEGG_A_GE8762.rds")


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
sum(rowSums(KEGG_A) == 0) # -> 1806, 1809, 1832 (p53), 1832 (ovarian 1,3), 1914 (leukemia), 1914 (GSE18842), 1832 (GSE15471), 1832 (GSE9476), 1832 (GSE42057), 1832 (GSE20291), 1832 (GSE5281_VCX), 1832 (GSE20189), 1832 (GSE25507), 1832 (GSE8762)
#saveRDS(KEGG_A, "KEGG_A_undir_GE13507.rds")
#saveRDS(KEGG_A, "KEGG_A_undir_leukemia.rds")
#saveRDS(KEGG_A, "KEGG_A_undir_GE18842.rds")
#saveRDS(KEGG_A, "KEGG_A_undir_GE15471.rds")
#saveRDS(KEGG_A, "KEGG_A_undir_GE9476.rds")
#saveRDS(KEGG_A, "KEGG_A_undir_GE42057.rds")
#saveRDS(KEGG_A, "KEGG_A_undir_GE20291.rds")
#saveRDS(KEGG_A, "KEGG_A_undir_GE5281VCX.rds")
#saveRDS(KEGG_A, "KEGG_A_undir_GE20189.rds")
#saveRDS(KEGG_A, "KEGG_A_undir_GE25507.rds")
saveRDS(KEGG_A, "KEGG_A_undir_GE8762.rds")

##########################
# Reading RDS file to CSV
##########################

am_rds <- readRDS("KEGG_A_undir_GE14764.rds")
write.csv(am_rds, "KEGG_A_undir_GE14764.csv", row.names=TRUE)
am_csv <- read.csv("KEGG_A_undir_GE14764.csv")
# View(am_csv)
# View(am_rds)

am_rds <- readRDS("KEGG_A_undir_GE31684.rds")
write.csv(am_rds, "KEGG_A_undir_GE31684.csv", row.names=TRUE)
am_csv <- read.csv("KEGG_A_undir_GE31684.csv")

am_rds <- readRDS("KEGG_A_undir_p53.rds")
write.csv(am_rds, "KEGG_A_undir_p53.csv", row.names=TRUE)
am_csv <- read.csv("KEGG_A_undir_p53.csv")

am_rds <- readRDS("KEGG_A_undir_GE17260.rds")
write.csv(am_rds, "KEGG_A_undir_GE17260.csv", row.names=TRUE)
am_csv <- read.csv("KEGG_A_undir_GE17260.csv")

am_rds <- readRDS("KEGG_A_undir_GE9891.rds")
write.csv(am_rds, "KEGG_A_undir_GE9891.csv", row.names=TRUE)
am_csv <- read.csv("KEGG_A_undir_GE9891.csv")

am_rds <- readRDS("KEGG_A_undir_leukemia.rds")
write.csv(am_rds, "KEGG_A_undir_leukemia.csv", row.names=TRUE)
am_csv <- read.csv("KEGG_A_undir_leukemia.csv")

am_rds <- readRDS("KEGG_A_undir_GE18842.rds")
write.csv(am_rds, "KEGG_A_undir_GE18842.csv", row.names=TRUE)
am_csv <- read.csv("KEGG_A_undir_GE18842.csv")

am_rds <- readRDS("KEGG_A_undir_GE15471.rds")
write.csv(am_rds, "KEGG_A_undir_GE15471.csv", row.names=TRUE)
am_csv <- read.csv("KEGG_A_undir_GE15471.csv")

am_rds <- readRDS("KEGG_A_undir_GE9476.rds")
write.csv(am_rds, "KEGG_A_undir_GE9476.csv", row.names=TRUE)
am_csv <- read.csv("KEGG_A_undir_GE9476.csv")

am_rds <- readRDS("KEGG_A_undir_GE42057.rds")
write.csv(am_rds, "KEGG_A_undir_GE42057.csv", row.names=TRUE)
am_csv <- read.csv("KEGG_A_undir_GE42057.csv")

am_rds <- readRDS("KEGG_A_undir_GE20291.rds")
write.csv(am_rds, "KEGG_A_undir_GE20291.csv", row.names=TRUE)
am_csv <- read.csv("KEGG_A_undir_GE20291.csv")

am_rds <- readRDS("KEGG_A_undir_GE5281VCX.rds")
write.csv(am_rds, "KEGG_A_undir_GE5281VCX.csv", row.names=TRUE)
am_csv <- read.csv("KEGG_A_undir_GE5281VCX.csv")

am_rds <- readRDS("KEGG_A_undir_GE20189.rds")
write.csv(am_rds, "KEGG_A_undir_GE20189.csv", row.names=TRUE)
am_csv <- read.csv("KEGG_A_undir_GE20189.csv")

am_rds <- readRDS("KEGG_A_undir_GE25507.rds")
write.csv(am_rds, "KEGG_A_undir_GE25507.csv", row.names=TRUE)
am_csv <- read.csv("KEGG_A_undir_GE25507.csv")

am_rds <- readRDS("KEGG_A_undir_GE8762.rds")
write.csv(am_rds, "KEGG_A_undir_GE8762.csv", row.names=TRUE)
am_csv <- read.csv("KEGG_A_undir_GE8762.csv")

#############################################
# Importing LCC CSV & checking AM correctness
#############################################

#lcc<- as.matrix(read.csv('./KEGG_A_undir_GE13507_LCC.csv', row.names = 1))
#lcc<- as.matrix(read.csv('./KEGG_A_undir_p53_LCC.csv', row.names = 1))
#lcc<- as.matrix(read.csv('./KEGG_A_undir_GE14764_LCC.csv', row.names = 1))
#lcc<- as.matrix(read.csv('./KEGG_A_undir_GE17260_LCC.csv', row.names = 1))
#lcc<- as.matrix(read.csv('./KEGG_A_undir_GE9891_LCC.csv', row.names = 1))
#lcc<- as.matrix(read.csv('./KEGG_A_undir_leukemia_LCC.csv', row.names = 1))
#lcc<- as.matrix(read.csv('./KEGG_A_undir_GE18842_LCC.csv', row.names = 1))
#lcc<- as.matrix(read.csv('./KEGG_A_undir_GE15471_LCC.csv', row.names = 1))
#lcc<- as.matrix(read.csv('./KEGG_A_undir_GE9476_LCC.csv', row.names = 1))
#lcc<- as.matrix(read.csv('./KEGG_A_undir_GE42057_LCC.csv', row.names = 1))
#lcc<- as.matrix(read.csv('./KEGG_A_undir_GE20291_LCC.csv', row.names = 1))
#lcc<- as.matrix(read.csv('./KEGG_A_undir_GE5281VCX_LCC.csv', row.names = 1))
#lcc<- as.matrix(read.csv('./KEGG_A_undir_GE20189_LCC.csv', row.names = 1))
#lcc<- as.matrix(read.csv('./KEGG_A_undir_GE25507_LCC.csv', row.names = 1))
lcc<- as.matrix(read.csv('./KEGG_A_undir_GE8762_LCC.csv', row.names = 1))

#View(lcc)
ncol(lcc)
nrow(lcc)
sum(rownames(lcc) != colnames(lcc))
sum(!rownames(lcc) %in% colnames(lcc))
sum(!colnames(lcc) %in% rownames(lcc))

rownames(lcc)[!rownames(lcc) %in% colnames(lcc)] <- sub('-', '.', rownames(lcc)[!rownames(lcc) %in% colnames(lcc)])
sum(rownames(lcc) != colnames(lcc))

#########################
#     ACO
#########################
#sd_not0 <- c()
#for(i in c(1:nrow(df))){ 
#  if(sd(df[i,]) != 0){ 
#    sd_not0 <- c(sd_not0, i)
#  }
#}
#df_filtered <- df[sd_not0,]
#KEGG_A <- readRDS("KEGG_A_undir_GE13507.rds")
#KEGG_A_LCC <- as.matrix(read.csv("KEGG_A_undir_GE13507_LCC.csv", row.names = 1))
KEGG_A_LCC <- lcc
source("ant_colony.R")
source("OSCC.R")
aco_modules <- optimize_modules(df, KEGG_A_LCC, n_ants = 5, n_iter = 5, starting_capacity = 1.0, 
                                n_cores = 1, alpha = 0.6, beta = 1.2, seed = 3, 
                                n_resamples = 100, class_labels = class_labels)

saveRDS(aco_modules, 'aco_beta_trial_9891_5.rds')

#########################
#     Parallelized ACO
#########################

library(parallel)
source("ant_colony.R")
source("OSCC.R")

max_cores <- detectCores()

KEGG_A_LCC <- lcc                                     
aco_modules <- optimize_modules(df, KEGG_A_LCC, n_ants = 40, n_iter = 5, starting_capacity = 1.0, 
                                n_cores = min(5, max_cores - 2), alpha = 0.6, beta = 1.2, 
                                seed = 1, n_resamples = 100, class_labels = class_labels, penalize_dist = TRUE)

saveRDS(aco_modules, 'aco_gamma_trial_8762_40_1.rds')
#gamma: n_iter = 5, starting_capacity = 1.0, n_cores = min(5, max_cores - 2), alpha = 0.6, beta = 1.2, n_resamples = 100, penalize_dist = TRUE
#delta: n_iter = 5, starting_capacity = 1.0, n_cores = min(5, max_cores - 2), alpha = 0.6, beta = 1.2, n_resamples = 1000, penalize_dist = TRUE
#omega: n_iter = 10, starting_capacity = 1.0, n_cores = min(5, max_cores - 2), alpha = 0.6, beta = 1.2, n_resamples = 100, penalize_dist = TRUE
#sigma: n_iter = 5, starting_capacity = 1.0, n_cores = min(5, max_cores - 2), alpha = 0.6, beta = 1.2, n_resamples = 100, penalize_dist = FALSE

#########################
# Preprocess ACO results
#########################

aco_modules <- readRDS('aco_gamma_trial_9476_25.rds')

limma_res <- acoCalcGeneTStats(df, class_labels, 100, 3)
hist(limma_res$complete_res$adj.P.Val)
limma_adj_pvals <- limma_res$complete_res$adj.P.Val
names(limma_adj_pvals) <- rownames(limma_res$complete_res)

all_stats <- data.frame()
min_stats <- data.frame()
tot_time <- c()
tot_time_par <- c()
min_oscc <- c()
min_pdecay <- c()
de_pvals <- c()
modules <- list()
filter_by <- 'p_decays'
cutoff <- 0.05
options(timeout=200)
for(gene_name in names(aco_modules)){
  gene <- aco_modules[[gene_name]]
  de_pvals <- c(de_pvals, limma_adj_pvals[gene_name])
  all_stats <- rbind(all_stats, gene$stats)
  min_stats <- rbind(min_stats, gene$stats[which(gene$stats[,filter_by] == min(gene$stats[,filter_by]))[1],])
  tot_time_par <- c(tot_time_par, gene$total_time)#sum(gene$stats[,'pDecay_times']))
  tot_time <- c(tot_time, sum(gene$stats[,'pDecay_times']))
  min_oscc <- c(min_oscc, gene$stats[which(gene$stats[,filter_by] == min(gene$stats[,filter_by]))[1],'osccs'])
  min_pdecay <- c(min_pdecay,gene$stats[which(gene$stats[,filter_by] == min(gene$stats[,filter_by]))[1],'p_decays'])
  modules[[gene_name]] <- gene$modules$modules[[which(gene$stats[,filter_by] == min(gene$stats[,filter_by]))[1]]]
}
colnames(min_stats) <- colnames(all_stats)
rownames(min_stats) <- names(aco_modules)
min_stats_new <- min_stats %>% mutate(p_fishers_log = -log10(p_fishers), 
                                      bonf_adj = p_fishers_log >= -log10(cutoff/5),
                                      p_fishers_bonf = vapply(p_fishers*5, FUN = min, 1),
                                      p_des = de_pvals,
                                      p_combs = pchisq(-2*(log(p_fishers) + log(p_des)), df = 4, lower.tail = FALSE),#pchisq(-2*(log(p_spheres) + log(p_decays) + log(p_des)), df = 6, lower.tail = FALSE),
                                      bonf_adj_combs = p_combs <= cutoff)
                                      #p_fishers_BH = p.adjust(p_fishers, method = "BH"),
                                      #p_fishers_bonf = p.adjust(p_fishers, method = "bonferroni"),
                                      #bonf_adj = p_fishers_bonf <= 0.05,
                                      #BH_adj = p_fishers_BH <= 0.05)

sum(min_stats_new$bonf_adj_combs)

aco <- rownames(min_stats_new %>% filter(p_fishers <= cutoff) %>% arrange(p_fishers))
aco <- rownames(min_stats_new %>% filter(bonf_adj_combs) %>% arrange(p_combs))
#tf <- treat(lmFit(df, model.matrix(~factor(class_labels))))
#aco2 <- aco[aco %in% rownames(tf$p.value[tf$p.value[,2] <= 0.05,])]

aco_sig_genes <- rownames(min_stats_new %>% filter(bonf_adj == TRUE) %>% arrange(p_fishers))
#saveRDS(aco_sig_genes, 'aco_sig_genes_42057_25.rds')
#write.csv(aco_sig_genes, 'aco_sig_genes_42057_25.csv', row.names = FALSE)


overlap_mtrx <- matrix(nrow = length(modules), ncol = length(modules), dimnames = list(rownames = names(modules), colnames = names(modules)))
#overlap_mtrx <- readRDS("overlap_matrix_18842.rds")
n <- length(modules)
for(id1 in 1:n){
  gene1 <- names(modules)[id1]
  module1 <- modules[[gene1]]
  len1 <- length(module1)
  for(id2 in 1:n){
    if(id2 >= id1){
      gene2 <- names(modules)[id2]
      module2 <- modules[[gene2]]
      len2 <- length(module2)
      overlap_mtrx[id1, id2] <- (sum(module1 %in% module2)**2)/(len1*len2)
    } else {
      overlap_mtrx[id1, id2] <- overlap_mtrx[id2, id1]
    }
  }
}

overlap_mtrx_sig <- overlap_mtrx[rownames(overlap_mtrx) %in% aco, colnames(overlap_mtrx) %in% aco]
#saveRDS(overlap_mtrx, "overlap_matrix_18842.rds")
#saveRDS(overlap_mtrx, "overlap_matrix_42057.rds")

heatmap(overlap_mtrx_sig,Rowv=NA,Colv=NA,reorderfun=NA,hclustfun=NA)

image(overlap_mtrx_sig,
      main = "GSE18842 module overlap")

#########################
#     GS
#########################

source("old_genesurrounder/GeneSurrounder.R")
source("old_genesurrounder/run_geneSurrounder.R")
gs_modules <- run_geneSurrounder(adj.matrix = KEGG_A_LCC,
                                 gene_scores_per_sample = df,
                                 class_labels = class_labels,
                                 nresamples = 100,
                                 seed = 121,
                                 num.Sphere.resamples = 100,
                                 gene.id = names(aco_modules),
                                 decay_only = FALSE,
                                 file_name = "gs_results_8762_1.csv",
                                 cores = 1 # Set to 1 for Windows
)
saveRDS(gs_modules, 'gs_results_8762_1.rds')

#########################
#  Preprocess GS results
#########################

gs_modules <- readRDS('gs_results_20291.rds')
diam <- 18
gs_modules_new <- gs_modules %>% mutate(p.Fisher_log = -log10(p.Fisher), 
                                        bonf_adj = p.Fisher_log >= -log10(0.05/diam))
                                        #p_fishers_bonf = p.adjust(p.Fisher, method = "bonferroni"),
                                        #p_fishers_BH = p.adjust(p.Fisher, method = "BH"),
                                        #bonf_adj = p_fishers_bonf <= 0.05,
                                        #BH_adj = p_fishers_BH <= 0.05)
sum(gs_modules_new$bonf_adj)

gs_sig_genes <- (gs_modules_new %>% filter(bonf_adj == TRUE) %>% dplyr::select(gene.id))[,1]

saveRDS(gs_sig_genes, 'gs_sig_genes_20291.rds')

saveRDS(gs_modules$gene.id, 'universe_genes_20291.rds')

#########################
#     Analysis
#########################


######p_decay vs score
#p_decay vs tau-b
pdecay_taub <- ggplot(gs_modules, aes(observed.tau_b, p.Decay)) +
  geom_point() +
  ggtitle("GS")

#p_decay vs oscc
pdecay_oscc <- ggplot(all_stats, aes(osccs, p_decays)) +
  geom_point() +
  ggtitle("ACO")

gridExtra::grid.arrange(pdecay_taub, pdecay_oscc, ncol = 2)

######p_decay distribution
p_dist <- ggplot(gs_modules, aes(p.Decay)) + 
  geom_histogram(bins = 50) +
  ggtitle("GS")

p_dist_zoom <- ggplot(gs_modules, aes(p.Decay)) + 
  geom_histogram(bins = 10) + 
  xlim(0, 0.1) + 
  ggtitle("GS")

p_dist2 <- ggplot(all_stats, aes(p_decays)) + 
  geom_histogram(bins = 50) + 
  ggtitle("ACO")

p_dist_zoom2 <- ggplot(all_stats, aes(p_decays)) + 
  geom_histogram(bins = 10) + 
  xlim(0, 0.1) +
  ggtitle("ACO")

gridExtra::grid.arrange(p_dist, p_dist_zoom, p_dist2, p_dist_zoom2, ncol = 2)

ggplot(data.frame(gs_pdecay = gs_modules[,'p.Decay'], aco_min_pdecay = min_pdecay), aes(gs_pdecay, aco_min_pdecay)) +
  geom_point() + 
  ggtitle("ACO minimum p-decay vs GS p-decay")

######Score distributions
#tau-b distribution
taubs <- ggplot(gs_modules, aes(observed.tau_b)) + 
  geom_histogram(bins = 20) + 
  ggtitle("GS")

#oscc distribution
osccs <- ggplot(all_stats, aes(osccs)) + 
  geom_histogram(bins = 20) + 
  ggtitle("ACO")

gridExtra::grid.arrange(taubs, osccs, ncol = 2)

ggplot(data.frame(tau_b = gs_modules[,'observed.tau_b'], min_oscc = min_oscc), aes(tau_b, min_oscc)) + 
  geom_point() + 
  ggtitle("tau-b vs minimum p-decay OSCC")

######ACO time/iteration
ggplot(all_stats, aes(pDecay_times)) + 
  geom_density() +
  xlab('times (seconds) / module') +
  ggtitle("ACO time/iteration")

######time per gene
gs_times <- ggplot(gs_modules, aes(time)) + 
  geom_density() +
  xlab('times (seconds) / gene') +
  ggtitle("GS time/gene")

aco_times <- ggplot(data.frame(pdecay_time = tot_time), aes(pdecay_time)) + 
  geom_density() +
  xlab('times (seconds) / gene (no parallelization)') +
  ggtitle("ACO time/gene")

aco_times_par <- ggplot(data.frame(pdecay_time = tot_time_par), aes(pdecay_time)) + 
  geom_density() +
  xlab('times (seconds) / gene (parallelized)') +
  ggtitle("Parallelized ACO time/gene")

gridExtra::grid.arrange(gs_times, aco_times, aco_times_par, ncol = 1)

gs_vs_aco <- ggplot(data.frame(GS = gs_modules[,'pdecay_time'], ACO = tot_time), aes(GS, ACO)) +
  geom_point() + 
  ggtitle("GS vs Unparallelized ACO time/gene")

gs_vs_paco <- ggplot(data.frame(GS = gs_modules[,'pdecay_time'], ACO = tot_time_par), aes(GS, ACO)) +
  geom_point() + 
  ggtitle("GS vs Parallelized ACO time/gene")

paco_vs_aco <- ggplot(data.frame(Unparallelized_ACO = tot_time, Parallelized_ACO = tot_time_par), aes(Parallelized_ACO, Unparallelized_ACO)) +
  geom_point() + 
  ggtitle("Unparallelized vs Parallelized ACO time/gene")

gridExtra::grid.arrange(gs_vs_aco, gs_vs_paco, paco_vs_aco, ncol = 1)

######Total time
print(paste("GS total time (seconds):", sum(gs_modules[,"pdecay_time"])))
print(paste("Unparallelized ACO total time (seconds):", sum(tot_time)))
print(paste("Parallelized ACO total time (seconds):", sum(tot_time_par)))

######Module Graphs
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
lean_modules <- runLEAN(df, KEGG_A_LCC, n_cores = NULL, seed = 3, 
                        n_resamples = 100, class_labels = class_labels)

saveRDS(lean_modules, 'lean_modules_53.rds')
