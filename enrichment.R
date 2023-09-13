source("ant_colony.R")
library(GSEABenchmarkeR) # GSEABenchmark data
library(org.Hs.eg.db)

geo2kegg <- loadEData("geo2kegg", preproc = TRUE)

ids <- c()
sizes <- c()
for(id in names(geo2kegg)){
  ids<- c(ids, id)
  idx <- which(names(geo2kegg) == id)
  eset <- geo2kegg[[idx]]
  sizes <- c(sizes, dim(eset)[2])
}
bm_ds_ids <- data.frame(ids, sizes) %>% dplyr::arrange(desc(sizes))

dataset_id <- "8762"
dataset <- paste("GSE", dataset_id, sep = "")
n_ants <- 40
seed <- 3
filter_by <- 'ADJ.PVAL_filter'
trial <- 'gamma'
cutoff <- 0.05
idx <- which(names(geo2kegg) == dataset)
aco_modules <- readRDS(paste('aco', trial, 'trial', dataset_id, n_ants, paste(seed,'rds', sep = '.'), sep = '_'))
#aco_modules <- readRDS('aco_gamma_trial_p53_40.rds')
gs_modules <- readRDS(paste('gs', 'results', dataset_id, paste(seed, 'rds', sep = '.'), sep = '_'))

# TP53 
data("p53DataSet")

df <- p53DataSet
class_labels <- c(rep("0", 17), rep("1", 33))

#GSEABenchmarkeR
eset <- geo2kegg[[idx]]
df <- eset@assays@data@listData$exprs
class_labels <- as.character(eset@colData@listData$GROUP)
uid_to_symbol <- lookUp(eset@NAMES, 'org.Hs.eg', 'SYMBOL')
gene_symbols <- unlist(uid_to_symbol)
rownames(df)<- gene_symbols

lcc <- as.matrix(read.csv(paste('./KEGG', 'A', 'undir', paste('GE', dataset_id, sep = ''), 'LCC.csv', sep = '_'), row.names = 1))
rownames(lcc)[!rownames(lcc) %in% colnames(lcc)] <- sub('-', '.', rownames(lcc)[!rownames(lcc) %in% colnames(lcc)])
adj_mat <- lcc
KEGG_A_LCC <- lcc
adj_mat <- adj_mat[rownames(adj_mat) %in% rownames(df), ]
adj_mat <- adj_mat[ ,colnames(adj_mat) %in% rownames(df)]
int_network <- igraph::graph_from_adjacency_matrix(adj_mat, mode = 'undirected')
diam <- diameter(int_network)

limma_res <- acoCalcGeneTStats(df, class_labels, 100, seed)
hist(limma_res$complete_res$adj.P.Val, main = "Histogram of Limma adjusted p-values for p53 dataset")
hist(limma_res$complete_res$t, main = "Histogram of Limma moderated t-statistics for p53 dataset")
limma_adj_pvals <- limma_res$complete_res$adj.P.Val
names(limma_adj_pvals) <- rownames(limma_res$complete_res)

aco_df <- df[rownames(df) %in% names(aco_modules),]

all_stats <- data.frame()
min_stats <- data.frame()
tot_time <- c()
tot_time_par <- c()
min_oscc <- c()
de_pvals <- c()
min_pdecay <- c()
adj_method <- 'BH'
get_min_by <- 'p_decays'
options(timeout=200)
for(gene_name in names(aco_modules)){
  gene <- aco_modules[[gene_name]]
  gene$stats <- as.matrix(data.frame(gene$stats) %>% mutate(adj_p_fishers = p.adjust(p_fishers, method = adj_method)))
  de_pvals <- c(de_pvals, limma_adj_pvals[gene_name])
  all_stats <- rbind(all_stats, gene$stats)
  min_stats <- rbind(min_stats, gene$stats[which(gene$stats[,get_min_by] == min(gene$stats[,get_min_by]))[1],])
  tot_time_par <- c(tot_time_par, gene$total_time)#sum(gene$stats[,'pDecay_times']))
  tot_time <- c(tot_time, sum(gene$stats[,'pDecay_times']))
  min_oscc <- c(min_oscc, gene$stats[which(gene$stats[,get_min_by] == min(gene$stats[,get_min_by]))[1],'osccs'])
  min_pdecay <- c(min_pdecay,gene$stats[which(gene$stats[,get_min_by] == min(gene$stats[,get_min_by]))[1],'p_decays'])
}
colnames(min_stats) <- colnames(all_stats)
rownames(min_stats) <- rownames(aco_df)

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

lean_modules <- runLEAN(df, KEGG_A_LCC, n_cores = NULL, seed = seed, 
                        n_resamples = 100, class_labels = class_labels)
#lean_modules <- LEANR::run.lean((min_stats_new %>% arrange(p_combs) %>% dplyr::select(p_combs))[,1], int_network, verbose = TRUE, ranked = FALSE, n_reps = 100, ncores = NULL)
lean <- lean_modules$restab[rownames(lean_modules$restab) %in% names(aco_modules),8]
lean <- lean[order(factor(names(lean), levels = names(aco_modules)))]

apply_comb_to <- 'adj_p_fishers' #'adj_p_fishers'
#k <- 1
min_stats_new <- min_stats %>% mutate(FC = apply(aco_df[,!as.logical(as.numeric(class_labels))], 1, mean) - apply(aco_df[,as.logical(as.numeric(class_labels))], 1, mean),
                                      p_des = de_pvals,
                                      p_lean = lean,
                                      fisher_stat = -2*(log(!!as.name(apply_comb_to)) + log(p_des)),
                                      p_fishers_log = -log10(p_fishers), 
                                      bonf_adj = p_fishers_log >= -log10(0.05/5),
                                      p_combs = pchisq(-2*(log(!!as.name(apply_comb_to)) + log(p_des)), df = 4, lower.tail = FALSE)
                                      #BH_adj_p = p.adjust(!!as.name(apply_BH_to), method = "BH", n = length(!!as.name(apply_BH_to))*k) #which adjustment to use?
)


####################################################
#######   Using DOSE and DO pathways        ########
####################################################

#DO enrichment
#install.packages("DOSE")
library(DOSE)

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("geneset")
#BiocManager::install("GSEABenchmarkeR")
#BiocManager::install("clusterProfiler")
#BiocManager::install("KEGG.db")

library(org.Hs.eg.db)

#aco <- readRDS(file = "aco_sig_genes_20291_25.rds")
aco <- rownames(min_stats_new %>% filter(adj_p_fishers <= cutoff) %>% arrange(adj_p_fishers))
gs <- gs_modules$gene.id[vapply(gs_modules$p.Fisher, FUN = function(x){min(x*diam,1)}, FUN.VALUE = 1) <= cutoff]
universe_genes <- names(aco_modules)
#gs <- readRDS(file = "gs_sig_genes_20291.rds")
#universe_genes <- readRDS(file = "universe_genes_20291.rds")


universe_genes

aco_entrez <- mapIds(org.Hs.eg.db,aco,'ENTREZID','SYMBOL')
gs_entrez <- mapIds(org.Hs.eg.db,gs,'ENTREZID','SYMBOL')
universe_entrez <- mapIds(org.Hs.eg.db,universe_genes,'ENTREZID','SYMBOL')

enr_aco <- enrichDO(
  aco_entrez,
  ont = "DO",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = universe_entrez,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  readable = TRUE
)

enr_gs <- enrichDO(
  gs_entrez,
  ont = "DO",
  pvalueCutoff = 0.05,
  pAdjustMethod = "bonferroni",
  universe = universe_entrez,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  readable = TRUE
)

enr_aco@result[["p.adjust"]]
enr_aco@result[["Description"]][1:sum(enr_aco@result[["p.adjust"]] <= 0.05)]
which(enr_aco@result[["ID"]] == "DOID:11054")

enr_gs@result[["p.adjust"]]
enr_gs@result[["Description"]][1:sum(enr_gs@result[["p.adjust"]] <= 0.05)]
which(enr_gs@result[["ID"]] == "DOID:11054")

####################################################
#######   Using DOSE and DGN pathways        #######
####################################################

enr_dgn__gs <- enrichDGN(
  gs_entrez,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = universe_entrez,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  readable = TRUE
)

enr_dgn__gs@result[["p.adjust"]]
enr_dgn__gs@result[["Description"]][1:min(200, sum(enr_dgn__gs@result[["p.adjust"]] <= 0.05))]
which(enr_dgn__gs@result[["ID"]] == "C2145472")


enr_dgn_aco <- enrichDGN(
  aco_entrez,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = universe_entrez,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  readable = TRUE
)

enr_dgn_aco@result[["p.adjust"]]
sum(enr_dgn_aco@result[["p.adjust"]] <= 0.05)
enr_dgn_aco@result[["Description"]][1:200]#sum(enr_dgn_aco@result[["p.adjust"]] <= 0.05)]
which(enr_dgn_aco@result[["ID"]] == "C2145472")

####################################################
#######   Using clusterProfiler and KEGG pathways ##
####################################################

library(clusterProfiler)

enr_kegg_aco <- clusterProfiler::enrichKEGG(
  gene = (aco_entrez)[!is.na(aco_entrez)],
  organism = "hsa",
  keyType = "kegg",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = universe_entrez,
  use_internal_data = FALSE
)

enr_kegg_aco@result[["p.adjust"]]
sum(enr_kegg_aco@result[["p.adjust"]] <= 0.05)
enr_kegg_aco@result[["Description"]][1:sum(enr_kegg_aco@result[["p.adjust"]] <= 0.05)]

which(enr_kegg_aco@result[["Description"]] == "Non-small cell lung cancer")
which(enr_kegg_aco@result[["Description"]] == "Acute myeloid leukemia")
which(enr_kegg_aco@result[["Description"]] == "p53 signaling pathway")

enr_kegg_gs <- clusterProfiler::enrichKEGG(
  gene = (gs_entrez)[!is.na(gs_entrez)],
  organism = "hsa",
  keyType = "ncbi-geneid",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = universe_entrez,
  use_internal_data = FALSE
)

enr_kegg_gs@result[["p.adjust"]]
sum(enr_kegg_gs@result[["p.adjust"]] <= 0.05)
enr_kegg_gs@result[["Description"]][1:sum(enr_kegg_gs@result[["p.adjust"]] <= 0.05)]

####################################################
#######   Using GSEABenchmarkeR and KEGG pathways ##
####################################################


library(GSEABenchmarkeR) # GSEABenchmark data
library(geneset) #KEGG pathway gene sets

res_df <- min_stats_new %>% mutate(limma.STAT = fisher_stat, PVAL = p_fishers, ADJ.PVAL = adj_p_fishers, ADJ.PVAL_filter = ADJ.PVAL <= cutoff, PVAL_filter = PVAL <= cutoff) %>% filter(!!as.name(filter_by)) %>% dplyr::select(c(FC, limma.STAT, PVAL, ADJ.PVAL))
df_sig <- df[rownames(df) %in% rownames(res_df),]

symbol_to_uid <- names(uid_to_symbol)
names(symbol_to_uid) <- gene_symbols
rownames(res_df) <- symbol_to_uid[rownames(res_df)]
rownames(df_sig) <- symbol_to_uid[rownames(df_sig)]

dim(res_df)

res_eset <- SummarizedExperiment(assays = list(exprs = df_sig), rowData = DataFrame(res_df), colData = colData(eset), metadata = metadata(eset))
#hsa_kegg <- KEGGREST::keggList("hsa")
#hsapathway <- EnrichmentBrowser::downloadPathways("hsa")
#kegg.gs <- EnrichmentBrowser::getGenesets('hsa', db = c('kegg'), cache = FALSE)
genesets_kegg <- geneset::getKEGG(org = "hsa", category = "pathway")

kegg_pathways <- genesets_kegg$geneset_name$name
names(kegg_pathways) <- genesets_kegg$geneset_name$id
kegg_gs <- list()
for(row in 1:nrow(genesets_kegg$geneset)){
  current_pathway_id <- genesets_kegg$geneset[row,1]
  current_gene <- genesets_kegg$geneset[row,2]
  kegg_gs[[kegg_pathways[current_pathway_id]]] <- c(kegg.gs[[kegg_pathways[current_pathway_id]]], current_gene)
}


res_EA <- runEA(res_eset, methods = 'gsa', gs = kegg_gs, padj.method = 'BH')

n_sig_pathways <- sum(res_EA$gsa[[dataset]]$ranking$ADJ.PVAL <= 0.05)
n_sig_pathways

res_EA$gsa[[dataset]]$ranking$GENE.SET[1:n_sig_pathways]

which(res_EA$padog[[dataset]]$ranking$GENE.SET == "Parkinson disease")

####################################################
#######   Using Wilcoxon Rank and KEGG pathways ####
####################################################

library(KEGGREST)
pathways.list <- KEGGREST::keggList("pathway", "hsa")
head(pathways.list)
pathway.codes <- sub("path:", "", names(pathways.list)) 
genes.by.pathway <- readRDS('genes_by_pathway.rds')
#genes.by.pathway <- sapply(pathway.codes,
#                           function(pwid){
#                             pw <- KEGGREST::keggGet(pwid)
#                             if (is.null(pw[[1]]$GENE)) return(NA)
#                             pw2 <- pw[[1]]$GENE[c(FALSE,TRUE)] # may need to modify this to c(FALSE, TRUE) for other organisms
#                             pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
#                             return(pw2)
#                           }
#)
#saveRDS(genes.by.pathway, 'genes_by_pathway.rds')

geneList <- min_stats_new$adj_p_fishers
names(geneList) <- rownames(min_stats_new)
geneList <- min_stats_new$p_lean
names(geneList) <- rownames(min_stats_new)
geneList <- min_stats_new$p_des
names(geneList) <- rownames(min_stats_new)
geneList <- vapply(gs_modules$p.Fisher, FUN = function(x){min(x*diam,1)}, FUN.VALUE = 1)
names(geneList) <- gs_modules$gene.id
# Wilcoxon test for each pathway
pVals.by.pathway <- t(sapply(names(genes.by.pathway),
                             function(pathway) {
                               pathway.genes <- genes.by.pathway[[pathway]]
                               list.genes.in.pathway <- intersect(names(geneList), pathway.genes)
                               list.genes.not.in.pathway <- setdiff(names(geneList), list.genes.in.pathway)
                               scores.in.pathway <- geneList[list.genes.in.pathway]
                               scores.not.in.pathway <- geneList[list.genes.not.in.pathway]
                               if (length(scores.in.pathway) > 0){
                                 p.value <- wilcox.test(scores.in.pathway, scores.not.in.pathway, alternative = "less")$p.value
                               } else{
                                 p.value <- NA
                               }
                               return(c(p.value = p.value, Annotated = length(list.genes.in.pathway)))
                             }
))

# Assemble output table
outdat <- data.frame(pathway.code = rownames(pVals.by.pathway))
outdat$pathway.name <- pathways.list[outdat$pathway.code]
outdat$p.value <- pVals.by.pathway[,"p.value"]
outdat$Annotated <- pVals.by.pathway[,"Annotated"]
outdat <- outdat[order(outdat$p.value),] %>% mutate(adj_p.value = p.adjust(p.value, method = 'BH'))
View(outdat)

write.csv(outdat, file = "Enrichment/gs_enrichment_8762_2.csv")

####################################################
#######   Using Coverage and KEGG pathways    ######
####################################################

comparative_enrich <- function(X0, scores = list(), pathways.list, pathway.codes, genes.by.pathway){
  result <- list()
  for(i in 1:length(genes.by.pathway)){
    pathway_result <- data.frame(matrix(nrow = length(seq (5, 500, 1)), 
                                        ncol = length(scores)+2))
    colnames(pathway_result) <- append(names(scores), c("limma", "cutoff"))
    count <- 0
    for(score in scores){
      count <- count + 1
      count2 <- 0
      genes <- sort(score)
      control_genes <- sort(X0) 
      for(cutoff in seq (5, 500, 1)){
        count2 <- count2 + 1
        pathway_result[count2,count] <- sum(names(genes[1:cutoff]) %in% genes.by.pathway[[i]])/length(genes.by.pathway[[i]])*100
        pathway_result[count2,ncol(pathway_result)] <- cutoff
        
        if(count == 1){
          pathway_result[count2,ncol(pathway_result)-1] <- sum(names(control_genes[1:cutoff]) %in% genes.by.pathway[[i]])/length(genes.by.pathway[[i]])*100
        }
      }
    }
    result[[pathways.list[[i]]]] <- pathway_result
  }
  return(result)
}

pathways.list <- KEGGREST::keggList("pathway", "hsa")
pathway.codes <- sub("path:", "", names(pathways.list))
genes.by.pathway <- readRDS('genes_by_pathway.rds')

limma_geneList <- min_stats_new$p_des
names(limma_geneList) <- rownames(min_stats_new)
lean_geneList <- min_stats_new$p_lean
names(lean_geneList) <- rownames(min_stats_new)
aco_geneList <- min_stats_new$adj_p_fishers
names(aco_geneList) <- rownames(min_stats_new)
aco_limma_geneList <- min_stats_new$p_combs
names(aco_limma_geneList) <- rownames(min_stats_new)
gs_geneList <- vapply(gs_modules$p.Fisher, FUN = function(x){min(x*diam,1)}, FUN.VALUE = 1)
names(gs_geneList) <- gs_modules$gene.id

sum(rownames(min_stats_new) != gs_modules$gene.id) #should be 0
gs_limma_geneList <- pchisq(-2*(log(gs_geneList) + log(limma_geneList)), df = 4, lower.tail = FALSE)
names(gs_limma_geneList) <- gs_modules$gene.id

comparative_enrichment <- comparative_enrich(limma_geneList, list(ACo = aco_geneList,
                                                                  ACo_limma = aco_limma_geneList,
                                                                  GS = gs_geneList,
                                                                  GS_limma = gs_limma_geneList,
                                                                  lean = lean_geneList),
                                             pathways.list, pathway.codes, genes.by.pathway)

#Check which pathways to include
bc_pathways <- c(#"Pancreatic cancer - Homo sapiens (human)", 
                 #"Non-small cell lung cancer - Homo sapiens (human)", 
                 #"Acute myeloid leukemia - Homo sapiens (human)", 
                 #"Epstein-Barr virus infection - Homo sapiens (human)", 
                 #"Viral carcinogenesis - Homo sapiens (human)", 
                 #"Pathways in cancer - Homo sapiens (human)"
                 "Parkinson disease - Homo sapiens (human)",
                 "Pathways of neurodegeneration - multiple diseases - Homo sapiens (human)",
                 "Alzheimer disease - Homo sapiens (human)",
                 "Huntington disease - Homo sapiens (human)"
                 )
other_pathways <- c(#"FoxO signaling pathway - Homo sapiens (human)",
                    #"Hippo signaling pathway - Homo sapiens (human)", 
                    #"NF-kappa B signaling pathway - Homo sapiens (human)", #AML
                    #"JAK-STAT signaling pathway - Homo sapiens (human)", #AML
                    #"EGFR tyrosine kinase inhibitor resistance - Homo sapiens (human)", #AML, NSCLC, PC
                    #"Hedgehog signaling pathway - Homo sapiens (human)", #AML, PC
                    #"Notch signaling pathway - Homo sapiens (human)", #AML, PC
                    #"PI3K-Akt signaling pathway - Homo sapiens (human)", #AML
                    #"Progesterone-mediated oocyte maturation - Homo sapiens (human)",
                    #"Oocyte meiosis - Homo sapiens (human)",
                    #"Rap1 signaling pathway - Homo sapiens (human)",
                    #"Ras signaling pathway - Homo sapiens (human)", #PC
                    #"Sphingolipid signaling pathway - Homo sapiens (human)",
                    #"TGF-beta signaling pathway - Homo sapiens (human)", #AML, PC
                    "Cell cycle - Homo sapiens (human)", #p53
                    "p53 signaling pathway - Homo sapiens (human)", #HD, NSCLC, PC
                    "Cellular senescence - Homo sapiens (human)", #p53
                    "Apoptosis - Homo sapiens (human)", #AD, PD, HD, p53
                    "Wnt signaling pathway - Homo sapiens (human)", #AML, PC, AD
                    "Proteasome - Homo sapiens (human)", #PD, AD, HD
                    "Tyrosine metabolism - Homo sapiens (human)", #PD
                    "Oxidative phosphorylation - Homo sapiens (human)", #PD, AD, HD
                    "Dopaminergic synapse - Homo sapiens (human)", #PD
                    "Protein processing in endoplasmic reticulum - Homo sapiens (human)", #PD, AD, HD
                    "Autophagy - animal - Homo sapiens (human)", #AD, HD
                    "Mitophagy - animal - Homo sapiens (human)", #PD
                    "AGE-RAGE signaling pathway in diabetic complications - Homo sapiens (human)", #AD
                    "Calcium signaling pathway - Homo sapiens (human)", #PD, AD, HD
                    "Insulin signaling pathway - Homo sapiens (human)", #AD
                    "Ubiquitin mediated proteolysis - Homo sapiens (human)", #PD
                    "RNA polymerase - Homo sapiens (human)", #HD
                    "Basal transcription factors - Homo sapiens (human)", #HD
                    "Glutamatergic synapse - Homo sapiens (human)", #HD
                    "Endocytosis - Homo sapiens (human)" #HD
                    )

for(pathway in append(bc_pathways, other_pathways)){
  plot <- ggplot(data = comparative_enrichment[[pathway]]) +
    geom_line(aes(cutoff, limma), color = "black") +
    geom_line(aes(cutoff, ACo), color = "red") +
    geom_line(aes(cutoff, lean), color = "blue") +
    geom_line(aes(cutoff, GS), color = "magenta") +
    #geom_line(aes(cutoff, GS_limma), color = "cyan") +
    ylab("Cumulative Coverage (%)") +
    xlab("Cutoff of sorted gene list") +
    ggtitle(pathway) 
  #subtitle = "Black: Limma, Blue: ACo-Limma, Red: ACo, Magenta: GS, Cyan: GS-Limma"
  print(plot)
}

