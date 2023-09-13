
########## Alzheimer's


# importing data

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

df1 <- eset1@assays@data@listData$exprs[common_genes,]
class_labels1 <- as.character(eset1@colData@listData$GROUP)
gene_symbols <- unlist(lookUp(eset1@NAMES, 'org.Hs.eg', 'SYMBOL'))
sum(rownames(df1) != names(gene_symbols))
rownames(df1)<- gene_symbols

df2 <- eset2@assays@data@listData$exprs[common_genes,]
class_labels2 <- as.character(eset2@colData@listData$GROUP)
gene_symbols <- unlist(lookUp(eset2@NAMES, 'org.Hs.eg', 'SYMBOL'))
sum(rownames(df2) != names(gene_symbols))
rownames(df2)<- gene_symbols

df3 <- eset3@assays@data@listData$exprs[common_genes,]
class_labels3 <- as.character(eset3@colData@listData$GROUP)
gene_symbols <- unlist(lookUp(eset3@NAMES, 'org.Hs.eg', 'SYMBOL'))
sum(rownames(df3) != names(gene_symbols))
rownames(df3)<- gene_symbols

lcc <- as.matrix(read.csv('./KEGG_A_undir_GE5281VCX_LCC.csv', row.names = 1))
rownames(lcc)[!rownames(lcc) %in% colnames(lcc)] <- sub('-', '.', rownames(lcc)[!rownames(lcc) %in% colnames(lcc)])
adj_mat <- lcc

# importing and processing results

aco_modules1 <- readRDS('aco_gamma_trial_5281VCX_40_3.rds')
aco_modules2 <- readRDS('aco_gamma_trial_5281HIP_40_3.rds')
aco_modules3 <- readRDS('aco_gamma_trial_5281EC_40_3.rds')

gs_modules1 <- readRDS('gs_results_5281VCX_3.rds')
gs_modules2 <- readRDS('gs_results_5281HIP_3.rds')
gs_modules3 <- readRDS('gs_results_5281EC_3.rds')

calcLimmaStats <- function(df, class_labels, n_resamples, seed){
  limma_res <- acoCalcGeneTStats(df, class_labels, n_resamples, seed)
  hist(limma_res$complete_res$adj.P.Val)
  limma_adj_pvals <- limma_res$complete_res$adj.P.Val
  names(limma_adj_pvals) <- rownames(limma_res$complete_res)
  return(limma_adj_pvals)
}

limma_adj_pvals1 <- calcLimmaStats(df1, as.factor(class_labels1), 100, 3)
limma_adj_pvals2 <- calcLimmaStats(df2, class_labels2, 100, 3)
limma_adj_pvals3 <- calcLimmaStats(df3, class_labels3, 100, 3)

processAcoStats <- function(aco_modules, df, class_labels, KEGG_A_LCC, limma_adj_pvals, get_min_by, cutoff, adj_method, seed){
  
  all_stats <- data.frame()
  min_stats <- data.frame()
  tot_time <- c()
  tot_time_par <- c()
  min_oscc <- c()
  min_pdecay <- c()
  de_pvals <- c()
  modules <- list()
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
  rownames(min_stats) <- names(aco_modules)
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
  lean <- lean_modules$restab[rownames(lean_modules$restab) %in% names(aco_modules),8]
  lean <- lean[order(factor(names(lean), levels = names(aco_modules)))]
  
  min_stats_new <- min_stats %>% mutate(p_fishers_log = -log10(p_fishers), 
                                        bonf_adj = p_fishers_log >= -log10(cutoff/5),
                                        p_fishers_bonf = vapply(p_fishers*5, FUN = min, 1),
                                        p_lean = lean,
                                        p_des = de_pvals,
                                        p_combs = pchisq(-2*(log(p_fishers) + log(p_des)), df = 4, lower.tail = FALSE),#pchisq(-2*(log(p_spheres) + log(p_decays) + log(p_des)), df = 6, lower.tail = FALSE),
                                        bonf_adj_combs = p_combs <= cutoff)
  #p_fishers_BH = p.adjust(p_fishers, method = "BH"),
  #p_fishers_bonf = p.adjust(p_fishers, method = "bonferroni"),
  #bonf_adj = p_fishers_bonf <= 0.05,
  #BH_adj = p_fishers_BH <= 0.05)
  
  sum(min_stats_new$bonf_adj_combs)
  return(min_stats_new)
}

min_stats_new1 <- processAcoStats(aco_modules1, df1, class_labels1, adj_mat, limma_adj_pvals1, 'p_decays', 0.05, 'BH', 3)
min_stats_new2 <- processAcoStats(aco_modules2, df2, class_labels2, adj_mat, limma_adj_pvals2, 'p_decays', 0.05, 'BH', 3)
min_stats_new3 <- processAcoStats(aco_modules3, df3, class_labels3, adj_mat, limma_adj_pvals3, 'p_decays', 0.05, 'BH', 3)
min_stats_new1 <- min_stats_new1[rownames(min_stats_new1) %in% rownames(min_stats_new2),]
min_stats_new2<- min_stats_new2[rownames(min_stats_new2) %in% rownames(min_stats_new1),]
min_stats_new3<- min_stats_new3[rownames(min_stats_new3) %in% rownames(min_stats_new1),]

adj_mat <- adj_mat[rownames(adj_mat) %in% names(aco_modules1), ]
adj_mat <- adj_mat[ ,colnames(adj_mat) %in% names(aco_modules1)]
int_network <- igraph::graph_from_adjacency_matrix(adj_mat, mode = 'undirected')
diam <- diameter(int_network)

gs_modules1 <-gs_modules1 %>% filter(gene.id %in% rownames(min_stats_new1))
gs_modules2 <-gs_modules2 %>% filter(gene.id %in% rownames(min_stats_new1))
gs_modules3 <-gs_modules3 %>% filter(gene.id %in% rownames(min_stats_new1))

# cross-study concordance

corr_aco_1_2 <- cor.test(x=min_stats_new1$p_fishers, y=min_stats_new2$p_fishers, method = 'spearman')
corr_aco_1_3 <- cor.test(x=min_stats_new1$p_fishers, y=min_stats_new3$p_fishers, method = 'spearman')
corr_aco_2_3 <- cor.test(x=min_stats_new2$p_fishers, y=min_stats_new3$p_fishers, method = 'spearman')

corr_gs_1_2 <- cor.test(x=gs_modules1$p.Fisher, y=gs_modules2$p.Fisher, method = 'spearman')
corr_gs_1_3 <- cor.test(x=gs_modules1$p.Fisher, y=gs_modules3$p.Fisher, method = 'spearman')
corr_gs_2_3 <- cor.test(x=gs_modules2$p.Fisher, y=gs_modules3$p.Fisher, method = 'spearman')

corr_lean_1_2 <- cor.test(x=min_stats_new1$p_lean, y=min_stats_new2$p_lean, method = 'spearman')
corr_lean_1_3 <- cor.test(x=min_stats_new1$p_lean, y=min_stats_new3$p_lean, method = 'spearman')
corr_lean_2_3 <- cor.test(x=min_stats_new2$p_lean, y=min_stats_new3$p_lean, method = 'spearman')

corr_limma_1_2 <- cor.test(x=min_stats_new1$p_des, y=min_stats_new2$p_des, method = 'spearman')
corr_limma_1_3 <- cor.test(x=min_stats_new1$p_des, y=min_stats_new3$p_des, method = 'spearman')
corr_limma_2_3 <- cor.test(x=min_stats_new2$p_des, y=min_stats_new3$p_des, method = 'spearman')

aco_concordance <- matrix(c(1,corr_aco_1_2$estimate, corr_aco_1_3$estimate, corr_aco_1_2$estimate, 1, corr_aco_2_3$estimate, corr_aco_1_3$estimate, corr_aco_2_3$estimate, 1), 
                          nrow = 3, ncol = 3, byrow = TRUE, 
                          dimnames = list(c("GSE5281_VCX", "GSE5281_HIP", "GSE5281_EC"), c("GSE5281_VCX", "GSE5281_HIP", "GSE5281_EC")))

gs_concordance <- matrix(c(1,corr_gs_1_2$estimate, corr_gs_1_3$estimate, corr_gs_1_2$estimate, 1, corr_gs_2_3$estimate, corr_gs_1_3$estimate, corr_gs_2_3$estimate, 1), 
                          nrow = 3, ncol = 3, byrow = TRUE, 
                          dimnames = list(c("GSE5281_VCX", "GSE5281_HIP", "GSE5281_EC"), c("GSE5281_VCX", "GSE5281_HIP", "GSE5281_EC")))

lean_concordance <- matrix(c(1,corr_lean_1_2$estimate, corr_lean_1_3$estimate, corr_lean_1_2$estimate, 1, corr_lean_2_3$estimate, corr_lean_1_3$estimate, corr_lean_2_3$estimate, 1), 
                         nrow = 3, ncol = 3, byrow = TRUE, 
                         dimnames = list(c("GSE5281_VCX", "GSE5281_HIP", "GSE5281_EC"), c("GSE5281_VCX", "GSE5281_HIP", "GSE5281_EC")))

limma_concordance <- matrix(c(1,corr_limma_1_2$estimate, corr_limma_1_3$estimate, corr_limma_1_2$estimate, 1, corr_limma_2_3$estimate, corr_limma_1_3$estimate, corr_limma_2_3$estimate, 1), 
                          nrow = 3, ncol = 3, byrow = TRUE, 
                          dimnames = list(c("GSE5281_VCX", "GSE5281_HIP", "GSE5281_EC"), c("GSE5281_VCX", "GSE5281_HIP", "GSE5281_EC")))

write.csv(aco_concordance, "Concordance/aco_gamma_40_5281_concordance_3.csv")
write.csv(gs_concordance, "Concordance/gs_5281_concordance_3.csv")
write.csv(lean_concordance, "Concordance/lean_5281_concordance_3.csv")
write.csv(limma_concordance, "Concordance/limma_5281_concordance_3.csv")