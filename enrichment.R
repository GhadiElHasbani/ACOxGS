source("ant_colony.R")

library(dplyr)
library(readr)
library(tibble)
library(ggplot2)
library(limma) # DEA tool / GS dependency
library(LEANR) # LEAN package
library(doMC) # LEAN dependency
library(annotate) #map NCBI Gene IDs to gene symbols
library(GSAR) # p53 data
library(GSEABenchmarkeR) # GSEABenchmark data
library(org.Hs.eg.db)
library(KEGGREST)

setwd("~/Desktop/LAU/Capstone/ACOxGS/")

get_processed_data <- function(geo2kegg, dataset_id, seed, n_ants = 40, trial = 'gamma', adj_method = 'BH', cutoff = 0.05, n_resamples = 100){
  
  ids <- c()
  sizes <- c()
  for(id in names(geo2kegg)){
    ids<- c(ids, id)
    idx <- which(names(geo2kegg) == id)
    eset <- geo2kegg[[idx]]
    sizes <- c(sizes, dim(eset)[2])
  }
  bm_ds_ids <- data.frame(ids, sizes) %>% dplyr::arrange(desc(sizes))
  
  dataset <- paste("GSE", dataset_id, sep = "")
  dataset_id <- sub("_", "", dataset_id)
  filter_by <- 'ADJ.PVAL_filter'
  idx <- which(names(geo2kegg) == dataset)
  aco_modules <- readRDS(paste(paste('Objects/Run ', seed,'/aco', sep=''), trial, 'trial', dataset_id, n_ants, paste(seed,'rds', sep = '.'), sep = '_'))
  #aco_modules <- readRDS('aco_gamma_trial_p53_40.rds')
  gs_modules <- readRDS(paste(paste('Objects/Run ', seed,'/gs', sep=''), 'results', dataset_id, paste(seed, 'rds', sep = '.'), sep = '_'))
  
  # TP53 
  #data("p53DataSet")
  
  #df <- p53DataSet
  #class_labels <- c(rep("0", 17), rep("1", 33))
  
  #GSEABenchmarkeR
  eset <- geo2kegg[[idx]]
  df <- eset@assays@data@listData$exprs
  class_labels <- as.character(eset@colData@listData$GROUP)
  uid_to_symbol <- lookUp(eset@NAMES, 'org.Hs.eg', 'SYMBOL')
  gene_symbols <- unlist(uid_to_symbol)
  rownames(df)<- gene_symbols
  
  lcc <- as.matrix(read.csv(paste('Networks/KEGG', 'A', 'undir', paste('GE', dataset_id, sep = ''), 'LCC.csv', sep = '_'), row.names = 1))
  rownames(lcc)[!rownames(lcc) %in% colnames(lcc)] <- sub('-', '.', rownames(lcc)[!rownames(lcc) %in% colnames(lcc)])
  adj_mat <- lcc
  KEGG_A_LCC <- lcc
  adj_mat <- adj_mat[rownames(adj_mat) %in% rownames(df), ]
  adj_mat <- adj_mat[ ,colnames(adj_mat) %in% rownames(df)]
  int_network <- igraph::graph_from_adjacency_matrix(adj_mat, mode = 'undirected')
  diam <- diameter(int_network)
  
  limma_res <- acoCalcGeneTStats(df, class_labels, n_resamples, seed, permutate = FALSE)
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
    
    return(LEANR::run.lean(gene_stats, int_network, verbose = FALSE, ranked = FALSE, n_reps = n_resamples, ncores = n_cores))
  }
  
  lean_modules <- runLEAN(df, KEGG_A_LCC, n_cores = NULL, seed = seed, 
                          n_resamples = n_resamples, class_labels = class_labels)
  #lean_modules <- LEANR::run.lean((min_stats_new %>% arrange(p_combs) %>% dplyr::select(p_combs))[,1], int_network, verbose = TRUE, ranked = FALSE, n_reps = 100, ncores = NULL)
  lean <- lean_modules$restab[rownames(lean_modules$restab) %in% names(aco_modules),8]
  lean <- lean[order(factor(names(lean), levels = names(aco_modules)))]
  
  #apply_comb_to <- 'adj_p_fishers' #'adj_p_fishers'
  #k <- 1
  #n_iter <- 5
  min_stats_new <- min_stats %>% mutate(FC = apply(aco_df[,!as.logical(as.numeric(class_labels))], 1, mean) - apply(aco_df[,as.logical(as.numeric(class_labels))], 1, mean),
                                        p_des = de_pvals,
                                        p_lean = lean,
                                        #fisher_stat = -2*(log(!!as.name(apply_comb_to)) + log(p_des))
                                        #p_fishers_log = -log10(p_fishers), 
                                        #bonf_adj = p_fishers_log >= -log10(cutoff/n_iter),
                                        #p_combs = pchisq(-2*(log(!!as.name(apply_comb_to)) + log(p_des)), df = 4, lower.tail = FALSE)
                                        #BH_adj_p = p.adjust(!!as.name(apply_BH_to), method = "BH", n = length(!!as.name(apply_BH_to))*k) #which adjustment to use?
  )
  return(list(min_stats_new, gs_modules, diam))
}


####################################################
#######   Using Wilcoxon Rank and KEGG pathways ####
####################################################

run_enrichment <- function(min_stats_new, gs_modules, diam, dataset_id, seed, n_ants = 40, trial = 'gamma', adj_method = 'BH'){
  
  pathways.list <- KEGGREST::keggList("pathway", "hsa")
  head(pathways.list)
  pathway.codes <- sub("path:", "", names(pathways.list)) 
  genes.by.pathway <- readRDS('Objects/genes_by_pathway.rds')
  #genes.by.pathway <- sapply(pathway.codes,
  #                           function(pwid){
  #                             pw <- KEGGREST::keggGet(pwid)
  #                             if (is.null(pw[[1]]$GENE)) return(NA)
  #                             pw2 <- pw[[1]]$GENE[c(FALSE,TRUE)] # may need to modify this to c(FALSE, TRUE) for other organisms
  #                             pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
  #                             return(pw2)
  #                           }
  #)
  #saveRDS(genes.by.pathway, 'Objects/genes_by_pathway.rds')
  
  tools <- c(paste("aco", trial, n_ants, sep = "_"), "gs", "lean", "limma")
  
  for(n_tool in 1:4){
    tool <- tools[n_tool]
    
    if(n_tool == 1){
      geneList <- min_stats_new$adj_p_fishers
      names(geneList) <- rownames(min_stats_new)
    }else if(n_tool == 2){
      geneList <- vapply(gs_modules$p.Fisher, FUN = function(x){min(x*diam,1)}, FUN.VALUE = 1)
      names(geneList) <- gs_modules$gene.id
    }else if(n_tool == 3){
      geneList <- min_stats_new$p_lean
      names(geneList) <- rownames(min_stats_new)
    }else{
      geneList <- min_stats_new$p_des
      names(geneList) <- rownames(min_stats_new)
    }
    
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
    outdat <- outdat[order(outdat$p.value),] %>% mutate(adj_p.value = p.adjust(p.value, method = adj_method))
    #View(outdat)
    write.csv(outdat, file = paste(paste("Enrichment", dataset_id, strsplit(tool, split = "_")[[1]][1], tool, sep='/'), "_enrichment_",  dataset_id, "_", seed, ".csv", sep = ''))
  }
}

#Example
geo2kegg <- loadEData("geo2kegg", preproc = TRUE)

dataset_id <- "20291"
seeds <- 1:10

for(seed in seeds){
  print(paste("================================= Run", seed, "========================================="))
  res <- get_processed_data(geo2kegg, dataset_id, seed)
  min_stats_new <- res[[1]]
  gs_modules <- res[[2]]
  diam <- res[[3]]
  run_enrichment(min_stats_new, gs_modules, diam, sub("_", "", dataset_id), seed)
  print("done")
  print("=================================================================================")
}

summarize_enrich <- function(dataset_id, trial = 'gamma', n_ants = 40, cutoff = 0.05, seeds = 1:6){
  for(tool in c(paste("aco", trial, n_ants, sep = "_"), "gs", "lean", "limma")){
    pattern <- ".*_["
    for(seed in seeds){
      pattern <- paste(pattern, seed, ",", sep = '')
    }
    dfs <- list.files(path=paste("Enrichment", dataset_id, strsplit(tool, split = "_")[[1]][1], sep = '/'), pattern = paste(tool, substr(temp, 1, nchar(temp) - 1), "].csv", sep = ''), full.names = TRUE) %>% 
      lapply(read_csv, col_types = cols(col_skip())) %>%
      lapply(tibble::rownames_to_column, var = "rank") %>%
      lapply(mutate, rank = as.integer(rank))
    
    #number of significant pathways
    n_sig_pathways <- dfs %>% lapply(summarise, n_sig_pathways = length(rank[adj_p.value <= cutoff])) %>% 
      bind_rows %>% 
      summarise(mean_n_sig_pathways = mean(n_sig_pathways), sd_n_sig_pathways = sd(n_sig_pathways))
    
    #organizing relevant pathways
    common_pathways <- c("Pathways of neurodegeneration - multiple diseases - Homo sapiens (human)",
                         "Apoptosis - Homo sapiens (human)",
                         "Proteasome - Homo sapiens (human)",
                         "Oxidative phosphorylation - Homo sapiens (human)",
                         "Protein processing in endoplasmic reticulum - Homo sapiens (human)",
                         "Calcium signaling pathway - Homo sapiens (human)")
    if(dataset_id == '20291'){
      main_pathway <- "Parkinson disease - Homo sapiens (human)"
      specific_pathways <- c("Tyrosine metabolism - Homo sapiens (human)",
                          "Dopaminergic synapse - Homo sapiens (human)",
                          "Mitophagy - animal - Homo sapiens (human)",
                          "Ubiquitin mediated proteolysis - Homo sapiens (human)")
    }else if(dataset_id == '5281VCX'){
      main_pathway <- "Alzheimer disease - Homo sapiens (human)"
      specific_pathways <- c("Autophagy - animal - Homo sapiens (human)",
                          "Wnt signaling pathway - Homo sapiens (human)",
                          "AGE-RAGE signaling pathway in diabetic complications - Homo sapiens (human)",
                          "Insulin signaling pathway - Homo sapiens (human)")
    }else if(dataset_id == '8762'){
      main_pathway <- "Huntington disease - Homo sapiens (human)"
      specific_pathways <- c("p53 signaling pathway - Homo sapiens (human)",
                          "Autophagy - animal - Homo sapiens (human)",
                          "RNA polymerase - Homo sapiens (human)",
                          "Basal transcription factors - Homo sapiens (human)",
                          "Glutamatergic synapse - Homo sapiens (human)",
                          "Endocytosis - Homo sapiens (human)")
    }
    relevant_pathways <- c(main_pathway, common_pathways, specific_pathways)
    
    #number of signficant relevant pathways
    sig_relevant_pathways <- dfs %>% 
      lapply(summarise, n_sig_relevant_paths = sum(relevant_pathways %in% pathway.name[adj_p.value <= cutoff])) %>% 
      bind_rows %>% 
      summarise(mean_n_sig_relevant_paths = mean(n_sig_relevant_paths), sd_n_sig_relevant_paths = sd(n_sig_relevant_paths))
    
    #rank of disease pathway
    rank_main_path <- dfs %>% 
      lapply(summarise, rank_disease_path = rank[pathway.name == main_pathway]) %>% 
      bind_rows %>%
      summarise(mean_rank_disease_path = mean(rank_disease_path), sd_rank_disease_path = sd(rank_disease_path))
    
    #write results to summary file
    file_path <- paste(paste("Enrichment", dataset_id, strsplit(tool, split = "_")[[1]][1], sep = '/'), "/", tool, "_enrichment_summary.txt", sep = '')
    write("Rank of disease pathway:\n", file_path, append = TRUE)
    write.table(rank_main_path, file = file_path, row.names = FALSE, append = TRUE)
    write("\n=============================================\n", file_path, append = TRUE)
    write("Number of significant pathways:\n", file_path, append = TRUE)
    write.table(n_sig_pathways, file = file_path, row.names = FALSE, append = TRUE)
    write("\n=============================================\n", file_path, append = TRUE)
    write("Number of signficant relevant pathways:\n", file_path, append = TRUE)
    write.table(sig_relevant_pathways, file = file_path, row.names = FALSE, append = TRUE)
  }
}

summarize_enrich("20291")
summarize_enrich("5281VCX")

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

#geo2kegg <- loadEData("geo2kegg", preproc = TRUE)
#seed <- 3
#dataset_id <- "20291"
#res <- get_processed_data(geo2kegg, dataset_id, seed)
#min_stats_new <- res[[1]]
#gs_modules <- res[[2]]
#diam <- res[[3]]

plot_comparative_enrich <- function(min_stats_new, gs_modules, diam, dataset_id, seed){
  
                                      
  pathways.list <- KEGGREST::keggList("pathway", "hsa")
  pathway.codes <- sub("path:", "", names(pathways.list))
  genes.by.pathway <- readRDS('genes_by_pathway.rds')
  
  limma_geneList <- min_stats_new$p_des
  names(limma_geneList) <- rownames(min_stats_new)
  lean_geneList <- min_stats_new$p_lean
  names(lean_geneList) <- rownames(min_stats_new)
  aco_geneList <- min_stats_new$adj_p_fishers
  names(aco_geneList) <- rownames(min_stats_new)
  gs_geneList <- vapply(gs_modules$p.Fisher, FUN = function(x){min(x*diam,1)}, FUN.VALUE = 1)
  names(gs_geneList) <- gs_modules$gene.id
  
  sum(rownames(min_stats_new) != gs_modules$gene.id) #should be 0
  
  comparative_enrichment <- comparative_enrich(limma_geneList, list(ACo = aco_geneList,
                                                                    GS = gs_geneList,
                                                                    lean = lean_geneList),
                                               pathways.list, pathway.codes, genes.by.pathway)
  
  #Check which pathways to include
  main_pathways <- c("Parkinson disease - Homo sapiens (human)",
                   "Pathways of neurodegeneration - multiple diseases - Homo sapiens (human)",
                   "Alzheimer disease - Homo sapiens (human)",
                   "Huntington disease - Homo sapiens (human)"
                   )
  other_pathways <- c("Cell cycle - Homo sapiens (human)", #p53
                      "p53 signaling pathway - Homo sapiens (human)", #HD
                      "Cellular senescence - Homo sapiens (human)", #p53
                      "Apoptosis - Homo sapiens (human)", #AD, PD, HD, p53
                      "Wnt signaling pathway - Homo sapiens (human)", #AD
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
  
  for(pathway in append(main_pathways, other_pathways)){
    plot <- ggplot(data = comparative_enrichment[[pathway]]) +
      geom_line(aes(cutoff, limma), color = "black") +
      geom_line(aes(cutoff, ACo), color = "red") +
      geom_line(aes(cutoff, lean), color = "blue") +
      geom_line(aes(cutoff, GS), color = "magenta") +
      ylab("Cumulative Coverage (%)") +
      xlab("Cutoff of sorted gene list") +
      ggtitle(pathway) 
    #subtitle = "Black: Limma, Blue: LEAN, Red: ACo, Magenta: GS"
    print(plot)
  }
}

