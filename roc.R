library(pROC)
library(dplyr)
library(ggplot2)
library(limma) # DEA tool / GS dependency
library(LEANR) # LEAN package
library(doMC) # LEAN dependency
library(annotate) #map NCBI Gene IDs to gene symbols
library(GSAR) # p53 data
library(GSEABenchmarkeR) # GSEABenchmark data
library(org.Hs.eg.db)

generate_ROC <- function(run, dataset, disease, variable_cutoff = FALSE){
    aco <- read.csv(paste("Enrichment/", dataset, "/aco/aco_gamma_40_enrichment_", dataset, "_", run,".csv", sep = "")) %>% arrange(pathway.name)
    gs <- read.csv(paste("Enrichment/", dataset, "/gs/gs_enrichment_", dataset, "_", run, ".csv", sep = "")) %>% arrange(pathway.name)
    lean <- read.csv(paste("Enrichment/", dataset, "/lean/lean_enrichment_", dataset, "_", run, ".csv", sep = "")) %>% arrange(pathway.name)
    limma <- read.csv(paste("Enrichment/", dataset, "/limma/limma_enrichment_", dataset, "_", run, ".csv", sep = "")) %>% arrange(pathway.name)
    
    relevant_paths_HD <- c("RNA polymerase - Homo sapiens (human)", #HD
    "Basal transcription factors - Homo sapiens (human)", #HD
    "Glutamatergic synapse - Homo sapiens (human)", #HD
    "Endocytosis - Homo sapiens (human)", #HD
    "Calcium signaling pathway - Homo sapiens (human)",
    "Protein processing in endoplasmic reticulum - Homo sapiens (human)", #PD, AD, HD
    "Autophagy - animal - Homo sapiens (human)", #AD, HD
    "Oxidative phosphorylation - Homo sapiens (human)", #PD, AD, HD
    "Proteasome - Homo sapiens (human)", #PD, AD, HD
    "Apoptosis - Homo sapiens (human)", #AD, PD, HD, p53
    "p53 signaling pathway - Homo sapiens (human)", #HD
    "Huntington disease - Homo sapiens (human)",
    "Pathways of neurodegeneration - multiple diseases - Homo sapiens (human)")
    
    relevant_paths_PD <- c("Apoptosis - Homo sapiens (human)", #AD, PD, HD, p53
    "Proteasome - Homo sapiens (human)", #PD, AD, HD
    "Tyrosine metabolism - Homo sapiens (human)", #PD
    "Oxidative phosphorylation - Homo sapiens (human)", #PD, AD, HD
    "Dopaminergic synapse - Homo sapiens (human)", #PD
    "Protein processing in endoplasmic reticulum - Homo sapiens (human)", #PD, AD, HD
    "Mitophagy - animal - Homo sapiens (human)", #PD
    "Calcium signaling pathway - Homo sapiens (human)", #PD, AD, HD
    "Ubiquitin mediated proteolysis - Homo sapiens (human)", #PD
    "Parkinson disease - Homo sapiens (human)",
    "Pathways of neurodegeneration - multiple diseases - Homo sapiens (human)")
    
    relevant_paths_AD <- c("Apoptosis - Homo sapiens (human)", #AD, PD, HD, p53
    "Proteasome - Homo sapiens (human)", #PD, AD, HD
    "Oxidative phosphorylation - Homo sapiens (human)", #PD, AD, HD
    "Protein processing in endoplasmic reticulum - Homo sapiens (human)", #PD, AD, HD
    "Calcium signaling pathway - Homo sapiens (human)", #PD, AD, HD
    "Autophagy - animal - Homo sapiens (human)", #AD, HD
    "Alzheimer disease - Homo sapiens (human)",
    "Insulin signaling pathway - Homo sapiens (human)", #AD
    "Wnt signaling pathway - Homo sapiens (human)", #AD
    "AGE-RAGE signaling pathway in diabetic complications - Homo sapiens (human)", #AD
    "Pathways of neurodegeneration - multiple diseases - Homo sapiens (human)")
    
    all_paths <- aco[,'pathway.name']
    
    get_y <- function(p_vals, alpha = 0.05){
        return(ifelse(is.na(p_vals), 0, as.numeric(p_vals <= alpha)))
    }
    
    if(disease == "PD"){
        true_y <- ifelse(all_paths %in% relevant_paths_PD, 1, 0)
    }else if(disease == "HD"){
        true_y <- ifelse(all_paths %in% relevant_paths_HD, 1, 0)
    }else{
        true_y <- ifelse(all_paths %in% relevant_paths_AD, 1, 0)
    }
    #resulting NA values have same index
    true_y <- true_y[!is.na(aco[,'adj_p.value'])]
    
    if(!variable_cutoff){ #also makes possible to combine all in one plot
        aco_y <- get_y(aco[,'adj_p.value'][!is.na(aco[,'adj_p.value'])])
        gs_y <- get_y(gs[,'adj_p.value'][!is.na(gs[,'adj_p.value'])])
        lean_y <- get_y(lean[,'adj_p.value'][!is.na(lean[,'adj_p.value'])])
        limma_y <- get_y(limma[,'adj_p.value'][!is.na(limma[,'adj_p.value'])])
        
        
        #plot_ROC <- function(y, true_y){
        #  return(roc(true_y,y,
        #                  smoothed = TRUE,
        #                  # arguments for plot
        #                  plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
        #                  print.auc=TRUE))
        #}
        
        #aco_roc <- plot_ROC(aco_y, true_y)
        #gs_roc <- plot_ROC(gs_y, true_y)
        #lean_roc <- plot_ROC(lean_y, true_y)
        #limma_roc <- plot_ROC(limma_y, true_y)
        aco_roc <- roc(aco_y, true_y)
        gs_roc <- roc(gs_y, true_y)
        lean_roc <- roc(lean_y, true_y)
        tryCatch(
            {
                limma_roc <- roc(limma_y, true_y)
                
                ggroc(list(ACO = aco_roc, GS = gs_roc, LEAN = lean_roc, limma = limma_roc))
            },
            error = function(cond) {
                message("limma produced no results for this experiment with this dataset. Excluding...")
                
                ggroc(list(ACO = aco_roc, GS = gs_roc, LEAN = lean_roc))
            }
        )
        
        
    }else{
        #### for variable cutoffs
        
        library(ROCR)
        
        aco_y <- aco[,'adj_p.value'][!is.na(aco[,'adj_p.value'])]
        gs_y <- gs[,'adj_p.value'][!is.na(gs[,'adj_p.value'])]
        lean_y <- lean[,'adj_p.value'][!is.na(lean[,'adj_p.value'])]
        limma_y <- limma[,'adj_p.value'][!is.na(limma[,'adj_p.value'])]
        
        
        aco_roc <- performance(prediction(aco_y, true_y),"tpr","fpr")
        gs_roc <- performance(prediction(gs_y, true_y),"tpr","fpr")
        lean_roc <- performance(prediction(lean_y, true_y),"tpr","fpr")
        
        plot(aco_roc)
        abline(a = 0, b = 1)
        
        plot(gs_roc)
        abline(a = 0, b = 1)
        
        plot(lean_roc)
        abline(a = 0, b = 1)
        
        tryCatch(
            {
                limma_roc <- performance(prediction(limma_y, true_y),"tpr","fpr")
                
                plot(limma_roc)
                abline(a = 0, b = 1)
            },
            error = function(cond) {
                message("limma produced no results for this experiment with this dataset. Excluding...")
            }
        )
        
        
    }
    
}

#Example 

#make sure to set the directory to be the cloned repository. something like: "~/Desktop/ACOxGS" with appropriate windows format
setwd(.......)

run <- 3 # can be any of 1 to 10 (int)
dataset <- '8762' #can be any of 20291, 5281VCX, 8762 (need to be matched to appropriate disease argument)
disease <- "HD" #can be any of "PD", "AD", "HD" in this order corresponding to dataset ID above

generate_ROC(run, dataset, disease)
