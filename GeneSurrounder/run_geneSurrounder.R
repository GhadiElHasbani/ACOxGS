run_geneSurrounder <- function(adj.matrix,
                               gene_scores_per_sample,
                               class_labels,
                               nresamples,
                               seed,
                               num.Sphere.resamples = 1,
                               gene.id,
                               decay_only = TRUE,
                               file_name = "gs_results.csv",
                               cores = 2 # Set to 1 for Windows
                               ){
  ##preprocess
  print("Preprocessing...")
  gene_scores_per_sample <- gene_scores_per_sample[rownames(gene_scores_per_sample) %in% rownames(adj.matrix),]
   #generate t-stats
  tstats <- calcGeneTStats(gene_scores_per_sample, class_labels, nresamples, seed)
  geneStats.observed <- tstats$observed
  perm.geneStats.matrix <- tstats$resampled
   #generate correlation matrix
  cor.matrix <- calcCorMatrix(data.frame(gene_scores_per_sample), corMethod = "pearson", exprName = "gene_scores_per_sample", useMethod = "everything")
   #generate distance matrix
  adj.matrix <- adj.matrix[rownames(adj.matrix) %in% rownames(gene_scores_per_sample), ]
  adj.matrix <- adj.matrix[ ,colnames(adj.matrix) %in% rownames(gene_scores_per_sample)]
  int_network <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
  diameter <- diameter(int_network)
  print(diameter)
  distance.matrix <- calcAllPairsDistances(int_network, directionPaths = "all", networkName = "int_network")
  distance.matrix[!is.finite(distance.matrix)] <- -1
  genes.assayedETnetwork <- intersect(rownames(distance.matrix), rownames(cor.matrix))
  
  print("Running GeneSurrounder...")
  gs_results <- data.frame()
  if(!decay_only){
    gs_results <- parallel::mclapply(1:length(gene.id), function(i){
      start_time <- Sys.time()
      print(paste("Run", i))
      
      res <- geneNIDG(distance.matrix = distance.matrix,
                     cor.matrix = cor.matrix,
                     geneStats.observed = geneStats.observed,
                     perm.geneStats.matrix = perm.geneStats.matrix,
                     genes.assayedETnetwork = genes.assayedETnetwork,
                     diameter = diameter, # diameter >= 8 # diameter < 8 gives an error due to geneid.d line 376 of GeneSurrounder.R
                     num.Sphere.resamples = num.Sphere.resamples,
                     gene.id = gene.id[i],
                     seed = seed
                      )
      gs <- res$res
      
      print(gs[which.min(gs$p.Fisher),])
      gs_results <- rbind(gs_results, gs[which.min(gs$p.Fisher),])
      print(gs_results)
      end_time <- Sys.time()
      print(paste("Time:", end_time - start_time, ""))
      print(paste("p-Decay Time:", res$pdecay_time, ""))
      return(gs_results %>% mutate(time = end_time - start_time, pdecay_time = res$pdecay_time))
    }, mc.cores = cores)
    
    gs_results <- bind_rows(gs_results)
    write.csv(gs_results, file_name)
    return(gs_results)
  }else{
    gs_results <- parallel::mclapply(1:length(gene.id), function(i){
      start_time <- Sys.time()
      print(paste("Run", i))
      distances <- distance.matrix[gene.id[i],
                                   genes.assayedETnetwork]
      
      
      sizes <- vapply(1:diameter,function(RADIUS){
        
        igenes.distances <- distances[distances <= RADIUS
                                      & distances > 0]
        
        length(igenes.distances)
        
        
      },
      numeric(1))
      observed.tau_b <- Observed.DecayDE(distance.matrix,
                                         gene.id[i],
                                         genes.assayedETnetwork,
                                         diameter,
                                         geneStats.observed)
      
      
      null.tau_b <- Resample.DecayDE(distance.matrix,
                                     gene.id[i],
                                     genes.assayedETnetwork,
                                     diameter,
                                     perm.geneStats.matrix,
                                     sizes)
      
      
      num.Decay.resamples <- nrow(perm.geneStats.matrix)
      
      # proportion of null taub \LEQ observed i.e more discordant
      p.Decay <- vapply(1:diameter,function(distance){
        
        observed.p <- observed.tau_b[distance]
        
        null.p <- null.tau_b[distance,1:num.Decay.resamples]
        
        return(length(null.p[null.p <= observed.p])/length(null.p))
        
      },
      numeric(1))
      
      p.Decay[p.Decay == 0] <- 1/(num.Decay.resamples+1)
      
      gs <- data.frame(gene.id = rep(gene.id[i],diameter),
                       radius = 1:diameter,
                       size = sizes,
                       observed.tau_b = observed.tau_b,
                       p.Decay = p.Decay)
      
      print(gs[which.min(gs$p.Decay),])
      gs_results <- gs[which.min(gs$p.Decay),]
      #gs_results <- rbind(gs_results, gs[which.min(gs$p.Decay),])
      end_time <- Sys.time()
      #time[i] <- end_time - start_time
      gs_results$time <- end_time - start_time
      print(paste("Time:", time[i], ""))
      return(gs_results)
    }, mc.cores = cores)
    #gs_results <- gs_results %>% mutate(time = time)
    gs_results <- bind_rows(gs_results)
    write.csv(gs_results, file_name)
    return(gs_results)
  }
}
