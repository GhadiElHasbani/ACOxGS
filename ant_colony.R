#install.packages("tidyverse")
#install.packages("mnd")
#setwd("C:/Users/Frederick/Desktop/GitHub_Repos/GS_optimization") #Fred
#setwd("~/Desktop/GeneSurrounder_optimization/GS_optimization/")#Ghadi

# Define Global vars
#normalized_distance <- 
#distance <- 
#gene_scores <-
#adj_mat <- 

#capacity <- normalize(log fold change)

#update_capacity <- function(cap, position, route) {
#  el <- get_energy_lost
#}

#correlation coeff + 1 

get_energy_lost <- function(possible_position, route, oscc_old, 
                            gene_scores, normalized_distance, distance, 
                            current_capacity, penalize_dist = TRUE) {
  
  #possible_position -> possibile neighbor being investigated
  #route -> list of previous gene ids traversed
  #route = [BRCA1, ...]
  
  current <- gene_scores[possible_position, 1]
  previous <- gene_scores[route, 1]
  current_distance <- distance[possible_position]
  previous_distance <- distance[route]
  
  current_info <- data.frame(score = current, dist =  current_distance)
  rownames(current_info) <- possible_position
  
  previous_info <- data.frame(score = previous, dist =  previous_distance)
  rownames(previous_info) <- route
  
  info <- rbind(current_info, previous_info)
  oscc_new <- calculate_OSCC(abs(info$score), info$dist)
  delta_oscc <- oscc_new - oscc_old
  normalized_delta_oscc <- (delta_oscc + 2)/4
  if(penalize_dist){
    el <- normalized_delta_oscc*normalized_distance[possible_position]
  }else{
    el <- normalized_delta_oscc
  }
  return(data.frame(oscc_current = oscc_new, 
              delta_oscc = normalized_delta_oscc,
              EL = el, 
              new_capacity = current_capacity - el, 
              favorability = 1 - el))
  
  #normalized distance/delta_pdecay  ->  [0,1]/[0,1]
  
  #normalized distance * p decay
}

move_ant <- function(current_position, route, oscc_old, 
                     gene_scores, normalized_distance, distance, adj_mat,
                     current_capacity, alpha, beta, pheromone_map, penalize_dist = TRUE, first = FALSE, remove_close = TRUE, verbose = 1, trial = 1, visited = c()){
  require(purrr)
  
  get_indv_probability <- function(alpha, beta, pheromone, attractiveness){
    return((pheromone**alpha)*((1+attractiveness)**beta)) #attractiveness [0,1]^beta will be smaller
  }
  
  #current_neighbors <- names(adj_mat[current_position, adj_mat[current_position, ] > 0])
  current_neighbors <- names(which(adj_mat[current_position, ] > 0))
  #print("Current neighbors: ", row.names = FALSE)
  #print(current_neighbors)
  unvisited_neighbors <- setdiff(current_neighbors, c(route, current_position, visited))
  #print("Unvisited neighbors: ", row.names =  FALSE)
  #print(unvisited_neighbors)
  
  if(remove_close){
    #close_neighbors <- names(distance[unvisited_neighbors][distance[unvisited_neighbors] <= distance[current_position]])
    close_neighbors <- names(which(distance[unvisited_neighbors] <= distance[current_position]))
    unvisited_neighbors <- setdiff(unvisited_neighbors, close_neighbors)
  }
  #print("Possible neighbors: ", row.names = FALSE)
  #print(unvisited_neighbors)
  
  if(!is_empty(unvisited_neighbors)){
    if(first){
      chosen_node <- sample(unvisited_neighbors, size = 1)
      possibility <- get_energy_lost(possible_position = chosen_node, 
                                     route = c(route, current_position), 
                                     oscc_old = oscc_old, 
                                     gene_scores = gene_scores, 
                                     normalized_distance = normalized_distance, 
                                     distance = distance, 
                                     current_capacity = current_capacity,
                                     penalize_dist = penalize_dist)
      #print(paste("Chosen Node: ", chosen_node))
      #print("Movement summary:", row.names = FALSE)
      #print(possibility)
      
      if(possibility$new_capacity >= 0){
        updated_pheromone_map <- pheromone_map
        updated_pheromone_map[current_position, chosen_node] <- updated_pheromone_map[current_position, chosen_node] + possibility$favorability
        if(verbose > 1){
            print(paste("Moved ant from", current_position, "to", chosen_node))
        }
        return(list(current_position = chosen_node, 
                    route = c(route, current_position), 
                    updated_pheromone_map = updated_pheromone_map,
                    current_capacity = possibility$new_capacity
              ))
      } else {
        if(verbose > 1){
            print("No more valid unvisited neighbors for this ant.")
        }
        
        return(move_ant(current_position, route, oscc_old, 
                        gene_scores, normalized_distance, distance, adj_mat,
                        current_capacity, alpha, beta, pheromone_map, penalize_dist, first, remove_close, verbose, visited = c(visited, chosen_node)))
      }
    } else {
      possibilities <- data.frame()
      probabilities <- vector()
      for(i in 1:length(unvisited_neighbors)){
        possibility <- unvisited_neighbors[i]
        possibilities <- rbind(possibilities, 
                               get_energy_lost(possible_position = possibility, 
                                               route = c(route, current_position), 
                                               oscc_old = oscc_old, 
                                               gene_scores = gene_scores, 
                                               normalized_distance = normalized_distance, 
                                               distance = distance, 
                                               current_capacity = current_capacity,
                                               penalize_dist = penalize_dist)
                               )
        if(possibilities[i,]$new_capacity >= 0){
          probabilities[i] <- get_indv_probability(alpha, beta, 
                                               pheromone_map[current_position,possibility], 
                                               possibilities[i,]$favorability)
        } else {
          probabilities[i] <- 0
        }
      }
      
      if(sum(probabilities) > 0 && 
         any(distance[unvisited_neighbors[probabilities > 0]] > distance[current_position])){
        
        normalized_probabilities <- probabilities / sum(probabilities)
        rownames(possibilities) <- unvisited_neighbors
        possibilities <- possibilities %>% mutate(prob = normalized_probabilities)
        chosen_node <- sample(rownames(possibilities), size = 1, prob = possibilities$prob)
        #print(paste("Chosen Node: ", chosen_node))
        updated_pheromone_map <- pheromone_map
        updated_pheromone_map[current_position, chosen_node] <- updated_pheromone_map[current_position, chosen_node] + possibilities[chosen_node,]$favorability
        #print(paste("Moved ant from", current_position, "to", chosen_node))
        #print(paste("New capacity:", possibilities[chosen_node,]$new_capacity))
        return(list(current_position = chosen_node, 
                    route = c(route, current_position), 
                    updated_pheromone_map = updated_pheromone_map,
                    current_capacity = possibilities[chosen_node,]$new_capacity,
                    oscc_current = possibilities[chosen_node,]$oscc_current
                      )
               )
      } else {
        if(verbose > 1){
            print("No more valid unvisited neighbors for this ant.")
        }
        return(FALSE)
      }
    }
  } else {
    if(verbose > 1){
        print("No more unvisited neighbors for this ant.")
    }
    return(FALSE)
  }
}

custom_max <- function(x){
  if(max(x) == 0){
    return(1)
  } else {
    return(max(x))
  }
}

optimize_module <- function(gene_id, 
                            gene_scores, 
                            adj_mat, 
                            distance = NULL, 
                            normalized_distance = NULL,
                            #capacities = NULL,
                            alpha, 
                            beta, 
                            starting_capacity = 1.0,
                            n_iter = 50, 
                            n_ants = NULL, 
                            n_cores = 1,
                            penalize_dist = TRUE,
                            seed = NULL,
                            verbose = 1)
  {
  require(parallel)
  
  if(is.null(distance) || is.null(normalized_distance)){
    require(mND)
    require(igraph)
    int_network <- graph_from_adjacency_matrix(adj_mat, mode = "undirected")
    distance <- calcAllPairsDistances(int_network, directionPaths = "all", networkName = "int_network")
    normalized_distance <- t(as.data.frame(t(distance)) %>% mutate_all(~ ./custom_max(.)))
    #normalized_distance <- normalize_adj_mat(distance)
  }
  
  #if(is.null(capacities)){
  #  capacities <- gene_scores/max(gene_scores)
    #starting_capacity <- gene_scores[gene_id, 1]/max(gene_scores)
  #}
  #starting_capacity <- capacities[gene_id, 1]
  starting_capacity <- starting_capacity
  
  if(is.null(n_ants)){
    # n_ants <- how to initialize number of ants
  }
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  if(sum(adj_mat[gene_id,]) == 0 | (sum(adj_mat[gene_id,]) == 1 & adj_mat[gene_id,gene_id] == 1)){
    print("Isolated node detected and skipped.")
    print("======================================================================")
    return(FALSE)
  }
  
  starting_pheromone_map <- data.frame(matrix(rep(1.0, nrow(adj_mat)*ncol(adj_mat)), 
                                     nrow = nrow(adj_mat), 
                                     ncol = ncol(adj_mat)), 
                              row.names = rownames(adj_mat))
  colnames(starting_pheromone_map) <- colnames(adj_mat)
  
  modules <- rep(list(starting_pheromone_map), n_iter)
  times  <- vector(length = n_iter)
  #if(verbose <= 1){
   #   print("======================================================================")
  #}
  results <- parallel::mclapply(1:n_iter, function(i){ # For each iteration
    start_time <- Sys.time()
    pheromone_map <- starting_pheromone_map
    if(verbose > 1){
        print("======================================================================")
        print(paste("Starting iteraton", i))
    }
    count_ants <- n_ants
    ants <- data.frame(can_move = rep(TRUE, n_ants), 
                       current_position = rep(gene_id, n_ants), 
                       oscc_old = rep(0, n_ants),
                       capacity = rep(starting_capacity, n_ants)
    )
    rownames(ants) <- 1:n_ants
    routes <- rep(list(c()), n_ants)
    
    while(count_ants > 0){ # While there are ants able to move
      shuffled_ants <- sample(1:n_ants, n_ants)
      for(ant in shuffled_ants){ # for all ants
        if(verbose > 1){
            print(paste("Considering ant", ant, "..."))
        }
        #ant_info <- ants[ant,]
        route <- routes[[ant]]
        if(ants[ant,]$can_move){ # if the ant can move
          if(ants[ant,]$capacity == 0){
            ants[ant,]$can_move <- FALSE
            if(is_empty(routes[[ant]])){
              routes[[ant]] <- c(ants[ant,]$current_position)
            }
            count_ants <- count_ants - 1
          } else {
            if(verbose > 1){
                print(paste("Moving ant", ant, "..."))
            }
            if(i == 1){ # if this is the first iteration, move ant randomly
              movement <- move_ant(current_position = ants[ant,]$current_position,
                                   route = route,
                                   oscc_old = ants[ant,]$oscc_old,
                                   gene_scores = gene_scores,
                                   normalized_distance = normalized_distance[gene_id, ],
                                   distance = distance[gene_id, ],
                                   adj_mat = adj_mat,
                                   current_capacity = ants[ant,]$capacity, 
                                   alpha = alpha,
                                   beta = beta,
                                   pheromone_map = pheromone_map,
                                   penalize_dist = penalize_dist,
                                   first = TRUE)
            } else { # if this is any iteration other than the first, move ant normally
              movement <- move_ant(current_position = ants[ant,]$current_position,
                                   route = route,
                                   oscc_old = ants[ant,]$oscc_old,
                                   gene_scores = gene_scores,
                                   normalized_distance = normalized_distance[gene_id, ],
                                   distance = distance[gene_id, ],
                                   adj_mat = adj_mat,
                                   current_capacity = ants[ant,]$capacity, 
                                   alpha = alpha,
                                   beta = beta,
                                   pheromone_map = pheromone_map,
                                   penalize_dist = penalize_dist)
            } 
            if(is.logical(movement)){ # if ant could not move
              if(verbose > 1){
                  print(paste("Ant", ant, "could not move."))
              }
              ants[ant,]$can_move <- FALSE
              #ants[ant,] <- ant_info
              if(is_empty(routes[[ant]])){
                routes[[ant]] <- c(ants[ant,]$current_position)
              }
              count_ants <- count_ants - 1
              #print(ants)
            } else { #if ant could move, update ant info
              ants[ant,]$oscc_old <- movement$oscc_current
              ants[ant,]$current_position <- movement$current_position
              ants[ant,]$capacity <- movement$current_capacity
              #ants[ant,] <- ant_info
              #print(ants)
              routes[[ant]] <- movement$route
              pheromone_map <- movement$updated_pheromone_map
            }
          }
        } else { # if ant cannot move
          if(verbose > 1){
              print(paste("Ant", ant, "cannot move."))
          }
          next()
        }
      } # end of for: done going through each ant once
    } # end of while: all ants cannot move
    
    #modules[[i]] <- unique(unlist(routes)) #changed from pheromone_map
    if(verbose > 1){
        print("No ants left to move.")
        print(paste("Iteration", i, "complete."))
    }
    end_time <- Sys.time()
    #times[i] <- end_time - start_time
    return(list(module = unique(unlist(routes)), time = end_time - start_time, iter = i))
  }, mc.cores = n_cores) # end of for: all iterations done
  
  print(paste("ACO complete."))
  for(i in 1:n_iter){
    modules[[results[[i]]$iter]] <- results[[i]]$module
    times[results[[i]]$iter] <- results[[i]]$time
  }
  #print("======================================================================")
  return(list(modules = modules, times = times))
}

calcAllPairsDistances <- function(network, 
                                  directionPaths="all",
                                  weightVector = NULL, 
                                  networkName)
  {
  # This a function from GS code
  # Args:
  #
  # Returns:
  #
  
  require(igraph)
  
  
  # weightVector = null => use weight attribute if it exists otherwise don't
  shortestPathsMatrix <- shortest.paths(network,
                                        v=V(network),
                                        to=V(network),
                                        mode = directionPaths,
                                        weights = weightVector)
  
  attr(shortestPathsMatrix,"networkName") <- networkName
  return(shortestPathsMatrix)
  
}

SumAbsCor <- function(gene.id,
                      cor.vector,
                      mods
                      ){
  # Calculate sum of cor.vector ----------------------------
  sum.abs.cor <- vapply(mods, function(mod){
    
    if(length(mod) > 1){
      res <- sum( abs( cor.vector[setdiff(mod, gene.id)] ), na.rm= TRUE)
      return(res)
    }else{
      return(0)
    }
    
    
  },
  numeric(1))
  return(sum.abs.cor)
  
}

Observed.SI <- function(gene.id,
                        cor.matrix,
                        mods){
  
  
  
  
  # Vector of cor between j and all other genes on network (excluding j)
  observed.cor <- SumAbsCor(gene.id,
                            cor.matrix[gene.id,],
                            mods)
  
  
}

Resample.SI <- function(gene.id,
                        cor.matrix,
                        mods,
                        num.Sphere.resamples=1000,
                        genes.assayedETnetwork,
                        seed){
  set.seed(seed)
  # Vector of cor between j and all other genes on network (excluding j)
  cor.with.j <- cor.matrix[gene.id,
                           setdiff(genes.assayedETnetwork,gene.id)]
  
  
  x <- rep(0,length(cor.with.j))
  names(x) <- setdiff(genes.assayedETnetwork,gene.id)
  A <- vapply(mods, function(mod){
    
    x[setdiff(mod, gene.id)] <- 1
    return(x)
    
    
  },
  numeric(length(x)) )
  
  A <- t(A)
  
  
  
  # For each reasampling calculating sum of abs cor at each distance
  # Under this new implementaiotn, cor added up in d are also included in d + 1
  # and we only resample 1000 times instead of 34000 times.
  resampled.sum.abs.cor <- replicate(num.Sphere.resamples,
                                     
                                     
                                     # Simplify to [1:34] from [1:34,1]
                                     c( A %*% abs( sample(cor.with.j,
                                                          size = length(cor.with.j),
                                                          replace = TRUE)
                                      ))
                                     
                                    )
  
}

calc_pSphere <- function(gene.id,
                         cor.matrix,
                         mods,
                         genes.assayedETnetwork,
                         num.Sphere.resamples,
                         seed
                         ){
  start_time <- Sys.time()
  observed.cor <- Observed.SI(gene.id,
                              cor.matrix,
                              mods)
  
  
  resampled.cor <- Resample.SI(gene.id,
                               cor.matrix,
                               mods,
                               num.Sphere.resamples,
                               genes.assayedETnetwork,
                               seed = seed)
  end_time <- Sys.time()
  time1 <- end_time - start_time
  time2 <- vector()
  # proportion of null sumabscor \GEQ
  p.Sphere <- vector()
  for(idx in 1:(length(mods))){
    start_time <- Sys.time()
    #observed.p <- observed.cor[[distance]]$cor
    observed.p <- observed.cor[idx]
    
    #null.p <- resampled.cor[[distance]]$cor
    null.p <- resampled.cor[idx,1:num.Sphere.resamples ]
    end_time <- Sys.time()
    time2[idx] <- end_time - start_time + time1
    p.Sphere[idx] <- length(null.p[null.p >= observed.p])/length(null.p)
    
  }
  
  p.Sphere[p.Sphere == 0] <- 1/(num.Sphere.resamples+1)
  return(list(p_spheres = p.Sphere, times = time2))
}

# Calc gene level statistics & a null set of gene level stats (shuffle phenotype labels)
acoCalcGeneTStats <- function(expr,
                           classLabels,
                           numResamples,
                           seed,
                           permutate = TRUE){
  # Calc gene level statistics
  # This code is being written using CurOvGradeKEGGnets[[2]]
  #
  # Args:
  #    expr: is a matrix of genes by samples
  #    classLabels: is a vector of class labels (e.g. high vs low)
  #    numrResamples: number of times the phenotype labels are shuffled
  # Returns:
  #    observedStats: a vector of observed moderated t-statistics
  #    permStats: a matrix of resampled moderated t statistics (resamplings are rows)
  
  # cf. d715_timecourse_contrasts, network_review_GSEAhat ?
  # Should I save the fit so I have the gene p values etc...?
  
  require(limma)
  print("Calculating observed gene stats using limma...")
  desMat <-  model.matrix(~factor(classLabels))
  # treat is a limma function
  # Given a microarray linear model fit, compute moderated t-statistic, etc
  fit <-  treat(lmFit(expr,desMat))
  observedStats = as.data.frame(fit$t[,2])
  
  if(permutate){
    print("Permutating labels and calculating resampled gene stats using limma...")
    if(!is.null(seed)){
      set.seed(seed)
    }
    permStats <- sapply(1:numResamples,function(resampleLoopIndex){
      
      # Shuffle the phenotype labels
      permLabels <- sample(classLabels,replace = FALSE)
      
      #Refit and recalculcate gene level statistics using permLabels
      permDesMat <-  model.matrix(~factor(permLabels))
      permFit <-  treat(lmFit(expr,permDesMat))
      #print(head(permFit$t[,2]))
      return(permFit$t[,2])
      
    })
    
    # Transpose permStats so the rows are resamplings
    permStats <- t(permStats)
    
    # List and return
    
    geneTStats <- list(observed=observedStats, resampled = permStats, complete_res = topTreat(fit, coef = 2, number = nrow(observedStats)))
    
    
  }else{
    geneTStats <- list(observed=observedStats, complete_res = topTreat(fit, coef = 2, number = nrow(observedStats)))
  }
  
  return(geneTStats)
}

optimize_modules <- function(gene_scores_per_sample,
                             adj_mat,
                             n_ants,
                             n_iter,
                             n_cores = 1,
                             starting_capacity = 1.0,
                             alpha,
                             beta,
                             seed,
                             n_resamples = 1000,
                             #rand_sample = FALSE,
                             class_labels,
                             penalize_dist = TRUE){
  
  source('OSCC.R')
  require(limma)
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  print("Preprocessing...")
  gene_scores_per_sample <- gene_scores_per_sample[rownames(gene_scores_per_sample) %in% rownames(adj_mat),]
  #if(rand_sample){
   #   gene_scores_per_sample <- gene_scores_per_sample[rownames(gene_scores_per_sample) %in% sample(rownames(gene_scores_per_sample), size = 500),]
  #}
  cor_mat <- calcCorMatrix(data.frame(gene_scores_per_sample), corMethod = "pearson", exprName = "gene_scores_per_sample", useMethod = "everything")
  adj_mat <- adj_mat[rownames(adj_mat) %in% rownames(gene_scores_per_sample), ]
  adj_mat <- adj_mat[ ,colnames(adj_mat) %in% rownames(gene_scores_per_sample)]
  int_network <- graph_from_adjacency_matrix(adj_mat, mode = "undirected")
  distance <- calcAllPairsDistances(int_network, directionPaths = "all", networkName = "int_network")
  distance[!is.finite(distance)] <- -1
  normalized_distance <- t(data.frame(distance) %>% mutate_all(~ ./custom_max(.)))
  rownames(normalized_distance) <- rownames(distance)
  #print(paste("ARE ALL 0?: ", sum(normalized_distance == 0, na.rm = TRUE)))
  
  calculate_pdecay <- function(observed_oscc,
                               resampled_scores,
                               distance,
                               chosen_nodes){
    null_oscc <- vector()
    for(i in 1:dim(resampled_scores)[1]){
      null_oscc[i] <- calculate_OSCC(as.numeric(distance[chosen_nodes]), abs(as.numeric(resampled_scores[i, chosen_nodes])))
    }
    return(sum(null_oscc <= observed_oscc)/length(null_oscc))
  }
  
  gene_stats <- acoCalcGeneTStats(gene_scores_per_sample, class_labels, n_resamples, seed)
  gene_scores <- gene_stats$observed
  permStats <- gene_stats$resampled
  limma_adj_pvals <- gene_stats$complete_res[,5]
  names(limma_adj_pvals) <- rownames(gene_stats$complete_res)
  
  results <- list()
  count <- 1
  for(gene_id in rownames(gene_scores)){
    time1 <- Sys.time()
    print(paste("Starting module optimization for gene #", count, "/", nrow(gene_scores), ": ", gene_id,"..."))
    count <- count + 1
    modules <- optimize_module(gene_id = gene_id,
                               gene_scores = gene_scores,
                               adj_mat = adj_mat,
                               distance = distance,
                               normalized_distance = normalized_distance,
                               alpha = alpha,
                               beta = beta,
                               starting_capacity = starting_capacity,
                               n_iter = n_iter,
                               n_ants = n_ants,
                               n_cores = n_cores,
                               penalize_dist = penalize_dist,
                               seed = seed)
    if(!is.logical(modules)){ #If gene is not isolated
      osccs <- vector()
      nums <- vector()
      denoms <- vector()
      common_terms <- vector()
      p_decays <- vector()
      pDecay_times <- vector()
      
      print("Calculating p-Decay...")
      pdecay_res <- mclapply(1:length(modules$modules), function(i){
        start_time <- Sys.time()
        chosen_nodes <- modules$modules[[i]]
        if(length(chosen_nodes) == 1){
          print(paste("WARNING: only 1 chosen node:", chosen_nodes))
          oscc <- 1
          num <- 1
          denom <- 1
          common_term <- 1
          p_decay <- 1
        } else {
            oscc_res <- calculate_OSCC(as.numeric(as.data.frame(distance)[gene_id, chosen_nodes]), abs(as.numeric(gene_scores[chosen_nodes, ])), full_res = TRUE)
            oscc <- oscc_res$oscc
            num <- oscc_res$numerator
            denom <- oscc_res$denominator
            common_term <- oscc_res$common_term
            p_decay <- calculate_pdecay(oscc, permStats, as.data.frame(distance)[gene_id, ], chosen_nodes)
        }
        
        pDecay_time <- Sys.time() - start_time + modules$times[i]
        return(list(oscc = oscc, num = num, denom = denom, common_term = common_term, p_decay = p_decay, pDecay_time = pDecay_time, iter = i))
      }, mc.cores = n_cores)
      
      for(i in 1:length(modules$modules)){
        osccs[pdecay_res[[i]]$iter] <- pdecay_res[[i]]$oscc
        nums[pdecay_res[[i]]$iter] <- pdecay_res[[i]]$num
        denoms[pdecay_res[[i]]$iter] <- pdecay_res[[i]]$denom
        common_terms[pdecay_res[[i]]$iter] <- pdecay_res[[i]]$common_term
        p_decays[pdecay_res[[i]]$iter] <- pdecay_res[[i]]$p_decay
        pDecay_times[pdecay_res[[i]]$iter] <- pdecay_res[[i]]$pDecay_time
      }
      
      print("Calculating p-Sphere...")
      psphere_res <- calc_pSphere(gene_id, cor_mat, modules$modules, rownames(adj_mat), n_resamples, seed)
      p_spheres <- psphere_res$p_spheres
      pSphere_times <- psphere_res$times
      
      print(paste("Done for gene", gene_id,":"))
      p_decays[p_decays == 0] <- 1/(n_resamples+1)
      stats <- cbind(osccs, nums, denoms, common_terms, p_decays, pDecay_times, p_spheres, pSphere_times, p_fishers = pchisq(-2*(log(p_spheres) + log(p_decays)), df = 4, lower.tail = FALSE))
      print(stats)
      results[[gene_id]] <- list(stats = stats, modules = modules, total_time = Sys.time() - time1)
      print("======================================================================")
    } else {
      #results[[gene_id]] <- "Isolated"
    }
  }
  print("Done for all genes")
  return(results)
}


#list(current_position = chosen_node, 
#route = c(route, current_position), 
#updated_pheromone_map = updated_pheromone_map,
#current_capacity = possibilities[chosen_node,]$new_capacity,
#oscc_current = possibilities[chosen_node,]$oscc_current)
