# ===========================================================================
# file: GeneSurrounder
# description:
# requires: 
# author: Sahil Shah <sahil.shah@u.northwestern.edu>
# ==========================================================================

#===========================================#
#===========================================#
require(pcaPP)

#' Pre-processing Step
#'
#' Before applying GeneSurrounder, the observed and resampled
#' differential expression of the genes should be calculated
#'
#' @param expr A matrix (genes by samples) of expression values.
#' @param classLabels A factor with levels corresponding to class labels.
#' @param numResamples defaults to 1000. The number of resamples when calculating resampled differential expression.
#'



# Calc gene level statistics & a null set of gene level stats (shuffle phenotype labels)
calcGeneTStats <- function(expr,classLabels,numResamples=1000, seed){
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
    set.seed(seed)
    require(limma)
    
    desMat =  model.matrix(~factor(classLabels))
    # treat is a limma function
    # Given a microarray linear model fit, compute moderated t-statistic, etc
    fit =  treat(lmFit(expr,desMat))
    observedStats = fit$t[,2]
    

    permStats <- sapply(1:numResamples,function(resampleLoopIndex){
        
        # Shuffle the phenotype labels
        permLabels <- sample(classLabels,replace = FALSE)
        
        #Refit and recalculcate gene level statistics using permLabels
        permDesMat =  model.matrix(~factor(permLabels))
        permFit =  treat(lmFit(expr,permDesMat))
        #print(head(permFit$t[,2]))
        return(permFit$t[,2])
        
    })

    # Transpose permStats so the rows are resamplings
    permStats = t(permStats)
    
    # List and return
    geneTStats = list(observed=observedStats,resampled = permStats)
    
    return(geneTStats)
    
}


#' Pre-processing Step
#'
#' Before applying GeneSurrounder, the distances on the
#' global network should be calculated
#'
#' @param network An igraph network on which to calculate the distances.
#' @param directionPaths defaults to "all". A string indicating how the distances should be calculated.
#' @param weightVector defaults to NULL. A numeric vector giving edge weights.
#' @param networkName A string containg the name of the network.
#'



calcAllPairsDistances <- function(network,directionPaths="all",weightVector = NULL,networkName){
    #
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


#' Pre-processing Step
#'
#' Before applying GeneSurrounder, the correlation
#' between the expression of the genes should be calculated
#'
#' @param exprMatrix A matrix (genes by samples) of expression values.
#' @param corMethod A string contating the correlation method to use.
#' @param exprName A string containg the name of the expression matrix.
#'


calcCorMatrix <- function(exprMatrix,corMethod,exprName,useMethod){
    #
    # Args:
    #
    # Returns:
    #
    
    #=======================================
    #=======================================
    
    # cor calculates correlations between *columns* => take transpose of exprMatrix
    corMatrix <- cor(t(exprMatrix),method=corMethod,use=useMethod)
    
    # Meta data for the for the corMatrix
    attr(corMatrix,"expr") <- exprName
    attr(corMatrix,"method") <- corMethod
    
    return(corMatrix)
    
}


#' Sphere of Influence Step
#'
#' Sphere of Influence assesses if a
#' candidate gene i meets the first criterion by testing if gene i
#' is more strongly correlated with its network neighbors than with a
#' random set of genes
#'
#' @param gene.id The name of the gene to which GeneSurrounder is applied
#' @param distance.matrix A matrix of the distances on the global network
#' @param cor.matrix A matrix of correlations between the expression of the genes
#' @param diameter The diameter of the global network
#' @param genes.assayedETnetwork The names of the genes that are assayed and on the network
#'



Observed.SI <- function(gene.id,
                        distance.matrix,
                        cor.matrix,
                        diameter,
                        genes.assayedETnetwork){




    # Vector of cor between j and all other genes on network (excluding j)
    cor.with.j <- cor.matrix[gene.id,
                             setdiff(genes.assayedETnetwork,gene.id)]

    distances.to.j <- distance.matrix[gene.id,
                                      setdiff(genes.assayedETnetwork,gene.id)]


    observed.cor <- SumAbsCor(cor.with.j,
                                diameter,
                                distances.to.j)


}



#' @param gene.id The name of the gene to which GeneSurrounder is applied
#' @param distance.matrix A matrix of the distances on the global network
#' @param cor.matrix A matrix of correlations between the expression of the genes
#' @param diameter The diameter of the global network
#' @param num.Sphere.resamples defaults to 1000. The number of resamples when running the Sphere of Influence Procedure
#' @param genes.assayedETnetwork The names of the genes that are assayed and on the network
#'


 Resample.SI <- function(gene.id,
                        distance.matrix,
                        cor.matrix,
                        diameter,
                        num.Sphere.resamples=1000,
                        genes.assayedETnetwork,
                        seed){
    set.seed(seed)
    # Vector of cor between j and all other genes on network (excluding j)
    cor.with.j <- cor.matrix[gene.id,
                             setdiff(genes.assayedETnetwork,gene.id)]

    distances.to.j <- distance.matrix[gene.id,
                                      setdiff(genes.assayedETnetwork,gene.id)]


    x <- rep(0,length(cor.with.j))

    A <- vapply(1:diameter, function(distance){

        x[distances.to.j <= distance] <- 1

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
                                                )
                                          )

                                        # SumAbsCor(sample(cor.with.j,
                                        #                     size = length(cor.with.j),
                                  #                        replace = TRUE),
                                        #           diameter,
                                        #           distances.to.j
                                        #           )

                                     )


    # apply(1:num.resamples, 1, function(resample.num){

    #     # Sample with replacement form cor.with.j ----------------------------
    #     resampled.cor <- sample(cor.with.j,
    #                             size = length(cor.with.j),
    #                             replace = TRUE)

    #     # Calculate sum of resampled cor ----------------------------
    #     sum.abs.cor <- apply(1:diameter, 1, function(distance){

    #         #
    #         sum( abs( resampled.cor[distances.to.j <= distance] ) )


    #         })


    #     return(sum.abs.cor)

        
    # })

}

#' Helper function for Sphere of Influence procedure
#'
#' \code{SumAbsCor} returns the total observed correlation between gene i and its neighbors
#'
#' @param cor.vector A vector of correlations
#' @param diameter The diameter of the global network
#' @param distances.to.j A vector of distances to the gene to which GeneSurrounder is applied
#'
#'

SumAbsCor <- function(cor.vector,diameter,distances.to.j){


    # Calculate sum of cor.vector ----------------------------
    sum.abs.cor <- vapply(1:diameter, function(distance){

        #
        sum( abs( cor.vector[distances.to.j <= distance] ), na.rm= TRUE)


    },
    numeric(1))

    return(sum.abs.cor)


}


#' Decay of Differential Expression step
#'
#' Decay of Differential Expression, tests whether the magnitude of 
#' differential expression of other genes j in the neighborhood is inversely 
#' related to the distance d(i,j) of gene j from gene i
#'
#' @param distance.matrix A matrix of the distances on the global network
#' @param gene.id The name of the gene to which GeneSurrounder is applied 
#' @param genes.assayedETnetwork The names of the genes that are assayed and on the network
#' @param diameter The diameter of the global network 
#' @param geneStats.observed A vector of the observed differential expression
#'

Observed.DecayDE <- function(distance.matrix,
							 gene.id,
							 genes.assayedETnetwork,
							 diameter,
							 geneStats.observed){



	distances <- distance.matrix[gene.id,
							 genes.assayedETnetwork]


	observed.tau_b <- vapply(1:diameter,function(RADIUS){
	  #Ghadi (1 line)
	  #idx <- which(distances <= RADIUS & distances > 0, arr.ind = TRUE)
    
		# Excludes gene j with distances > 0 
		igenes.distances <- distances[distances <= RADIUS
									 & distances > 0]
	  igenes.names <- names(igenes.distances)
	  
	  #Ghadi (2 lines)
		#igenes.names <- unique(rownames(idx))
		#igenes.distances <- distances[igenes.names,igenes.names]


		return(

			cor.fk(abs(geneStats.observed[igenes.names]), igenes.distances)

			)

	},
	numeric(1))

}


#' @param distance.matrix A matrix of the distances on the global network
#' @param gene.id The name of the gene to which GeneSurrounder is applied
#' @param genes.assayedETnetwork The names of the genes that are assayed and on the network
#' @param diameter The diameter of the global network
#' @param perm.geneStats.matrix A matrix of the resampled differential expression
#' @param sizes The number of genes assayed and on the network in each neighborhood
#'


Resample.DecayDE <- function(distance.matrix,
                             gene.id,
                             genes.assayedETnetwork,
                             diameter,
                             perm.geneStats.matrix,
                             sizes){


    distances <- distance.matrix[gene.id,
                             genes.assayedETnetwork]





    num.genes <- length(genes.assayedETnetwork) - 1
    geneid.d <- which(sizes == num.genes)[1] #This outputs NA when diameter < 8, sizes gives number close to gene numbers but not quite, which() cannot find it

    # null.tau b matrix of 1000 row 34 columns b/c stacks columns
    #null.tau_b <- vapply(1:geneid.d,function(RADIUS){
    null.tau_b <- vapply(1:diameter,function(RADIUS){ #Ghadi

        # Excludes gene j with distances > 0
        igenes.distances <- distances[distances <= RADIUS
                                      & distances > 0]
        igenes.names <- names(igenes.distances)





        null <- vapply(1:nrow(perm.geneStats.matrix),function(resample.index){

            return ( cor.fk( abs(perm.geneStats.matrix[resample.index,igenes.names]) ,
                            igenes.distances)

                   )

        },
        numeric(1))

        # X <- cbind(igenes.distances,
        #       t( abs(geneStats$resampled[,igenes.names]) )
        #       )

        # Y <- cor.fk(X)


    },
    numeric(nrow(perm.geneStats.matrix))) # numeric is length of one of the columns vapply stacks

    # 2:genid.d means number of rows in null.tau_b is not geneid.d
    
    null.tau_b <- cbind(null.tau_b,
                    matrix(
                        #rep(null.tau_b[, geneid.d],diameter - geneid.d),
                        rep(null.tau_b[, geneid.d],diameter), #Ghadi
                        nrow = nrow(perm.geneStats.matrix))

                        )



    # cf. ?vapply
    # If FUN.VALUE is not an array, the result is a matrix with
    # length(FUN.VALUE) rows and length(X) columns,

    # return matrix distance by resample number a la Resample.SI
    return(t(null.tau_b))


}


# 2: gene id. : cbind first row.
# preallocate matrix size instead of cbind?


#' GeneSurrounder wrapper
#'
#' The GeneSurrounder method consists of two tests that are run
#' independently of each other and then combined to determine if
#' the putative disease gene is a “disruptive” candidate disease gene
#' meeting both criteria.
#'
#' @param distance.matrix A matrix of the distances on the global network
#' @param cor.matrix A matrix of correlations between the expression of the genes
#' @param geneStats.observed A vector of the observed differential expression
#' @param perm.geneStats.matrix A matrix of the resampled differential expression
#' @param genes.assayedETnetwork The names of the genes that are assayed and on the network
#' @param diameter The diameter of the global network
#' @param num.Sphere.resamples The number of resamples when running the Sphere of Influence Procedure
#' @param gene.id The name of the gene to which GeneSurrounder is applied


geneNIDG <- function(distance.matrix,
                     cor.matrix,
                     geneStats.observed,
                     perm.geneStats.matrix,
                     genes.assayedETnetwork,
                     diameter,
                     num.Sphere.resamples,
                     gene.id,
                     seed){


    # size -------------------------------------------------------------------

    start_time <- Sys.time()
    distances <- distance.matrix[gene.id,
                             genes.assayedETnetwork]


    sizes <- vapply(1:diameter,function(RADIUS){

        igenes.distances <- distances[distances <= RADIUS
                                      & distances > 0]

        length(igenes.distances)


    },
    numeric(1))



    # gene.DecayDE --------------------------------------------------------------

    # vector of observed tau b and null tau b a la observed.cor and resampled.cor

    observed.tau_b <- Observed.DecayDE(distance.matrix,
                                         gene.id,
                                         genes.assayedETnetwork,
                                         diameter,
                                         geneStats.observed)


    null.tau_b <- Resample.DecayDE(distance.matrix,
                                 gene.id,
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

    end_time <- Sys.time()
    tot_time <- end_time - start_time
    # gene.SI ----------------------------------------------------------------


    observed.cor <- Observed.SI(gene.id,
                                    distance.matrix,
                                    cor.matrix,
                                    diameter,
                                    genes.assayedETnetwork)


    resampled.cor <- Resample.SI( gene.id,
                                      distance.matrix,
                                      cor.matrix,
                                      diameter,
                                      num.Sphere.resamples,
                                      genes.assayedETnetwork,
                                  seed = seed)

    # proportion of null sumabscor \GEQ
    p.Sphere <- vapply(1:diameter,function(distance){

         #observed.p <- observed.cor[[distance]]$cor
         observed.p <- observed.cor[distance]
        
        #null.p <- resampled.cor[[distance]]$cor
        null.p <- resampled.cor[distance,1:num.Sphere.resamples ]

        return(length(null.p[null.p >= observed.p])/length(null.p))


     },
     numeric(1))

    p.Sphere[p.Sphere == 0] <- 1/(num.Sphere.resamples+1)


    # no base R combine p values ? combine after running --------------------

    # ------------------------------------------------------------------------

    # ?adply
    # The most unambiguous behaviour is achieved when .fun returns a data frame
    
    return(list(res = data.frame(gene.id = rep(gene.id,diameter),
                      radius = 1:diameter,
                      size = sizes,
                      observed.tau_b = observed.tau_b,
                      p.Decay = p.Decay,
                      observed.cor = observed.cor,
                      p.Sphere = p.Sphere,
                      p.Fisher = pchisq(-2*(log(p.Sphere) + log(p.Decay)), df = 4, lower.tail = FALSE)
                      ),
                pdecay_time = tot_time
                )
           )


}


#' Plot results of GeneSurrounder
#'
#' \code{plotRadiusVS} plots the results of GeneSurrounder using the output of the geneNIDG function
#
#' @param geneNIDG.gene A dataframe that is the output of the geneNIDG function
#'


plotRadiusVS <- function(geneNIDG.gene){

geneNIDG.hsa4171 <- geneNIDG.gene

par(mfcol=c(4,1),mar=c(1,4,2,1), oma=c(2,2,2,2),cex.axis = 1.5)


plot(geneNIDG.hsa4171$radius,
-log10(geneNIDG.hsa4171$p.Sphere),
type = 'b',
cex = 2,
pch = 21,
lwd = 2,
bg = 'lightblue',
xlab = '',
ylab = '-log10(p.Sphere)',
ylim = c(0,5),
main = 'MCM2 (DNA replication factor) identified as "mechanistic driver" ',
xaxt='n')

abline(h = -log10(0.01),lty=3,lwd=2)
abline(h = -log10(0.05),lty=4,lwd=2)

# http://seananderson.ca/2013/10/21/panel-letters.html

mtext("(A)",
      side = 3,
      adj = 0.05,
      line = -2.25,
      cex = 1)

# axis(4)
# mtext("(A)", side=4, line=3)


# radius vs p.Decay -------------------------------------------------------

plot(geneNIDG.hsa4171$radius,
-log10(geneNIDG.hsa4171$p.Decay),
type = 'b',
cex = 2,
pch = 21,
lwd = 2,
bg = 'lightblue',
xlab = '',
ylab = '-log10(p.Decay)',
ylim = c(0,5),
xaxt='n')


abline(h = -log10(0.01),lty=3,lwd=2)
abline(h = -log10(0.05),lty=4,lwd=2)

# http://seananderson.ca/2013/10/21/panel-letters.html

mtext("(B)",
      side = 3,
      adj = 0.05,
      line = -2.25,
      cex = 1)



# radius vs p.Fisher -------------------------------------------------------

plot(geneNIDG.hsa4171$radius,
-log10(geneNIDG.hsa4171$p.Fisher),
type = 'b',
cex = 2,
pch = 21,
lwd = 2,
bg = 'lightblue',
xlab = '',
ylab = '-log10(p.NIDG)',
# ylim = c(0,5),
xaxt='n')


abline(h = -log10(0.01),lty=3,lwd=2)
abline(h = -log10(0.05),lty=4,lwd=2)

# http://seananderson.ca/2013/10/21/panel-letters.html

mtext("(C)",
      side = 3,
      adj = 0.05,
      line = -2.25,
      cex = 1)


# radius vs size      -------------------------------------------------------


plot(geneNIDG.hsa4171$radius,
 geneNIDG.hsa4171$size,
type = 'b',
cex = 2,
pch = 21,
bg = 'lightblue',
xlab = 'Neighborhood Radius',
ylab = 'Number of Assayed Genes',
lwd = 2)


abline(h=2708)

# http://seananderson.ca/2013/10/21/panel-letters.html


mtext("(D)",
      side = 3,
      adj = 0.05,
      line = -2.25,
      cex = 1)


}

#################
### from: https://github.com/sahildshah1/gene-surrounder
###
### Shah, S. D., & Braun, R. (2019). GeneSurrounder: network-based identification of disease genes in expression data. BMC bioinformatics, 20(1), 1-12.
