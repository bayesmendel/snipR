## slice algorithm to calculate GMD
gmd <- function(R, S, fs, fm){
  fams <- sort(unlist(R))
  R.sizes <- as.numeric(sapply(R, length))
  
  ## M = map from each family to cluster number in R
  M <- data.frame(FamID = fams, Cluster = rep(0, length(fams)))
  g <- rep(seq_along(R), R.sizes)
  M$Cluster <- g[match(fams, unlist(R))]
  
  ## Compute cost
  cost <- 0
  for(i in 1:length(S)){
    # pMap = list of clusters in R that contain the families in 
    # the i-th cluster of S, along with their counts
    pMap <- setNames(data.frame(matrix(nrow = 0, ncol = 2)), c("Cluster", "Count"))
    for(r in S[[i]]){
      # if we haven't seen this cluster before, add it to pMap
      if(!(M$Cluster[M$FamID == r] %in% pMap$Cluster)){
        pMap[nrow(pMap) + 1, ] <- list(NA, 0)
        pMap$Cluster[nrow(pMap)] <- M$Cluster[M$FamID == r]
      }
      # increment count
      pMap$Count[pMap$Cluster == M$Cluster[M$FamID == r]] <-
        pMap$Count[pMap$Cluster == M$Cluster[M$FamID == r]] + 1
    }
    
    # compute cost to generate i-th cluster of S
    SiCost <- 0
    totalRecs <- 0
    for(j in 1:nrow(pMap)){
      # add cost to split the j-th cluster of R in pMap
      size.Rj <- length(which(M$Cluster == pMap$Cluster[j]))
      if(size.Rj > pMap$Count[j]){
        SiCost <- SiCost + fs(pMap$Count[j], size.Rj - pMap$Count[j])
      }
      size.Rj <- size.Rj - pMap$Count[j]
      if(totalRecs != 0){
        # cost to merge into i-th cluster of S
        SiCost <- SiCost + fm(pMap$Count[j], totalRecs)
      }
      totalRecs <- totalRecs + pMap$Count[j]
    }
    cost <- cost + SiCost
  }
  return(cost)
}

## Pairwise precision
pairPrecision <- function(R, S){
  pairsR <- clevr::clusters_to_pairs(R)
  pairsS <- clevr::clusters_to_pairs(S)
  if(nrow(pairsR) > 0){
    if(nrow(pairsS) > 0){
      return(clevr::precision_pairs(pairsS, pairsR))
    } else{
      return(0)
    }
  } else{
    warning("There are no pairs in the proposed clustering.")
    return(NaN)
  }
}

## Pairwise recall
pairRecall <- function(R, S){
  pairsR <- clevr::clusters_to_pairs(R)
  pairsS <- clevr::clusters_to_pairs(S)
  if(nrow(pairsS) > 0){
    if(nrow(pairsR) > 0){
      return(clevr::recall_pairs(pairsS, pairsR))
    } else{
      return(0)
    }
  } else{
    warning("There are no pairs in the true clustering.")
    return(NaN)
  }
}

## Pairwise F1
pairF1 <- function(R, S){
  pairsR <- clevr::clusters_to_pairs(R)
  pairsS <- clevr::clusters_to_pairs(S)
  if(nrow(pairsR) > 0 & nrow(pairsS) > 0){
    return(clevr::f_measure_pairs(pairsS, pairsR))
  } else{
    warning("There are no pairs in either the proposed or true clustering.")
    return(NaN)
  }
}



# list2Cluster <- function(R){
#   fams <- unlist(R)
#   R.cluster <- length(fams)
#   for(i in 1:length(fams)){
#     R.cluster[i] <- which(sapply(R, `%in%`, x = fams[i]))
#   }
#   return(R.cluster)
# }

## Cluster precision
clusterPrecision <- function(R, S){
  num.both <- length(Reduce(intersect, list(R, S)))
  return(num.both / length(R))
}

## Cluster recall
clusterRecall <- function(R, S){
  num.both <- length(Reduce(intersect, list(R, S)))
  return(num.both / length(S))
}

## Cluster F1
clusterF1 <- function(R, S, prec = NULL, rec = NULL){
  if(!is.null(prec) & !is.null(rec)){
    return(2 * prec * rec / (prec + rec))
  } else{
    prec <- clusterPrecision(R, S)
    rec <- clusterRecall(R, S)
    return(2 * prec * rec / (prec + rec))
  }
}

#' Evaluates performance metrics on the deduplicated data.
#'
#' @param object Object of class \code{\link{Duplicates}}
#' @param requestID Character string denoting the variable that stores the request ID
#' @param famID Character string denoting the variable that stores the family ID
#' @param thresh.i Vector of the thresholds for the intersection method.
#' Use \code{NULL} if the intersection method was not used.
#' @param thresh.g Vector of the thresholds for the greedy scoring method.
#' Use \code{NULL} if the greedy scoring method was not used.
#' @param fs Cost function for splitting for the GMD metric
#' @param fm Cost function for splitting for the GMD metric
#' @return A data frme consisting of the performance metrics for each threshold.
#' @details Outputs 7 metrics: pairwise precision, pairwise recall, pairwise
#' F1, cluster precision, cluster recall, cluster F1, and generalized merge distance
#' (GMD). See Menestrina et al. (2010) for details.
#' 
#' @references Menestrina, D., Whang, S. E., & Garcia-Molina, H. (2010).
#' Evaluating entity resolution results. Proceedings of the VLDB Endowment, 3(1-2), 208-219.
#'
#' @export
summaryMetrics <- function(object, requestID, famID,
                           thresh.i = NULL, thresh.g = NULL,
                           fs = function(x, y) 1,
                           fm = function(x, y) 1){
  
  # true clustering of the families
  S <- as.list(unique(object$rawData[, eval(famID)]))
  dups <- setdiff(object$rawData[, eval(requestID)], object$rawData[, eval(famID)])
  if(length(dups) > 0){
    for(i in 1:length(dups)){
      cluster <- which(sapply(S, `%in%`, x = object$rawData[object$rawData[, eval(requestID)] == dups[i], eval(famID)][1]))
      S[[cluster]] <- c(S[[cluster]], dups[i])
    }
  }
  
  cnames <- c("threshold", "threshType", "pairPrecision", "pairRecall", "pairF1",
              "clusterPrecision", "clusterRecall", "clusterF1", "GMD")
  metrics <- setNames(data.frame(matrix(NA, length(thresh.i) + length(thresh.g),
                                        length(cnames))), cnames)
  metrics$threshold <- c(thresh.i, thresh.g)
  metrics$threshType <- c(rep("intersection", length(thresh.i)),
                          rep("greedy", length(thresh.g)))
  
  if(!is.null(thresh.i)){
    for(i in 1:length(thresh.i)){
      R <- object$dupsList$intersection[[paste0("thresh_", thresh.i[i])]]
      if(length(R) / length(S) > 0.5){
        metrics$pairPrecision[metrics$threshType == "intersection" &
                                metrics$threshold == thresh.i[i]] <- pairPrecision(R, S)
        metrics$pairRecall[metrics$threshType == "intersection" &
                             metrics$threshold == thresh.i[i]] <- pairRecall(R, S)
        metrics$pairF1[metrics$threshType == "intersection" & 
                         metrics$threshold == thresh.i[i]] <- pairF1(R, S)
        metrics$clusterPrecision[metrics$threshType == "intersection" &
                                   metrics$threshold == thresh.i[i]] <- clusterPrecision(R, S)
        metrics$clusterRecall[metrics$threshType == "intersection" &
                                metrics$threshold == thresh.i[i]] <- clusterRecall(R, S)
        metrics$clusterF1[metrics$threshType == "intersection" &
                            metrics$threshold == thresh.i[i]] <- clusterF1(R, S,
                                                                           prec = metrics$clusterPrecision[metrics$threshType == "intersection" &
                                                                                                             metrics$threshold == thresh.i[i]],
                                                                           rec = metrics$clusterRecall[metrics$threshType == "intersection" &
                                                                                                         metrics$threshold == thresh.i[i]])
        metrics$GMD[metrics$threshType == "intersection" &
                      metrics$threshold == thresh.i[i]] <- gmd(R, S, fs, fm)
      }
    }
  }
  if(!is.null(thresh.g)){
    for(i in 1:length(thresh.g)){
      R <- object$dupsList$greedy[[paste0("thresh_", thresh.g[i])]]
      if(length(R) / length(S) > 1){
        metrics$pairPrecision[metrics$threshType == "greedy" &
                                metrics$threshold == thresh.g[i]] <- pairPrecision(R, S)
        metrics$pairRecall[metrics$threshType == "greedy" &
                             metrics$threshold == thresh.g[i]] <- pairRecall(R, S)
        metrics$pairF1[metrics$threshType == "greedy" &
                         metrics$threshold == thresh.g[i]] <- pairF1(R, S)
        metrics$clusterPrecision[metrics$threshType == "greedy" &
                                   metrics$threshold == thresh.g[i]] <- clusterPrecision(R, S)
        metrics$clusterRecall[metrics$threshType == "greedy" &
                                metrics$threshold == thresh.g[i]] <- clusterRecall(R, S)
        metrics$clusterF1[metrics$threshType == "greedy" &
                            metrics$threshold == thresh.g[i]] <- clusterF1(R, S,
                                                                           prec = metrics$clusterPrecision[metrics$threshType == "greedy" &
                                                                                                             metrics$threshold == thresh.g[i]],
                                                                           rec = metrics$clusterRecall[metrics$threshType == "greedy" &
                                                                                                         metrics$threshold == thresh.g[i]])
        metrics$GMD[metrics$threshType == "greedy" &
                      metrics$threshold == thresh.g[i]] <- gmd(R, S, fs, fm)
      }
    }
  }
  return(metrics)
}