#' Reveals the performance of the deduplication process.
#'
#' @param object Object of class \code{\link{Duplicates}} or \code{\link{Scores}}
#' @param requestID Character string denoting the variable that stores the request ID
#' @param famID Character string denoting the variable that stores the family ID
#' @return An object of class \code{\link{Truth}}.
#'
#' @export
revealTruth <- function(object, thresh, requestID, famID) {
  UseMethod("revealTruth", object)
}

#' @rdname revealTruth
#' @export
revealTruth.default <- function(object, thresh, requestID, famID) {
  print("Truth revealable on Duplicates objects only.")
  return(NULL)
}

#' @rdname revealTruth
#' @export
revealTruth.Duplicates <- function(object, requestID, famID) {
  tdat1 <- object$rawData %>%
    as_tibble() %>%
    select(!!sym(famID), !!sym(requestID)) %>%
    rename(famID1 = !!sym(famID), ID1 = !!sym(requestID))
  tdat2 <- object$rawData %>%
    as_tibble() %>%
    select(!!sym(famID), !!sym(requestID)) %>%
    rename(famID2 = !!sym(famID), ID2 = !!sym(requestID))
  
  dups <- object$rawData %>%
    mutate(Duplicate = as.numeric(!!sym(famID) != !!sym(requestID))) %>%
    rename(famID = !!sym(famID), requestID = !!sym(requestID)) %>%
    group_by(famID) %>%
    slice_max(Duplicate, with_ties = FALSE) %>%
    ungroup() %>%
    select(famID, requestID, Duplicate)
  
  ids <- object$rawData %>%
    select(!!sym(famID), !!sym(requestID)) %>%
    rename(famID = !!sym(famID), requestID = !!sym(requestID))
  
  ## comparing the requestID and famID
  if(object$method == "both"){ # if there are both intersection and greedy scores
    nbs.truth <- inner_join(object$Neighbors$intersection %>%
                              as_tibble() %>%
                              mutate(ID1 = as.numeric(ID1), ID2 = as.numeric(ID2)) %>%
                              left_join(tdat1, by = 'ID1') %>%
                              left_join(tdat2, by = 'ID2') %>%
                              rename(matchScore.i = matchScore),
                            object$Neighbors$greedy %>%
                              as_tibble() %>%
                              mutate(ID1 = as.numeric(ID1), ID2 = as.numeric(ID2)) %>%
                              left_join(tdat1, by = 'ID1') %>%
                              left_join(tdat2, by = 'ID2') %>%
                              rename(matchScore.g = matchScore),
                            by = c("ID1", "ID2", "famID1", "famID2"))
    fams.truth <- data.frame(famID = rep(unique(dups$famID),
                                         each = length(c(object$details$intersection,
                                                         object$details$greedy))),
                             threshold = rep(c(object$details$intersection, object$details$greedy),
                                             length(unique(dups$famID))),
                             threshType = rep(c(rep("intersection", length(object$details$intersection)),
                                                rep("greedy", length(object$details$greedy))),
                                              length(unique(dups$famID))))
    fams.truth <- left_join(fams.truth, dups %>% select(famID, Duplicate), by = "famID")
    fams.truth <- fams.truth %>% mutate(Represented = 0, Extra = 0)
  } else{
    nbs.truth <- object$Neighbors %>%
      as_tibble() %>%
      mutate(ID1 = as.numeric(ID1), ID2 = as.numeric(ID2)) %>%
      left_join(tdat1, by = 'ID1') %>%
      left_join(tdat2, by = 'ID2')
    
    if(object$method == "intersection"){
      fams.truth <- data.frame(famID = rep(unique(dups$famID),
                                           each = length(object$details$intersection)),
                               threshold = rep(object$details$intersection,
                                               length(unique(dups$famID))),
                               threshType = rep(rep("intersection", length(object$details$intersection)),
                                                length(unique(dups$famID))))
      fams.truth <- left_join(fams.truth, dups %>% select(famID, Duplicate), by = "famID")
      fams.truth <- fams.truth %>% mutate(Represented = 0)
    } else if(object$method == "greedy"){
      fams.truth <- data.frame(famID = rep(unique(dups$famID),
                                           each = length(object$details$greedy)),
                               threshold = rep(object$details$greedy,
                                               length(unique(dups$famID))),
                               threshType = rep(rep("greedy", length(object$details$greedy)),
                                                length(unique(dups$famID))))
      fams.truth <- left_join(fams.truth, dups %>% select(famID, Duplicate), by = "famID")
      fams.truth <- fams.truth %>% mutate(Represented = 0)
    }
  }
  
  ## categorizing the neighbors
  if(object$method == "both"){
    for(i in 1:length(object$details$intersection)){
      nbs.truth.i <- nbs.truth %>%
        mutate(threshold = object$details$intersection[i], threshType = "intersection",
               dupType = case_when(famID1 == famID2 & matchScore.i >= threshold ~ 'Called Duplicate',
                                   famID1 == famID2 & matchScore.i < threshold ~ 'Missed Duplicate',
                                   is.na(ID1) | is.na(ID2) ~ 'Missed Duplicate',
                                   matchScore.i >= threshold ~ 'False Positive',
                                   matchScore.i < threshold ~ 'True Negative'))
      if(i == 1){
        nbs.truth.tot <- nbs.truth.i
      } else{
        nbs.truth.tot <- rbind(nbs.truth.tot, nbs.truth.i)
      }
      
      ## families who are represented in the duplicate entities, along with the
      ## number of extra representatives (duplicates that weren't deduplicated)
      dupsReps.i <- ids$famID[match(object$dupsReps$intersection[[i]], ids$requestID)]
      fams.truth$Represented[which(fams.truth$threshold == object$details$intersection[i] &
                                     fams.truth$threshType == "intersection" &
                                     fams.truth$famID %in% dupsReps.i)] <- 1
      tab <- table(dupsReps.i)
      if(max(as.numeric(names(table(tab)))) > 1){
        for(j in 2:max(as.numeric(names(table(tab))))){
          fams.truth$Extra[which(fams.truth$threshold == object$details$intersection[i] &
                                   fams.truth$threshType == "intersection" &
                                   fams.truth$famID %in% as.numeric(names(tab[which(tab == j)])))] <- j - 1
        }
      }
    }
    for(i in 1:length(object$details$greedy)){
      nbs.truth.i <- nbs.truth %>%
        mutate(threshold = object$details$greedy[i], threshType = "greedy",
               dupType = case_when(famID1 == famID2 & matchScore.g >= threshold ~ 'Called Duplicate',
                                   famID1 == famID2 & matchScore.g < threshold ~ 'Missed Duplicate',
                                   is.na(ID1) | is.na(ID2) ~ 'Missed Duplicate',
                                   matchScore.g >= threshold ~ 'False Positive',
                                   matchScore.g < threshold ~ 'True Negative'))
      nbs.truth.tot <- rbind(nbs.truth.tot, nbs.truth.i)
      
      ## families who are represented in the duplicate entities, along with the
      ## number of extra representatives (duplicates that weren't deduplicated)
      dupsReps.i <- ids$famID[match(object$dupsReps$greedy[[i]], ids$requestID)]
      fams.truth$Represented[which(fams.truth$threshold == object$details$greedy[i] &
                                     fams.truth$threshType == "greedy" &
                                     fams.truth$famID %in% dupsReps.i)] <- 1
      tab <- table(dupsReps.i)
      if(max(as.numeric(names(table(tab)))) > 1){
        for(j in 2:max(as.numeric(names(table(tab))))){
          fams.truth$Extra[which(fams.truth$threshold == object$details$greedy[i] &
                                   fams.truth$threshType == "greedy" &
                                   fams.truth$famID %in% as.numeric(names(tab[which(tab == j)])))] <- j - 1
        }
      }
    }
  } else if(object$method == "intersection"){
    for(i in 1:length(object$details$intersection)){
      nbs.truth.i <- nbs.truth %>%
        mutate(threshold = thresh[i], threshType = "intersection",
               dupType = case_when(famID1 == famID2 & matchScore >= threshold ~ 'Called Duplicate',
                                   famID1 == famID2 & matchScore < threshold ~ 'Missed Duplicate',
                                   is.na(ID1) | is.na(ID2) ~ 'Missed Duplicate',
                                   matchScore >= threshold ~ 'False Positive',
                                   matchScore < threshold ~ 'True Negative'))
      if(i == 1){
        nbs.truth.tot <- nbs.truth.i
      } else{
        nbs.truth.tot <- rbind(nbs.truth.tot, nbs.truth.i)
      }
      
      ## families who are represented in the duplicate entities, along with the
      ## number of extra representatives (duplicates that weren't deduplicated)
      dupsReps.i <- ids$famID[match(object$dupsReps[[i]], ids$requestID)]
      fams.truth$Represented[which(fams.truth$threshold == object$details$intersection[i] &
                                     fams.truth$threshType == "intersection" &
                                     fams.truth$famID %in% dupsReps.i)] <- 1
      tab <- table(dupsReps.i)
      if(max(as.numeric(names(table(tab)))) > 1){
        for(j in 2:max(as.numeric(names(table(tab))))){
          fams.truth$Extra[which(fams.truth$threshold == object$details$intersection[i] &
                                   fams.truth$threshType == "intersection" &
                                   fams.truth$famID %in% as.numeric(names(tab[which(tab == j)])))] <- j - 1
        }
      }
    }
  } else if(object$method == "greedy"){
    for(i in 1:length(object$details$greedy)){
      thresh.quant <- quantile(nbs.truth$matchScore, thresh[i])
      nbs.truth.i <- nbs.truth %>%
        mutate(threshold = thresh[i], threshType = "greedy",
               dupType = case_when(famID1 == famID2 & matchScore >= thresh.quant ~ 'Called Duplicate',
                                   famID1 == famID2 & matchScore < thresh.quant ~ 'Missed Duplicate',
                                   is.na(ID1) | is.na(ID2) ~ 'Missed Duplicate',
                                   matchScore >= thresh.quant ~ 'False Positive',
                                   matchScore < thresh.quant ~ 'True Negative'))
      if(i == 1){
        nbs.truth.tot <- nbs.truth.i
      } else{
        nbs.truth.tot <- rbind(nbs.truth.tot, nbs.truth.i)
      }
      
      ## families who are represented in the duplicate entities, along with the
      ## number of extra representatives (duplicates that weren't deduplicated)
      dupsReps.i <- ids$famID[match(object$dupsReps[[i]], ids$requestID)]
      fams.truth$Represented[which(fams.truth$threshold == object$details$greedy[i] &
                                     fams.truth$threshType == "greedy" &
                                     fams.truth$famID %in% dupsReps.i)] <- 1
      tab <- table(dupsReps.i)
      if(max(as.numeric(names(table(tab)))) > 1){
        for(j in 2:max(as.numeric(names(table(tab))))){
          fams.truth$Extra[which(fams.truth$threshold == object$details$greedy[i] &
                                   fams.truth$threshType == "greedy" &
                                   fams.truth$famID %in% as.numeric(names(tab[which(tab == j)])))] <- j - 1
        }
      }
    }
  }
  
  # fams.truth <- fams.truth %>%
  #   mutate(repResult = case_when(Duplicate == 1 & Represented == 1 ~ "Represented Duplicate",
  #                                Duplicate == 1 & Represented == 0 ~ "Unrepresented Duplicate",
  #                                Duplicate == 0 & Represented == 1 ~ "Represented Non-duplicate",
  #                                Duplicate == 0 & Represented == 0 ~ "Unrepresented Non-duplicate"))
  
  TruthObj <- Truth(object, as.data.frame(nbs.truth.tot), as.data.frame(fams.truth))
  
  return(TruthObj)
  
  # scored_tibble <- object$Neighbors %>%
  #   as_tibble() %>%
  #   separate(ID1,into=c('ID1','type1'),sep='_') %>%
  #   separate(ID2,into=c('ID2','type2'),sep='_') %>%
  #   mutate(isInjected = case_when(type1 == 'injected' | type2 == 'injected' ~ 1,
  #                                 TRUE ~ 0)) %>%
  #   select(ID1, ID2, background, matchScore, cluster, isInjected)
  # truth_tibble <- object$rawData %>% # all families with true duplicates
  #   as_tibble() %>%
  #   separate(requestId, into=c('requestId', 'injectiontype')) %>%
  #   filter(Duplicate == 1 | injectiontype == 'injected') %>%
  #   distinct(requestId, Duplicate, injectiontype) %>%
  #   rename(isDuplicate = Duplicate) %>%
  #   mutate(isInjected = case_when(injectiontype=='injected' ~ 1,
  #                                 TRUE ~ 0)) %>%
  #   mutate(isDuplicate = case_when(isDuplicate==1 | isInjected == 1 ~ 1,
  #                                  TRUE ~ 0)) %>%
  #   select(requestId, isDuplicate, isInjected)
  # merge1 <- scored_tibble %>%
  #   left_join(truth_tibble, by=c('ID1'='requestId','isInjected')) # finding out if the family with ID1 has a duplicate
  # merge2 <- scored_tibble %>%
  #   left_join(truth_tibble, by=c('ID2'='requestId','isInjected')) # finding out if the family with ID2 has a duplicate
  # scored <- bind_rows(merge1, merge2) %>%
  #   replace_na(list(isDuplicate=0)) %>%
  #   mutate(isDiscovered = 1,
  #          isCalled = as.numeric(as.numeric(matchScore) > .5)) %>%
  #   rename(isBackground = background) %>%
  #   select(-cluster)
  # not_discovered <- truth_tibble %>%
  #   filter(!(requestId %in% c(scored$ID1, scored$ID2)))
  # if (nrow(not_discovered) > 0) {
  #   scored <- not_discovered %>%
  #     rename(ID1 = requestId) %>%
  #     mutate(isBackground = 0, matchScore = 0, isDiscovered=0,isCalled=0) %>%
  #     bind_rows(scored)
  # }
  # return(scored)
}

#' Reveals
#'
#' @param object \code{\link{Scores}} object
#' @param requestID Character value for the column that identifies the apparent family ID
#' @param famID Character value for the column that identifies the true family ID
#' @param thresh Vector of thresholds at which to classify scores as duplicates. These can
#' either be integers from 1 to 7 (for the intersection score) or numbers between 0 and 1,
#' representing quantiles (for the greedy match score). If the \code{Scores} object has scores
#' for both the intersection and greedy match scores, then \code{thresh} is a list of 2 vectors,
#' where the first element is named "intersection" and the second element is named "greedy".
#' @rdname revealTruth
#' @export
revealTruth.Scores <- function(object, requestID, famID, thresh) {
  tdat1 <- object$rawData %>%
    as_tibble() %>%
    select(!!sym(famID), !!sym(requestID)) %>%
    rename(famID1 = !!sym(famID), ID1 = !!sym(requestID))
  tdat2 <- object$rawData %>%
    as_tibble() %>%
    select(!!sym(famID), !!sym(requestID)) %>%
    rename(famID2 = !!sym(famID), ID2 = !!sym(requestID))
  
  ## comparing the requestID and famID
  if(object$method == "both"){
    dups <- inner_join(object$Neighbors$intersection %>%
                         as_tibble() %>%
                         mutate(ID1 = as.numeric(ID1), ID2 = as.numeric(ID2)) %>%
                         left_join(tdat1, by = 'ID1') %>%
                         left_join(tdat2, by = 'ID2') %>%
                         rename(matchScore.i = matchScore),
                       object$Neighbors$greedy %>%
                         as_tibble() %>%
                         mutate(ID1 = as.numeric(ID1), ID2 = as.numeric(ID2)) %>%
                         left_join(tdat1, by = 'ID1') %>%
                         left_join(tdat2, by = 'ID2') %>%
                         rename(matchScore.g = matchScore),
                       by = c("ID1", "ID2", "famID1", "famID2"))
  } else{
    dups <- object$Neighbors %>%
      as_tibble() %>%
      mutate(ID1 = as.numeric(ID1), ID2 = as.numeric(ID2)) %>%
      left_join(tdat1, by = 'ID1') %>%
      left_join(tdat2, by = 'ID2')
  }
  
  # if(is.list(res.scores$Neighbors)){
  #   thresh <- list(intersection = sort(unique(as.numeric(unlist(lapply(strsplit(names(select(res.scores$Neighbors$intersection,
  #                                                                                            starts_with("guessDup"))), "_"),
  #                                                                      function(x) x[[2]]))))),
  #                  greedy = sort(unique(as.numeric(unlist(lapply(strsplit(names(select(res.scores$Neighbors$greedy,
  #                                                                                      starts_with("guessDup"))), "_"),
  #                                                                function(x) x[[2]]))))))
  # } else{
  #   thresh <- sort(unique(as.numeric(unlist(lapply(strsplit(names(select(res.scores$Neighbors,
  #                                                                        starts_with("guessDup"))), "_"),
  #                                                  function(x) x[[2]])))))
  # }
  
  ## Categorizing the neighbor as "called duplicate", "false positive", "missed duplicate",
  ## or "false negative".
  if(object$method == "both"){
    for(i in 1:length(thresh$intersection)){
      dups.i <- dups %>%
        mutate(threshold = thresh$intersection[i], threshType = "intersection",
               dupType = case_when(famID1 == famID2 & matchScore.i >= threshold ~ 'Called Duplicate',
                                   famID1 == famID2 & matchScore.i < threshold ~ 'Missed Duplicate',
                                   is.na(ID1) | is.na(ID2) ~ 'Missed Duplicate',
                                   matchScore.i >= threshold ~ 'False Positive',
                                   matchScore.i < threshold ~ 'True Negative'))
      if(i == 1){
        dups.tot <- dups.i
      } else{
        dups.tot <- rbind(dups.tot, dups.i)
      }
    }
    for(i in 1:length(thresh$greedy)){
      dups.i <- dups %>%
        mutate(threshold = thresh$greedy[i], threshType = "greedy",
               dupType = case_when(famID1 == famID2 & matchScore.g >= threshold ~ 'Called Duplicate',
                                   famID1 == famID2 & matchScore.g < threshold ~ 'Missed Duplicate',
                                   is.na(ID1) | is.na(ID2) ~ 'Missed Duplicate',
                                   matchScore.g >= threshold ~ 'False Positive',
                                   matchScore.g < threshold ~ 'True Negative'))
      dups.tot <- rbind(dups.tot, dups.i)
    }
  } else if(object$method == "intersection"){
    for(i in 1:length(thresh)){
      dups.i <- dups %>%
        mutate(threshold = thresh[i], threshType = "intersection",
               dupType = case_when(famID1 == famID2 & matchScore >= threshold ~ 'Called Duplicate',
                                   famID1 == famID2 & matchScore < threshold ~ 'Missed Duplicate',
                                   is.na(ID1) | is.na(ID2) ~ 'Missed Duplicate',
                                   matchScore >= threshold ~ 'False Positive',
                                   matchScore < threshold ~ 'True Negative'))
      if(i == 1){
        dups.tot <- dups.i
      } else{
        dups.tot <- rbind(dups.tot, dups.i)
      }
    }
  } else if(object$method == "greedy"){
    for(i in 1:length(thresh)){
      thresh.quant <- quantile(dups$matchScore, thresh[i])
      dups.i <- dups %>%
        mutate(threshold = thresh[i], threshType = "greedy",
               dupType = case_when(famID1 == famID2 & matchScore >= thresh.quant ~ 'Called Duplicate',
                                   famID1 == famID2 & matchScore < thresh.quant ~ 'Missed Duplicate',
                                   is.na(ID1) | is.na(ID2) ~ 'Missed Duplicate',
                                   matchScore >= thresh.quant ~ 'False Positive',
                                   matchScore < thresh.quant ~ 'True Negative'))
      if(i == 1){
        dups.tot <- dups.i
      } else{
        dups.tot <- rbind(dups.tot, dups.i)
      }
    }
  }
  return(dups.tot)
}
