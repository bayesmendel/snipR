#' deDuplicate scored neighbors matrix.
#' 
#' Outputs the \code{\link{Duplicates}} object, providing the duplicate entities and
#' the representatives for each duplicate entity.
#'
#' @param object \code{\link{Scores}} object containing pre-computed Neighbors and match scores.
#' @param thresh Vector of thresholds at which to classify scores as duplicates. These can
#' either be integers from 1 to 7 (for the intersection score) or numbers between 0 and 1,
#' representing quantiles (for the greedy match score).
#' @param priority A list of structure (var = 'Varx', min = TRUE) with 'Varx'
#' being a character value corresponding to a column in rawData.
#' This parameter determines how to sort the duplicates.
#' If \code{min = TRUE}, then we use the minimum value of 'Varx' for each duplicate entity.
#' Otherwise, we use the maximum value.
#' @param isProband A character value indicating the column indicator column that denotes the proband.
#' @param dateFormat Character string of the format of the date. This is only used
#' if the priority variable is a date. The format should match the formats of
#' class \code{POSIXlt} used in the \code{base::strptime} function.
#' @param requestID Column that has the ID for the family.
#' @return An object of class \code{\link{Duplicates}} containing the duplicate entities
#' and representatives for each duplicate entity (including singletons without duplicates)..
#'
#' @export
deDuplicate <- function(object, thresh, priority, isProband,
                        dateFormat = NULL, requestID) {
  UseMethod("deDuplicate", object)
}

#' @rdname deDuplicate
#' @export
deDuplicate.default <- function(object, thresh, priority, isProband,
                                dateFormat = NULL, requestID) {
  print("deDuplication allowable on Scores objects only.")
  return(NULL)
}

#' @rdname deDuplicate
#' @export
deDuplicate.Scores <- function(object, thresh, priority, isProband,
                               dateFormat = NULL, requestID) {
  # probands
  dat.pro <- object$rawData %>%
    filter(!!sym(isProband) == 1)
  
  if(object$method == "both"){ ## if there are both intersection and greedy match scores
    
    dupsList <- dupsReps <- setNames(vector("list", 2), c("intersection", "greedy"))
    dupsList$intersection <- dupsReps$intersection <- setNames(vector("list", length(thresh$intersection)),
                                                               paste0("thresh_", thresh$intersection))
    dupsList$greedy <- dupsReps$greedy <- setNames(vector("list", length(thresh$greedy)),
                                                   paste0("thresh_", thresh$greedy))
    
    for(i in 1:length(thresh$intersection)){
      # neighbors with the match score
      nbs <- object$Neighbors$intersection %>%
        filter(matchScore >= thresh$intersection[i]) %>%
        transmute(ID1, ID2)
      
      if (nrow(nbs)==0) {
        print(paste0('0 duplicates found for threshold ', thresh[i], '. Try lowering the threshold.'))
        dupsList$intersection[[i]] <- as.list(unique(object$rawData[, eval(requestID)]))
        names(dupsList$intersection[[i]]) <- 1:length(dupsList$intersection[[i]])
        dupsReps$intersection[[i]] <- unique(object$rawData[, eval(requestID)])
      } else{
        
        # adding extra variables that include both the request ID and the priority variable
        nbs <- nbs %>% mutate(idVar1 = paste(ID1, dat.pro[match(ID1, dat.pro[, object$ID]), priority$var],
                                             sep = "_"),
                              idVar2 = paste(ID2, dat.pro[match(ID2, dat.pro[, object$ID]), priority$var],
                                             sep = "_"))
        
        ## Obtaining duplicate entities (equivalence classes)
        dups <- igraph::components(igraph::graph_from_data_frame(nbs %>% select(idVar1, idVar2)))
        members <- dups$membership
        id.members <- as.numeric(unlist(strsplit(names(members), "_"))[seq(1, 2*length(members), 2)])
        dupsList$intersection[[i]] <- tapply(id.members, members, sort)
        famids.singleton <- unique(setdiff(object$rawData %>%
                                             filter(!!sym(isProband) == 1) %>%
                                             select(!!sym(object$ID)) %>%
                                             unlist() %>%
                                             as.numeric(),
                                           unlist(dupsList$intersection[[i]])))
        dupsList$intersection[[i]] <- c(dupsList$intersection[[i]], famids.singleton)
        
        # Obtaining the representative for each duplicate entity
        if(!is.null(dateFormat)){
          dupsReps$intersection[[i]] <- tapply(names(members), members, function(x){
            dates <- strptime(unlist(strsplit(x, "_"))[seq(2, 2*length(x), 2)],
                              dateFormat)
            return(as.numeric(unlist(strsplit(x, "_"))[seq(1, 2*length(x), 2)])[
              which(dates == min(dates))[1]
            ])
          })
        } else{
          dupsReps$intersection[[i]] <- tapply(names(members), members, function(x){
            var.priority <- unlist(strsplit(x, "_"))[seq(2, 2*length(x), 2)]
            return(as.numeric(unlist(strsplit(x, "_"))[seq(1, 2*length(x), 2)])[
              which(var.priority == min(var.priority))[1]
            ])
          })
        }
        dupsReps$intersection[[i]] <- as.numeric(unlist(dupsReps$intersection[[i]]))
        dupsReps$intersection[[i]] <- c(dupsReps$intersection[[i]], famids.singleton)
      }
    }
    
    for(i in 1:length(thresh$greedy)){
      # neighbors with the match score
      nbs <- object$Neighbors$greedy %>%
        filter(matchScore >= quantile(matchScore, thresh$greedy[i])) %>%
        transmute(ID1, ID2)
      
      if (nrow(nbs)==0) {
        print(paste0('0 duplicates found for threshold ', thresh[i], '. Try lowering the threshold.'))
        dupsList$greedy[[i]] <- as.list(unique(object$rawData[, eval(requestID)]))
        names(dupsList$greedy[[i]]) <- 1:length(dupsList$greedy[[i]])
        dupsReps$greedy[[i]] <- unique(object$rawData[, eval(requestID)])
      } else{
        
        # adding extra variables that include both the request ID and the priority variable
        nbs <- nbs %>% mutate(idVar1 = paste(ID1, dat.pro[match(ID1, dat.pro[, object$ID]), priority$var],
                                             sep = "_"),
                              idVar2 = paste(ID2, dat.pro[match(ID2, dat.pro[, object$ID]), priority$var],
                                             sep = "_"))
        
        ## Obtaining duplicate entities (equivalence classes)
        dups <- igraph::components(igraph::graph_from_data_frame(nbs %>% select(idVar1, idVar2)))
        members <- dups$membership
        id.members <- as.numeric(unlist(strsplit(names(members), "_"))[seq(1, 2*length(members), 2)])
        dupsList$greedy[[i]] <- tapply(id.members, members, sort)
        famids.singleton <- unique(setdiff(object$rawData %>%
                                             filter(!!sym(isProband) == 1) %>%
                                             select(!!sym(object$ID)) %>%
                                             unlist() %>%
                                             as.numeric(),
                                           unlist(dupsList$greedy[[i]])))
        dupsList$greedy[[i]] <- c(dupsList$greedy[[i]], famids.singleton)
        
        # Obtaining the representative for each duplicate entity
        if(!is.null(dateFormat)){
          dupsReps$greedy[[i]] <- tapply(names(members), members, function(x){
            dates <- strptime(unlist(strsplit(x, "_"))[seq(2, 2*length(x), 2)],
                              dateFormat)
            return(as.numeric(unlist(strsplit(x, "_"))[seq(1, 2*length(x), 2)])[
              which(dates == min(dates))[1]
            ])
          })
        } else{
          dupsReps$greedy[[i]] <- tapply(names(members), members, function(x){
            var.priority <- unlist(strsplit(x, "_"))[seq(2, 2*length(x), 2)]
            return(as.numeric(unlist(strsplit(x, "_"))[seq(1, 2*length(x), 2)])[
              which(var.priority == min(var.priority))[1]
            ])
          })
        }
        dupsReps$greedy[[i]] <- as.numeric(unlist(dupsReps$greedy[[i]]))
        dupsReps$greedy[[i]] <- c(dupsReps$greedy[[i]], famids.singleton)
      }
    }
  } else{
    dupsList <- dupsReps <- setNames(vector("list", length(thresh)), paste0("thresh_", thresh))
    
    for(i in 1:length(thresh)){
      # neighbors with the match score
      if(thresh[i] %in% 1:7){
        nbs <- object$Neighbors %>%
          filter(matchScore >= thresh[i]) %>%
          transmute(ID1, ID2)
      } else if(thresh[i] > 0 & thresh[i] < 1){
        nbs <- object$Neighbors %>%
          filter(matchScore >= quantile(matchScore, thresh[i])) %>%
          transmute(ID1, ID2)
      }
      
      if (nrow(nbs)==0) {
        print(paste0('0 duplicates found for threshold ', thresh[i], '. Try lowering the threshold.'))
        dupsList[[i]] <- as.list(unique(object$rawData[, eval(requestID)]))
        names(dupsList[[i]]) <- 1:length(dupsList[[i]])
        dupsReps[[i]] <- unique(object$rawData[, eval(requestID)])
      } else{
        
        
        # adding extra variables that include both the request ID and the priority variable
        nbs <- nbs %>% mutate(idVar1 = paste(ID1, dat.pro[match(ID1, dat.pro[, object$ID]), priority$var],
                                             sep = "_"),
                              idVar2 = paste(ID2, dat.pro[match(ID2, dat.pro[, object$ID]), priority$var],
                                             sep = "_"))
        
        ## Obtaining duplicate entities (equivalence classes)
        dups <- igraph::components(igraph::graph_from_data_frame(nbs %>% select(idVar1, idVar2)))
        members <- dups$membership
        id.members <- as.numeric(unlist(strsplit(names(members), "_"))[seq(1, 2*length(members), 2)])
        dupsList[[i]] <- tapply(id.members, members, sort)
        famids.singleton <- unique(setdiff(object$rawData %>%
                                             filter(!!sym(isProband) == 1) %>%
                                             select(!!sym(object$ID)) %>%
                                             unlist() %>%
                                             as.numeric(),
                                           unlist(dupsList[[i]])))
        dupsList[[i]] <- c(dupsList[[i]], famids.singleton)
        
        # Obtaining the representative for each duplicate entity
        if(!is.null(dateFormat)){
          dupsReps[[i]] <- tapply(names(members), members, function(x){
            dates <- strptime(unlist(strsplit(x, "_"))[seq(2, 2*length(x), 2)],
                              dateFormat)
            return(as.numeric(unlist(strsplit(x, "_"))[seq(1, 2*length(x), 2)])[
              which(dates == min(dates))[1]
            ])
          })
        } else{
          dupsReps[[i]] <- tapply(names(members), members, function(x){
            var.priority <- unlist(strsplit(x, "_"))[seq(2, 2*length(x), 2)]
            return(as.numeric(unlist(strsplit(x, "_"))[seq(1, 2*length(x), 2)])[
              which(var.priority == min(var.priority))[1]
            ])
          })
        }
        dupsReps[[i]] <- as.numeric(unlist(dupsReps[[i]]))
        dupsReps[[i]] <- c(dupsReps[[i]], famids.singleton)
      }
    }
  }
  
  ## Duplicates object
  if(object$method == "both"){
    object <- Duplicates(object, dupsList, dupsReps,
                         c(thresh, priority))
  } else if(object$method == "intersection"){
    object <- Duplicates(object, dupsList, dupsReps,
                         c(list(intersection = thresh), priority))
  } else{
    object <- Duplicates(object, dupsList, dupsReps,
                         c(list(greedy = thresh), priority))
  }
  return(object)
}
