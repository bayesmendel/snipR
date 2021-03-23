#' Run entire deduplication algorithm
#'
#' @param pedigrees Pedigree data to deduplicate.
#' @param requestID Column that has the ID for the family.
#' @param isProband Column that indicates the proband.
#' @param keyVars Character vector of column names for the variables in the sort key
#' @param keyVars.male Optional character vector of column names for the variables in the sort key
#' that are specific to males
#' @param keyVars.female Optional character vector of column names for the variables in the sort key
#' that are specific to females
#' @param keyWt Numeric vector of weights assigned to variables in the sort key, corresponding
#' to \code{keyVars}. If \code{NULL}, the standard deviations of the
#' variables in the data will be used as weights.
#' @param blockVar Vector of column names for the blocking variables, where families in different
#' blocks will not be considered when searching for duplicates.
#' @param repSN Number of iterations when sorting neighbors according to the sort key
#' @param windowSN Integer representing the size of the sliding window to use during sorted neighbors.
#' @param keyLength Numeric vector representing the number of key variables
#' (out of \code{c(keyVar.bin, keyVar.cont)}) to concatenate per sort key. If missing, the key
#' lengths will be randomly generated. The length of the vector should be \code{repSN}.
#' @param method If "intersection", we use the intersection score. If "greedy", we use the
#' greedy match score. If "both", we use both.
#' @param thresh Vector of thresholds. If \code{method = "intersection"}, then a pair is
#' considered neighbors if the intersection match score is greater than the threshold. If
#' \code{method = "greedy"}, then the threshold is treated as a percentile, and a pair is
#' considered neighbors if the greedy match score is greater than the percentile. If
#' \code{method = "both"}, then the user should provide a list, such as
#' \code{list(intersection = 1:7, threshold = c(0.8, 0.9))}.
#' @param priority A list of structure (var = 'Varx', min = TRUE) with 'Varx'
#' being a character value corresponding to a column in rawData.
#' This parameter determines how to sort the duplicates.
#' If \code{min = TRUE}, then we use the minimum value of 'Varx' for each duplicate entity.
#' Otherwise, we use the maximum value.
#' @param dateFormat Character string of the format of the date. This is only used
#' if the priority variable is a date. The format should match the formats of
#' class \code{POSIXlt} used in the \code{base::strptime} function.
#' @param printRuntime If TRUE, will print the runtime
#' @param seed Seed
#' @return An object of class \code{\link{Duplicates}} containing the duplicate entities
#' and representatives for each duplicate entity (including singletons without duplicates).
#'
#' @export

bsn <- function(pedigrees, requestID, isProband, keyVars,
                keyVars.male = NULL, keyVars.female = NULL,
                keyWt = NULL, blockVar = NULL, repSN = 1, windowSN = 10,
                keyLength = length(keyVars), method = "intersection",
                thresh = 1:7, priority, dateFormat = NULL,
                printRuntime = TRUE, seed = NULL){
  
  start <- Sys.time()
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  self <- pedigrees %>%
    filter(!!sym(isProband) == 1) %>% 
    as.data.frame()
  mothers <- pedigrees %>% 
    group_by(!!sym(requestID)) %>%
    arrange(desc(!!sym(isProband)), .by_group=TRUE) %>%
    filter(ID==MotherID[1]) %>%
    ungroup %>% as.data.frame()
  fathers <- pedigrees %>% 
    group_by(!!sym(requestID)) %>%
    arrange(desc(!!sym(isProband)), .by_group=TRUE) %>%
    filter(ID==FatherID[1]) %>%
    ungroup %>% as.data.frame()
  mothersmothers <- pedigrees %>% 
    group_by(!!sym(requestID)) %>%
    arrange(desc(!!sym(isProband)), .by_group=TRUE) %>% 
    mutate(tempID = MotherID[1]) %>% 
    filter(n()>1) %>%
    filter(ID == MotherID[ID==tempID]) %>% 
    ungroup %>% as.data.frame()
  mothersfathers <- pedigrees %>% 
    group_by(!!sym(requestID)) %>%
    arrange(desc(!!sym(isProband)), .by_group=TRUE) %>% 
    mutate(tempID = MotherID[1]) %>% 
    filter(n()>1) %>%
    filter(ID == FatherID[ID==tempID]) %>% 
    ungroup %>% as.data.frame()
  fathersmothers <- pedigrees %>% 
    group_by(!!sym(requestID)) %>%
    arrange(desc(!!sym(isProband)), .by_group=TRUE) %>% 
    mutate(tempID = FatherID[1]) %>% 
    filter(n()>1) %>%
    filter(ID == MotherID[ID==tempID]) %>% 
    ungroup %>% as.data.frame()
  fathersfathers <- pedigrees %>% 
    group_by(!!sym(requestID)) %>%
    arrange(desc(!!sym(isProband)), .by_group=TRUE) %>% 
    mutate(tempID = FatherID[1]) %>% 
    filter(n()>1) %>%
    filter(ID == FatherID[ID==tempID]) %>% 
    ungroup %>% as.data.frame()
  
  ## Set key variables
  
  # Define the key variables below for the sorted neighbors. These should correspond to columns in the dataset under consideration for which we eventually wish to include in the duplicate pair match scoring.
  # variables to be used in generating the sort keys
  # keyLength.female <- keyLength.male <- keyLength
  if(!is.null(keyVars.male)){
    vars.female <- setdiff(keyVars, keyVars.male)
    # keyLength.female[keyLength.female > length(vars.female)] <- length(vars.female)
  } else{
    vars.female <- keyVars
  }
  # ind.female <- match(vars.female, keyVars)
  if(!is.null(keyVars.female)){
    vars.male <- setdiff(keyVars, keyVars.female)
    # keyLength.male[keyLength.male > length(vars.male)] <- length(vars.male)
  } else{
    vars.male <- keyVars
  }
  # ind.male <- match(vars.male, keyVars)
  # vars <- c(keyVar.bin, keyVar.cont, blockVar)
  
  ## Weight each variable in the scoring by its variability
  # if(is.null(keyWt)){
  #   wt <- rep(NA, length(vars))
  #   for (i in 1:length(vars)) {
  #     wt[i] <- sd(as.numeric(as.character(pedigrees[, vars[i]])), na.rm=TRUE)
  #     if (is.na(wt[i])) {wt[i] <- 1}
  #     # selfs.neighbsObj[["keyVars"]][kvar,"keyWt"] <- wt
  #   }
  # }
  # 
  # if(!is.null(keyWt)){
  #   wt.female <- keyWt[match(vars.female, keyVars)]
  #   wt.male <- keyWt[match(vars.male, keyVars)]
  #   wt.all <- keyWt
  # } else{
  #   wt.all <- rep(NA, length(keyVars))
  #   for (i in 1:length(keyVars)) {
  #     wt.all[i] <- sd(as.numeric(as.character(pedigrees[, keyVars[i]])), na.rm = TRUE)
  #     if (is.na(wt.all[i])) {wt.all[i] <- 1}
  #   }
  #   wt.male <- wt.all[ind.male]
  #   wt.female <- wt.all[ind.female]
  # }
  
  ## Initialize a Neighbors object
  selfs.neighbsObj <- Neighbors(self, ID = requestID, keyVars = keyVars)
  # repeat for all family members
  mothers.neighbsObj <- Neighbors(mothers, ID = requestID, keyVars = vars.female)
  fathers.neighbsObj <- Neighbors(fathers, ID = requestID, keyVars = vars.male)
  mothersmothers.neighbsObj <- Neighbors(mothersmothers, ID = requestID, keyVars = vars.female)
  mothersfathers.neighbsObj <- Neighbors(mothersfathers, ID = requestID, keyVars = vars.male)
  fathersmothers.neighbsObj <- Neighbors(fathersmothers, ID = requestID, keyVars = vars.female)
  fathersfathers.neighbsObj <- Neighbors(fathersfathers, ID = requestID, keyVars = vars.male)
  
  # perform an iteration of algorithm, specifying the number of repetitions and window size
  selfs.blocksObj <- blockedSN(selfs.neighbsObj, blockVar = blockVar,
                               repSN = repSN, windowSN = windowSN,
                               keyLength = replace(keyLength, keyLength > nrow(selfs.neighbsObj$keyVars),
                                                   nrow(selfs.neighbsObj$keyVars)))
  # repeat for all family members
  mothers.blocksObj <- blockedSN(mothers.neighbsObj, blockVar = blockVar,
                                 repSN = repSN, windowSN = windowSN,
                                 keyLength = replace(keyLength, keyLength > nrow(mothers.neighbsObj$keyVars),
                                                     nrow(mothers.neighbsObj$keyVars)))
  fathers.blocksObj <- blockedSN(fathers.neighbsObj, blockVar = blockVar,
                                 repSN = repSN, windowSN = windowSN,
                                 keyLength = replace(keyLength, keyLength > nrow(fathers.neighbsObj$keyVars),
                                                     nrow(fathers.neighbsObj$keyVars)))
  mothersmothers.blocksObj <- blockedSN(mothersmothers.neighbsObj, blockVar = blockVar,
                                        repSN = repSN, windowSN = windowSN,
                                        keyLength = replace(keyLength, keyLength > nrow(mothersmothers.neighbsObj$keyVars),
                                                            nrow(mothersmothers.neighbsObj$keyVars)))
  mothersfathers.blocksObj <- blockedSN(mothersfathers.neighbsObj, blockVar = blockVar,
                                        repSN = repSN, windowSN = windowSN,
                                        keyLength = replace(keyLength, keyLength > nrow(mothersfathers.neighbsObj$keyVars),
                                                            nrow(mothersfathers.neighbsObj$keyVars)))
  fathersmothers.blocksObj <- blockedSN(fathersmothers.neighbsObj, blockVar = blockVar,
                                        repSN = repSN, windowSN = windowSN,
                                        keyLength = replace(keyLength, keyLength > nrow(fathersmothers.neighbsObj$keyVars),
                                                            nrow(fathersmothers.neighbsObj$keyVars)))
  fathersfathers.blocksObj <- blockedSN(fathersfathers.neighbsObj, blockVar = blockVar,
                                        repSN = repSN, windowSN = windowSN,
                                        keyLength = replace(keyLength, keyLength > nrow(fathersfathers.neighbsObj$keyVars),
                                                            nrow(fathersfathers.neighbsObj$keyVars)))
  
  # ## Generating background pairs
  # selfs.blocksObj <- sampleBackground(selfs.blocksObj, num = backgroundProp)
  # 
  ## We may leverage the pedigree information by intersecting the found pairs across defined family member types (father, mother, etc.).
  self.pairs <- cbind(selfs.blocksObj$Neighbors, member = "self")
  mothers.pairs <- cbind(mothers.blocksObj$Neighbors, member = "mothers")
  fathers.pairs <- cbind(fathers.blocksObj$Neighbors, member = "fathers")
  mothersmothers.pairs <- cbind(mothersmothers.blocksObj$Neighbors, member = "mothersmothers")
  mothersfathers.pairs <- cbind(mothersfathers.blocksObj$Neighbors, member = "mothersfathers")
  fathersmothers.pairs <- cbind(fathersmothers.blocksObj$Neighbors, member = "fathersmothers")
  fathersfathers.pairs <- cbind(fathersfathers.blocksObj$Neighbors, member = "fathersfathers")
  all.pairs <- rbind(self.pairs, mothers.pairs, fathers.pairs, mothersmothers.pairs,
                     mothersfathers.pairs, fathersmothers.pairs, fathersfathers.pairs)
  unique.pairs <- as_tibble(all.pairs) %>% 
    group_by(ID1, ID2) %>%
    mutate(found = n()) %>%
    filter(row_number() == 1) %>% 
    select(ID1, ID2, background, found) %>%
    ungroup() %>%
    as.data.frame()
  
  blocksObj <- selfs.blocksObj
  blocksObj$Neighbors <- as.matrix(unique.pairs)
  
  # get match score -- either intersection score or greedy match score
  if(method == "intersection"){
    scoresObj <- scoreNeighbors(blocksObj, method = "intersection")
    nbs <- scoresObj$Neighbors[, c("ID1", "ID2", "matchScore")]
    nbs <- as.data.frame(nbs)
    nbs <- nbs %>% mutate(matchScore = as.numeric(matchScore))
    
    # for(i in 1:length(thresh)){
    #   nbs$temp <- 0
    #   nbs$temp[nbs$matchScore >= thresh[i]] <- 1
    #   names(nbs)[which(names(nbs) == "temp")] <- paste0("guessDup_", thresh[i])
    # }
    scoresObj$Neighbors <- nbs
  } else if(method == "greedy"){
    scoresObj <- scoreNeighbors(blocksObj, method = "greedy")
    nbs <- scoresObj$Neighbors[, c("ID1", "ID2", "matchScore")]
    nbs <- as.data.frame(nbs)
    nbs <- nbs %>% mutate(matchScore = as.numeric(matchScore))
    
    # for(i in 1:length(thresh)){
    #   nbs$temp <- 0
    #   nbs$temp[nbs$matchScore >= quantile(nbs$matchScore, thresh[i])] <- 1
    #   names(nbs)[which(names(nbs) == "temp")] <- paste0("guessDup_", thresh[i])
    # }
    scoresObj$Neighbors <- nbs
  } else if(method == "both"){
    scoresObj.i <- scoreNeighbors(blocksObj, method = "intersection")
    nbs.i <- scoresObj.i$Neighbors[, c("ID1", "ID2", "matchScore")]
    nbs.i <- as.data.frame(nbs.i)
    nbs.i <- nbs.i %>% mutate(matchScore = as.numeric(matchScore))
    scoresObj.g <- scoreNeighbors(blocksObj, method = "greedy")
    nbs.g <- scoresObj.g$Neighbors[, c("ID1", "ID2", "matchScore")]
    nbs.g <- as.data.frame(nbs.g)
    nbs.g <- nbs.g %>% mutate(matchScore = as.numeric(matchScore))
    
    # for(i in 1:length(thresh$intersection)){
    #   nbs.i$temp <- 0
    #   nbs.i$temp[nbs.i$matchScore >= thresh$intersection[i]] <- 1
    #   names(nbs.i)[which(names(nbs.i) == "temp")] <- paste0("guessDup_", thresh$intersection[i])
    # }
    # for(i in 1:length(thresh$greedy)){
    #   nbs.g$temp <- 0
    #   nbs.g$temp[nbs.g$matchScore >= quantile(nbs.g$matchScore, thresh$greedy[i])] <- 1
    #   names(nbs.g)[which(names(nbs.g) == "temp")] <- paste0("guessDup_", thresh$greedy[i])
    # }
    scoresObj <- scoresObj.i
    scoresObj$Neighbors <- list(intersection = nbs.i, greedy = nbs.g)
    scoresObj$method <- "both"
  }
  
  ## Deduplicate neighbors
  dupObj <- deDuplicate(object = scoresObj, thresh = thresh, priority = priority,
                        isProband = isProband, dateFormat = dateFormat,
                        requestID = requestID)
  
  if(printRuntime){
    print(difftime(Sys.time(), start, units = "secs"))
  }
  
  return(dupObj)
}
