#' Run deduplication algorithm, separated by relative type.
#' 
#' This only runs a subset of the entire algorithm, only looking at one specific
#' relative type out of the 7. It also stops after obtaining the neighbors, and 
#' doesn't score them or split into duplicate entities. This is for the purpose of running the
#' algorithm in parallel.
#'
#' @param pedigrees Pedigree data to deduplicate.
#' @param relative Relative type
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
#' @param printRuntime If TRUE, will print the runtime
#' @param seed Seed
#' @return An object of class \code{\link{Blocks}} containing the scored neighbor pairs.
#'
#' @export

bsnRelative <- function(pedigrees, relative, requestID, isProband, keyVars,
                        keyVars.male = NULL, keyVars.female = NULL,
                        keyWt = NULL, blockVar = NULL, repSN = 1, windowSN = 2, keyLength,
                        printRuntime = TRUE, seed = NULL){
  
  start <- Sys.time()
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
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
  
  if(relative == "self"){
    rel <- pedigrees %>%
      filter(!!sym(isProband) == 1) %>% 
      as.data.frame()
    vars <- keyVars
  } else if(relative == "mother"){
    rel <- pedigrees %>% 
      group_by(!!sym(requestID)) %>%
      arrange(desc(!!sym(isProband)), .by_group=TRUE) %>%
      filter(ID==MotherID[1]) %>%
      ungroup %>% as.data.frame()
    vars <- vars.female
  } else if(relative == "father"){
    rel <- pedigrees %>% 
      group_by(!!sym(requestID)) %>%
      arrange(desc(!!sym(isProband)), .by_group=TRUE) %>%
      filter(ID==FatherID[1]) %>%
      ungroup %>% as.data.frame()
    vars <- vars.male
  } else if(relative == "mothersmother"){
    rel <- pedigrees %>% 
      group_by(!!sym(requestID)) %>%
      arrange(desc(!!sym(isProband)), .by_group=TRUE) %>% 
      mutate(tempID = MotherID[1]) %>% 
      filter(n()>1) %>%
      filter(ID == MotherID[ID==tempID]) %>% 
      ungroup %>% as.data.frame()
    vars <- vars.female
  } else if(relative == "mothersfather"){
    rel <- pedigrees %>% 
      group_by(!!sym(requestID)) %>%
      arrange(desc(!!sym(isProband)), .by_group=TRUE) %>% 
      mutate(tempID = MotherID[1]) %>% 
      filter(n()>1) %>%
      filter(ID == FatherID[ID==tempID]) %>% 
      ungroup %>% as.data.frame()
    vars <- vars.male
  } else if(relative == "fathersmother"){
    rel <- pedigrees %>% 
      group_by(!!sym(requestID)) %>%
      arrange(desc(!!sym(isProband)), .by_group=TRUE) %>% 
      mutate(tempID = FatherID[1]) %>% 
      filter(n()>1) %>%
      filter(ID == MotherID[ID==tempID]) %>% 
      ungroup %>% as.data.frame()
    vars <- vars.female
  } else if(relative == "fathersfather"){
    rel <- pedigrees %>% 
      group_by(!!sym(requestID)) %>%
      arrange(desc(!!sym(isProband)), .by_group=TRUE) %>% 
      mutate(tempID = FatherID[1]) %>% 
      filter(n()>1) %>%
      filter(ID == FatherID[ID==tempID]) %>% 
      ungroup %>% as.data.frame()
    vars <- vars.male
  }
  
  ## Initialize a Neighbors object
  neighbsObj <- Neighbors(rel, ID = requestID, keyVars = vars)
  
  # perform an iteration of algorithm, specifying the number of repetitions and window size
  blocksObj <- blockedSN(neighbsObj, blockVar = blockVar,
                         repSN = repSN, windowSN = windowSN,
                         keyLength = replace(keyLength, keyLength > nrow(neighbsObj$keyVars),
                                             nrow(neighbsObj$keyVars)))
  
  if(printRuntime){
    print(difftime(Sys.time(), start, units = "secs"))
  }
  
  return(blocksObj)
}
