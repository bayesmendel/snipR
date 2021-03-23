#' Blocked sorted neighbors generic.
#'
#' \code{blockedSN} returns a \code{\link{Blocks}} object created from a single
#' iteration of the blocked sorted neighbors algorithm.
#'
#' @param object \code{\link{Neighbors}} or \code{\link{Blocks}} object to perform BSN algorithm on.
#' @return An object of class \code{\link{Blocks}} containing the neighbors found
#' and keys used during the blocked sorted neighbors iteration. 
#'
#' @export
blockedSN <- function(object, blockVar, repSN, windowSN, keyLength) {
    UseMethod("blockedSN", object)
}

#' @rdname blockedSN
#' @export 
blockedSN.default <- function(object, blockVar, repSN, windowSN, keyLength) {
    print("Please first initialize a Neighbors object.")
    return(NULL)
}

#' @rdname blockedSN
#' @param blockVar Character value.  Blocking variable.
#' @param repSN Integer value. How much iterations of sorted neighbors to perform?
#' @param windowSN Integer value. Size of sliding window to use during sorted neighbors.
#' @param keyLength Numeric value or vector. How many keyVars to concatenate per sort key.
#' @export 
blockedSN.Neighbors <- function(object, blockVar, repSN, windowSN, keyLength) {
    keyVars <- object[["keyVars"]][, "keyVars"]
    keyVars <- keyVars[!(keyVars %in% blockVar)]
    keyWt <- as.numeric(object[["keyVars"]][!(object[["keyVars"]][, "keyVars"] %in% blockVar), "keyWt"])
    # generate pseudo random key lengths if none provided
    if (missing(keyLength)) {
        stats::runif(1)
        keyLength <- abs(.Random.seed[1:repSN])%%10
        idx <- keyLength > length(keyVars)
        keyLength[idx] <- length(keyVars)
    }
    if (length(keyLength) != repSN) {
        keyLength <- rep_len(keyLength, repSN)
    }
    # generate keys to sort on
    keys <- lapply(1:repSN, function(sn) {
        keyVars[sample(length(keyVars), keyLength, prob = keyWt)]
    })
    # get unique values to block on
    
    if(!is.null(blockVar)){
        blockVec <- object[["rawData"]][, blockVar]
        iterval <- unique(blockVec)
        iter <- 1:length(iterval)
        # partition input data to separate blocks
        subdat <- lapply(iter, function(it) {
            object[["rawData"]][blockVec == iterval[it], ]
        })
    } else{
        subdat <- list(object[["rawData"]])
        iter <- 1
    }
    names(subdat) <- iter
    rez <- lapply(iter, function(it) {
        #print(it)
        print(paste0("iter: ", it, " out of ", length(iter)))
        if (nrow(subdat[[it]])==1) {return(NULL)}
        if (nrow(subdat[[it]])==2) {
            nei <- cbind(subdat[[it]][1,object[["ID"]]],subdat[[it]][2,object[["ID"]]])
            colnames(nei) <- c("sorted.ids", "")
            return(nei)
        }
        nei <- lapply(1:repSN, function(sn) {
            if (nrow(subdat[[it]]) > windowSN) {
                sortedNeighbors(subdat[[it]], keys[[sn]], windowSN = windowSN, ID = object[["ID"]])
            } else {
                sortedNeighbors(subdat[[it]], keys[[sn]], windowSN = nrow(subdat[[it]])-1, ID = object[["ID"]])
            }
        })
        do.call(rbind, nei)
    })
    rez <- do.call(rbind, rez)
    rez <- t(apply(rez, 1, sort))
    rez <- unique(rez)
    rez <- cbind(rez, rep(0, nrow(rez)))
    colnames(rez) <- c("ID1", "ID2", "background")
    # promote to Blocks object
    object <- Blocks(object, rez, list(blockVar = blockVar, sortKeys = keys))
    return(object)
}

#' @rdname blockedSN
#' @export 
blockedSN.Blocks <- function(object, blockVar, repSN, windowSN, keyLength) {
    keyVars <- object[["keyVars"]][, "keyVars"]
    keyVars <- keyVars[keyVars != blockVar]
    # generate pseudo random key lengths if none provided
    if (missing(keyLength)) {
        stats::runif(1)
        keyLength <- abs(.Random.seed[1:repSN])%%10
        idx <- keyLength > length(keyVars)
        keyLength[idx] <- length(keyVars)
    }
    if (length(keyLength != repSN)) {
        keyLength <- rep_len(keyLength, repSN)
    }
    # generate keys to sort on
    keys <- lapply(1:repSN, function(sn) {
        keyVars[sample(length(keyVars), keyLength)]
    })
    # get unique values to block on
    blockVec <- object[["rawData"]][, blockVar]
    iterval <- unique(blockVec)
    iter <- 1:length(iterval)
    # partition input data to separate blocks
    subdat <- lapply(iter, function(it) {
        object[["rawData"]][blockVec == iterval[it], ]
    })
    names(subdat) <- iter
    rez <- lapply(iter, function(it) {
        if (nrow(subdat[[it]])==1) {return(NULL)}
        if (nrow(subdat[[it]])==2) {
            nei <- cbind(subdat[[it]][1,object[["ID"]]],subdat[[it]][2,object[["ID"]]])
            colnames(nei) <- c("sorted.ids", "")
            return(nei)
        }
        nei <- lapply(1:repSN, function(sn) {
            if (nrow(subdat[[it]]) > windowSN) {
                sortedNeighbors(subdat[[it]], keys[[sn]], windowSN = windowSN, ID = object[["ID"]])
            } else {
                sortedNeighbors(subdat[[it]], keys[[sn]], windowSN = nrow(subdat[[it]])-1, ID = object[["ID"]])
            }
        })
        do.call(rbind, nei)
    })
    rez <- do.call(rbind, rez)
    rez <- t(apply(rez, 1, sort))
    rez <- cbind(rez, rep(0, nrow(rez)))
    colnames(rez) <- c("ID1", "ID2", "background")
    object <- set(object, "Neighbors", rez)
    object <- set(object, "keysUsed", list(blockVar = blockVar, sortKeys = keys))
    return(object)
}

#' @rdname blockedSN
#' @export 
blockedSN.Scores <- function(object, blockVar, repSN, windowSN, keyLength) {
    keyVars <- object[["keyVars"]][, "keyVars"]
    keyVars <- keyVars[keyVars != blockVar]
    # generate pseudo random key lengths if none provided
    if (missing(keyLength)) {
        stats::runif(1)
        keyLength <- abs(.Random.seed[1:repSN])%%10
        idx <- keyLength > length(keyVars)
        keyLength[idx] <- length(keyVars)
    }
    if (length(keyLength != repSN)) {
        keyLength <- rep_len(keyLength, repSN)
    }
    # generate keys to sort on
    keys <- lapply(1:repSN, function(sn) {
        keyVars[sample(length(keyVars), keyLength)]
    })
    # get unique values to block on
    blockVec <- object[["rawData"]][, blockVar]
    iterval <- unique(blockVec)
    iter <- 1:length(iterval)
    # partition input data to separate blocks
    subdat <- lapply(iter, function(it) {
        object[["rawData"]][blockVec == iterval[it], ]
    })
    names(subdat) <- iter
    rez <- lapply(iter, function(it) {
        if (nrow(subdat[[it]])==1) {return(NULL)}
        if (nrow(subdat[[it]])==2) {
            nei <- cbind(subdat[[it]][1,object[["ID"]]],subdat[[it]][2,object[["ID"]]])
            colnames(nei) <- c("sorted.ids", "")
            return(nei)
        }
        nei <- lapply(1:repSN, function(sn) {
            if (nrow(subdat[[it]]) > windowSN) {
                sortedNeighbors(subdat[[it]], keys[[sn]], windowSN = windowSN, ID = object[["ID"]])
            } else {
                sortedNeighbors(subdat[[it]], keys[[sn]], windowSN = nrow(subdat[[it]])-1, ID = object[["ID"]])
            }
        })
        do.call(rbind, nei)
    })
    rez <- do.call(rbind, rez)
    rez <- t(apply(rez, 1, sort))
    rez <- cbind(rez, rep(0, nrow(rez)), rep(NA, nrow(rez)))
    colnames(rez) <- c("ID1", "ID2", "background", "matchScore")
    object <- set(object, "Neighbors", rez)
    object <- set(object, "keysUsed", list(blockVar = blockVar, sortKeys = keys))
    return(object)
}
