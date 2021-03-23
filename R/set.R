#' Set Neighbors attribute generic.
#'
#' \code{set} sets the attribute for the \code{\link{Neighbors}} object.
#'
#' @param object object to be updated.
#' @param attr Character input. The attribute to update (e.g. 'rawData', 'ID') 
#' @param val The new values to add or update. For attributes 'Neighbors' and 'keysUsed',
#' \code{set} will concatenate the values. Otherwise old values will be replaced.
#' @return object of same class (\code{\link{Neighbors}}, \code{\link{Blocks}}, etc.).
#'
#' @export
set <- function(object, attr, val) {
    UseMethod("set", object)
}

#' @rdname set
#' @export 
set.Neighbors <- function(object, attr, val) {
    if (attr %in% c("keyType", "keyWt", "keyVars")) {
        if (ncol(as.matrix(val)) > 1) {
            object[["keyVars"]] <- val
        } else {
            object[["keyVars"]][, attr] <- val
        }
    } else if (attr %in% c("ID", "rawData")) {
        object[[attr]] <- val
    } else {
        print("Unrecognized attribute")
        return(NULL)
    }
    return(object)
}

#' @rdname set
set.Blocks <- function(object, attr, val) {
    if (attr == "keysUsed") {
        oldKeys <- object[["keysUsed"]]
        numPrevious <- length(oldKeys)
        str <- paste0("c(oldKeys, iter", numPrevious + 1, "=list(val))")
        newKeys <- eval(parse(text = str))
        object[["keysUsed"]] <- newKeys
    } else if (attr == "Neighbors") {
        oldNeighbors <- object[["Neighbors"]]
        if (is.null(oldNeighbors)) {
            print("You may not initialize this manually. Call a new Blocks object.")
            return(NULL)
        } else {
            rez <- rbind(oldNeighbors, val)
            rez <- rez[order(rez[, "background"]), ]
            rez <- rez[!duplicated(rez[, 1:2]), ]
            object[["Neighbors"]] <- rez
        }
    } else if (attr == "scoreVec") {
        oldNeighbors <- object[["Neighbors"]]
        rez <- cbind(oldNeighbors, matchScore = val)
        object[["Neighbors"]] <- rez
    } else if (attr %in% c("keyType", "keyWt", "keyVars")) {
        if (ncol(as.matrix(val)) > 1) {
            object[["keyVars"]] <- val
        } else {
            object[["keyVars"]][, attr] <- val
        }
    } else if (attr %in% c("ID", "rawData")) {
        print("Initialize a new Neighbors object to modify this attribute.")
        return(object)
    } else {
        print("Unrecognized attribute")
        return(NULL)
    }
    return(object)
}

#' @rdname set
set.Scores <- function(object, attr, val) {
    if (attr == "keysUsed") {
        oldKeys <- object[["keysUsed"]]
        numPrevious <- length(oldKeys)
        str <- paste0("c(oldKeys, iter", numPrevious + 1, "=list(val))")
        newKeys <- eval(parse(text = str))
        object[["keysUsed"]] <- newKeys
    } else if (attr == "Neighbors") {
        oldNeighbors <- object[["Neighbors"]]
        if (is.null(oldNeighbors)) {
            print("You may not initialize this manually. Call a new Blocks object.")
            return(NULL)
        } else {
            rez <- rbind(oldNeighbors, val)
            rez <- rez[order(rez$background), ]
            rez <- rez[!duplicated(rez[, 1:2]), ]
            object[["Neighbors"]] <- rez
        }
    } else if (attr == "scoreVec") {
        oldNeighbors <- object[["Neighbors"]]
        oldNeighbors[, "matchScore"] <- val
        object[["Neighbors"]] <- oldNeighbors
    } else if (attr %in% c("keyType", "keyWt", "keyVars")) {
        if (ncol(as.matrix(val)) > 1) {
            object[["keyVars"]] <- val
        } else {
            object[["keyVars"]][, attr] <- val
        }
    } else if (attr %in% c("ID", "rawData")) {
        print("Initialize a new Neighbors object to modify this attribute.")
        return(object)
    } else {
        print("Unrecognized attribute")
        return(NULL)
    }
    return(object)
}
