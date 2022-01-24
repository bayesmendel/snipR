#' Calculate match score generic.
#'
#' \code{scoreNeighbors} returns a \code{\link{Scores}} object with the \code{Neighbors}
#' slot populated with a match score for each potential duplicate found from the SNIP algorithm.
#'
#' @param method Character value.  If "intersection", then uses the count of times the match was found.
#' If "greedy", then uses the greedy match score.
#' @param object object containing pre-computed Neighbors from the SNIP algorithm.
#' @return An object of class \code{\link{Scores}} containing the scored neighbor pairs.
#'
#' @export
scoreNeighbors <- function(object, method='intersection') {
    UseMethod("scoreNeighbors", object)
}

#' @rdname scoreNeighbors
#' @export
scoreNeighbors.default <- function(object, method='intersection') {
    print("Scoring allowable on Blocks and Scores objects only.")
    print("Please first block the Neighbors object.")
    return(NULL)
}

#' @param method Character value.  If intersection, then uses the normalized count of times the match was found.
#' If "greedy", then uses the greedy match score.
#' @rdname scoreNeighbors
#' @export
scoreNeighbors.Blocks <- function(object, method='intersection') {
    if (method == 'intersection') {
        object <- Scores(object, object$Neighbors[,4], method)
        return(object)
    } else {
        numericVars <- object$keyVars[object$keyVars[, "keyType"] == "numeric", "keyVars"]
        numericWts <- as.numeric(object$keyVars[object$keyVars[, "keyType"] == "numeric", "keyWt"])
        matchVars <- object$keyVars[object$keyVars[, "keyType"] != "numeric", "keyVars"]
        matchWts <- as.numeric(object$keyVars[object$keyVars[, "keyType"] != "numeric", "keyWt"])
        ID1 <- object$Neighbors[, 1]
        ID2 <- object$Neighbors[, 2]
        merged1 <- merge(object$rawData, ID1, by.x = object$ID, by.y = 1, all.y = TRUE)
        merged2 <- merge(object$rawData, ID2, by.x = object$ID, by.y = 1, all.y = TRUE)
        numericData1 <- merged1[, numericVars]
        numericData2 <- merged2[, numericVars]
        matchData1 <- merged1[, matchVars]
        matchData2 <- merged2[, matchVars]
        # match score for numeric variables
        numericScore <- data.matrix(numericData1) - data.matrix(numericData2)
        numericScore <- ifelse(numericScore == 0, 1, (1/(abs(numericScore)+1)))
        numericScore <- numericScore * matrix(rep(numericWts, each = nrow(numericScore)), nrow = nrow(numericScore))
        scoreVec <- rowSums(numericScore, na.rm = TRUE)
        # match score for match variables
        matchScore <- data.matrix(matchData1 == matchData2) * 1
        matchScore <- matchScore * matrix(rep(matchWts, each = nrow(matchScore)), nrow = nrow(matchScore))
        scoreVec <- scoreVec + rowSums(matchScore, na.rm = TRUE)
        # take the mean
        scoreVec <- scoreVec/nrow(object$keyVars)
        # promote to Scores object
        object <- Scores(object, scoreVec, method)
        return(object)
    }
}

#' @param method Character value.  If intersection, then uses the count of times the match was found.
#' If "greedy", then uses the greedy match score.
#' @rdname scoreNeighbors
#' @export
scoreNeighbors.Scores <- function(object, method='intersection') {
    if (method == 'intersection') {
        object <- set(object, "scoreVec", object$Neighbors[,4])
        return(object)
    } else {
        numericVars <- object$keyVars[object$keyVars[, "keyType"] == "numeric", "keyVars"]
        numericWts <- as.numeric(object$keyVars[object$keyVars[, "keyType"] == "numeric", "keyWt"])
        matchVars <- object$keyVars[object$keyVars[, "keyType"] != "numeric", "keyVars"]
        matchWts <- as.numeric(object$keyVars[object$keyVars[, "keyType"] != "numeric", "keyWt"])
        idx <- is.na(object$Neighbors[, "matchScore"])
        ID1 <- object$Neighbors[idx, 1]
        ID2 <- object$Neighbors[idx, 2]
        merged1 <- merge(object$rawData, ID1, by.x = object$ID, by.y = 1, all.y = TRUE)
        merged2 <- merge(object$rawData, ID2, by.x = object$ID, by.y = 1, all.y = TRUE)
        numericData1 <- merged1[, numericVars]
        numericData2 <- merged2[, numericVars]
        matchData1 <- merged1[, matchVars]
        matchData2 <- merged2[, matchVars]
        # match score for numeric variables
        numericScore <- data.matrix(numericData1) - data.matrix(numericData2)
        numericScore <- ifelse(numericScore == 0, 1, (1/(abs(numericScore)+1)))
        numericScore <- numericScore * matrix(rep(numericWts, each = nrow(numericScore)), nrow = nrow(numericScore))
        scoreVec <- rowSums(numericScore, na.rm = TRUE)
        # match score for match variables
        matchScore <- data.matrix(matchData1 == matchData2) * 1
        matchScore <- matchScore * matrix(rep(matchWts, each = nrow(matchScore)), nrow = nrow(matchScore))
        scoreVec <- scoreVec + rowSums(matchScore, na.rm = TRUE)
        # take the mean
        scoreVec <- scoreVec/nrow(object$keyVars)
        new <- rep_len(NA, nrow(object$Neighbors))
        new[idx] <- scoreVec
        new[!idx] <- object$Neighbors[, "matchScore"]
        scoreVec <- new
        object <- set(object, "scoreVec", scoreVec)
        return(object)
    }
}
