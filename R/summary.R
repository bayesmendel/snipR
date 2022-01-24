#' Summarize object generic.
#'
#' \code{summary} summarizes the contents of the object.
#'
#' @param object object to be summarized.
#' @param ... other arguments to be passed.
#' @rdname summary
#' @export 
summary.Neighbors <- function(object, ...) {
    ans <- list()
    ans[["Class"]] <- "Object of class Neighbors"
    ans[["Current number of unique observations: "]] <- nrow(object$rawData)
    ans[["Possible pairwise comparisons: "]] <- nrow(object$rawData) * (nrow(object$rawData) - 1)/2
    ans[["Pairwise comparisons made: "]] <- 0
    ans[["Identifying variable: "]] <- object$ID
    ans[["Binary keys: "]] <- object$keyVars[object$keyVars[, 2] == "binary", ]
    ans[["String keys: "]] <- object$keyVars[object$keyVars[, 2] == "string", ]
    ans[["Numeric keys: "]] <- object$keyVars[object$keyVars[, 2] == "numeric", ]
    ans
}

#' @rdname summary
#' @export 
summary.Blocks <- function(object, ...) {
    ans <- list()
    ans[["Class"]] <- "Object of class Blocks"
    ans[["Current number of unique observations: "]] <- nrow(object$rawData)
    ans[["Possible pairwise comparisons: "]] <- nrow(object$rawData) * (nrow(object$rawData) - 1)/2
    ans[["Pairwise comparisons made: "]] <- nrow(object$Neighbors)
    ans[["Percent comparisons considered: "]] <- nrow(object$Neighbors)/(nrow(object$rawData) * (nrow(object$rawData) - 1)/2) * 100
    ans
}

#' @rdname summary
#' @export 
summary.Scores <- function(object, ...) {
    ans <- list()
    ans[["Class"]] <- "Object of class Scores"
    ans[["Current number of unique observations: "]] <- nrow(object$rawData)
    ans[["Possible pairwise comparisons: "]] <- nrow(object$rawData) * (nrow(object$rawData) - 1)/2
    ans[["Pairwise comparisons made: "]] <- nrow(object$Neighbors)
    ans[["Percent comparisons considered: "]] <- nrow(object$Neighbors)/(nrow(object$rawData) * (nrow(object$rawData) - 1)/2) * 100
    ans[["Match score distribution"]] <- summary(as.numeric(object$Neighbors[, "matchScore"]))
    ans
}

#' @rdname summary
#' @export 
summary.Duplicates <- function(object, ...) {
    ans <- list()
    ans[["Class"]] <- "Object of class Duplicates"
    ans[["Duplicates ordered on: "]]
    if (object$details[[4]]) {
        ans[["Duplicates thresholded at: "]] <- paste0(object$details[[1]], " percentile (", object$details[[2]], ") in descending order on ", 
            object$details[[3]])
    } else {
        ans[["Duplicates thresholded at: "]] <- paste0(object$details[[1]], " percentile (", object$details[[2]], ") in ascending order on ", 
            object$details[[3]])
        
    }
    ans[["Current number of duplicate entities: "]] <- length(object$dupsList)
    maxSize <- 2
    for (i in 1:length(object$dupsList)) {
        maxSize <- max(maxSize, nrow(object$dupsList[[i]]))
    }
    ans[["Largest entity size: "]] <- maxSize
    ans
}

