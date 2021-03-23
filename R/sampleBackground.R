#' Populate Neighbors with random pairs.
#'
#' @param object object containing pre-computed Neighbors from the BSN algorithm.
#' @param num Number or proportion of current Neighbors if <1 to populate with.
#' @return An object of class \code{\link{Neighbors}} containing the possible duplicate pairs
#' sampled at random. 
#'
#' @export
sampleBackground <- function(object, num) {
    UseMethod("sampleBackground", object)
}

#' @rdname sampleBackground
#' @export 
sampleBackground.default <- function(object, num) {
    print("Sampling allowable on Neighbors, Blocks, and Scores objects only.")
    return(NULL)
}

#' @rdname sampleBackground
#' @export 
sampleBackground.Neighbors <- function(object, num) {
    if (num < 1) {
        print("Provide a number >1. Defaulting to 100.")
        num <- 100
    }
    a1 <- sample(object[["rawData"]][, object[["ID"]]], size = num, replace = TRUE)
    a2 <- sample(object[["rawData"]][, object[["ID"]]], size = num, replace = TRUE)
    rez <- cbind(a1, a2)
    rez <- t(apply(rez, 1, sort))
    rez <- unique(rez)
    rez <- cbind(rez, rep(1, nrow(rez)))
    colnames(rez) <- c("ID1", "ID2", "background")
    # promote to Blocks object
    object <- Blocks(object, rez, list(blockVar = "Background Sampling", sortKeys = paste0("N=", num)))
    return(object)
}

#' @rdname sampleBackground
#' @export 
sampleBackground.Blocks <- function(object, num) {
    if (num <= 1) {
        num <- round(abs(num) * length(object[["Neighbors"]][, 3] == 0))
    }
    if (num <= 1) {
        num <- round(abs(num) * nrow(object[["Neighbors"]]))
    }
    a1 <- sample(object[["rawData"]][, object[["ID"]]], size = num, replace = TRUE)
    a2 <- sample(object[["rawData"]][, object[["ID"]]], size = num, replace = TRUE)
    rez <- cbind(a1, a2)
    rez <- t(apply(rez, 1, sort))
    rez <- unique(rez)
    rez <- rez[rez[, 1] != rez[, 2], ]
    rez <- cbind(rez, rep(1, nrow(rez)))
    colnames(rez) <- c("ID1", "ID2", "background")
    object <- set(object, "Neighbors", rez)
    return(object)
}

#' @rdname sampleBackground
#' @export 
sampleBackground.Scores <- function(object, num) {
    if (num < 1) {
        num <- round(abs(num) * length(object[["Neighbors"]][, 3] == 0))
    }
    if (num < 1) {
        num <- round(abs(num) * nrow(object[["Neighbors"]]))
    }
    a1 <- sample(object[["rawData"]][, object[["ID"]]], size = num, replace = TRUE)
    a2 <- sample(object[["rawData"]][, object[["ID"]]], size = num, replace = TRUE)
    rez <- cbind(a1, a2)
    rez <- t(apply(rez, 1, sort))
    rez <- unique(rez)
    rez <- rez[rez[, 1] != rez[, 2], ]
    rez <- cbind(rez, rep(1, nrow(rez)))
    colnames(rez) <- c("ID1", "ID2", "background")
    # create temporary objects to score new background samplings
    tempobject <- Neighbors(object[["rawData"]], object[["ID"]], object[["keyVars"]])
    tempobject <- Blocks(tempobject, rez, list(blockVar = "Background Sampling", sortKeys = paste0("N=", num)))
    tempobject <- scoreNeighbors(tempobject)
    newScores <- tempobject[["Neighbors"]]
    # and consolidate it
    object <- set(object, "Neighbors", newScores)
    return(object)
}
