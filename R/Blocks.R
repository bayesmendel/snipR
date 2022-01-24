#' Blocks object.
#'
#' Stores the neighbors and parameters from calling \code{blockedSN}.
#'
#' @param NeighborsObj \code{\link{Neighbors}} object to perform SNIP.
#' @param Neighbors Data frame with rows as potential matches and
#' columns for potential match pair IDs. Sorted.
#' @param keysUsed List containing blocking variables and sort keys
#' for each iteration of the algorithm.
#' @return An object of class \code{\link{Blocks}}. 
#'
#' @export
Blocks <- function(NeighborsObj, Neighbors, keysUsed) {
    BlocksObj <- NeighborsObj
    BlocksObj[["Neighbors"]] <- Neighbors
    BlocksObj[["keysUsed"]] <- c(iter1 = list(keysUsed))
    class(BlocksObj) <- "Blocks"
    return(BlocksObj)
}
