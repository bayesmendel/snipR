#' Scores object.
#'
#' Stores the scored neighbors and parameters from calling 
#' \code{\link{scoreNeighbors}}.
#'
#' @param BlocksObj \code{\link{Blocks}} object to score.
#' @param Scores Numeric vector of match scores.
#' @param method The type of score. Can be \code{"intersection"}, \code{"greedy"},
#' or \code{"both"}.
#' @return An object of class \code{\link{Scores}}. 
#'
#' @export
Scores <- function(BlocksObj, Scores, method) {
    ScoresObj <- set(BlocksObj, "scoreVec", Scores)
    ScoresObj$method <- method
    class(ScoresObj) <- "Scores"
    return(ScoresObj)
}
