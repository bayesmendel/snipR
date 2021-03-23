#' Duplicates object.
#'
#' Stores the results from calling \code{\link{deDuplicate}}.
#'
#' @param ScoresObj \code{\link{Scores}} object.

#' @param dupsList A list of the duplicate entities. Each duplicate entity is
#' a vector of request IDs, where the corresponding families are all considered
#' duplicates of each other.
#' @param dupsReps A vector of duplicate ID representatives for each duplicate entity.
#' The chosen representative is based on the priority variable used in \code{\link{deDuplicate}}.
#' @param details Details to be printed regarding the thresholding and ordering
#' @return An object of class \code{\link{Duplicates}}. 
#'
#' @export
Duplicates <- function(ScoresObj, dupsList, dupsReps, details) {
    DupsObj <- ScoresObj
    DupsObj[["dupsList"]] <- dupsList
    DupsObj[["dupsReps"]] <- dupsReps
    # DupsObj[["datDedup"]] <- datDedup
    DupsObj[["details"]] <- details
    class(DupsObj) <- "Duplicates"
    return(DupsObj)
}
