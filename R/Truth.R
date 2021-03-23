#' Truth object.
#'
#' Stores the results from calling \code{\link{revealTruth}}.
#'
#' @param dupObj \code{\link{Duplicates}} object.
#' @param neighborsTruth A data frame providing the match scores and duplication truth
#' outcome (true positive, true negative, false positive, false negative) for each
#' pair of neighbors and threshold.
#' @param famsTruth A data frame providing the representation status for each family
#' (whether they are still included in the deduplicated data) as well as the number of
#' extra representatives of the family.
#'
#' @export
Truth <- function(dupObj, neighborsTruth, famsTruth) {
  TruthObj <- dupObj
  TruthObj[["neighborsTruth"]] <- neighborsTruth
  TruthObj[["famsTruth"]] <- famsTruth
  class(TruthObj) <- "Truth"
  return(TruthObj)
}
