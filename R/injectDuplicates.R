#' Inject duplicates.
#'
#' \code{injectDuplicates} returns a data frame with known duplicates added to the rows.
#'
#' @param orig.data Data
#' @param pDup Numeric value. Proportion of dataset number of rows to inject duplicates with.
#' @param errRt Parameter of the Poisson distribution for the number of errors in each duplicate.
#' @return An object of class \code{\link{Blocks}} containing the neighbors found
#' and keys used during the blocked sorted neighbors iteration. 
#'
#' @export
injectDuplicates.Neighbors <- function(orig.data, pDup, errRt) {
  nInjections <- pDup*nrow(object[["rawData"]])
  browser()
  for (i in 1:nInjections) {
    if (fams[[id]][1,"Duplicate"]>0) {
      fams[[id]]$Duptype <- "original"
      for (j in 1:fams[[id]][1,"nDuplicates"]) {
        tfam2 <- simulateErrors(trimfam, errRt=errRt)
        errs <- tfam2[[2]]
        tfam2 <- tfam2[[1]]
        tid <- paste(sample(seq(0,9,1),10,replace=TRUE),collapse="")
        tfam2$requestId <- tid
        tfam2$Duptype <- "isDup"
        tfam2$nErrors <- errs
        fams[[tid]] <- tfam2
      }
    } else {fams[[id]]$Duptype <- "noDup"}
  }
  return(object)
}
