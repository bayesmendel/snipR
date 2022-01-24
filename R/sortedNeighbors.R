#' Perform a single iteration of sorted neighbors.
#'
#' \code{sortedNeighbors} sorts the input data frame on the specified sort key
#' and returns all neighbors within a specified window.
#'
#' @param dat Data frame. Data used for sorting and discovering neighbors.
#' @param key Integer or character vector. The columns on which to construct the sort key.
#' @param windowSN Integer value. Size of sliding window to use during sorted neighbors.
#' @param ID Integer or character value. Identifiable column in dat.
#' @return A data frame with 2 columns and n rows, each row containing neighbor pairs.
#'
#' @export
sortedNeighbors <- function(dat, key, windowSN, ID) {
    # save.image('temp.Rda')
    str <- paste0("dat[order(", paste(paste0("dat$", key), collapse = ","), "),]")
    tempdat <- eval(parse(text = str))
    sorted.ids <- tempdat[, ID]
    str <- "cbind(sorted.ids"
    if (windowSN >= 1) {
        for (i in 1:(windowSN)) {
            str <- paste0(str, ",c(tail(sorted.ids, -", i, "), rep(NA,", i, "))")
        }
    }
    str <- paste0(str, ")")
    famsMat <- eval(parse(text = str))
    str <- "rbind(famsMat[,c(1,2)]"
    if (windowSN > 1) {
        for (i in 2:windowSN) {
            str <- paste0(str, ",famsMat[,c(1,", i + 1, ")]")
        }
    }
    str <- paste0(str, ")")
    famsout <- eval(parse(text = str))
    famsout <- famsout[stats::complete.cases(famsout), ]
    return(famsout)
}
