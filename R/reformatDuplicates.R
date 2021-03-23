#' format deduplicated data back to a list for exporting.
#'
#' @param object \code{\link{Duplicates}} object. 
#' @return An object of class list containing data matrices ordered according
#' to their duplicate status (1=unique, 2=duplicates1, 3=duplicates2, etc.)
#'
#' @export
reformatDuplicates <- function(object) {
    UseMethod("reformatDuplicates", object)
}

#' @rdname reformatDuplicates
#' @export 
reformatDuplicates.default <- function(object) {
    print("reformatting allowable on Duplicates objects only.")
    return(NULL)
}

#' @rdname reformatDuplicates
#' @export 
reformatDuplicates.Duplicates <- function(object) {
    hasDup <- do.call(c, lapply(object$dupsList, function(i) {
        i[, object[["ID"]]]
    }))
    oldIDs <- unique(object[["rawData"]][, object[["ID"]]])
    noDuplicates <- oldIDs[!(oldIDs %in% hasDup)]
    uniqueData <- merge(object[["rawData"]], noDuplicates, by.x = object$ID, by.y = 1, all.y = TRUE)
    uniqueData <- cbind(uniqueData, uniqueID = uniqueData[, object$ID])
    maxCount <- 2
    for (i in 1:length(object$dupsList)) {
        maxCount <- max(maxCount, nrow(object$dupsList[[i]]))
    }
    cleanList <- lapply(1:maxCount, function(i) {
        curr <- lapply(object$dupsList, function(k) {
            if (nrow(k) >= i) {
                k[i, ]
            }
        })
        do.call(rbind, curr)
    })
    names(cleanList) <- paste0("isDup", 1:maxCount)
    cleanList <- c(list(rawData = object[["rawData"]], noDups = uniqueData), cleanList)
    class(cleanList) <- "list"
    return(cleanList)
}
