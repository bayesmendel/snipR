#' Plot Neighbors object generic.
#'
#' \code{plot} plots the contents of the object.
#'
#' @param x code{\link{Neighbors}} or \code{\link{Neighbors}} oject to be summarized.
#' @param type Char, How to display the plot (eg. 'box', 'density')
#' @param ... other arguments to be passed.
#' @rdname plot
#' @export
plot.Scores <- function(x, type, ...) {
    scores <- x[["Neighbors"]][, "matchScore"]
    bg <- x[["Neighbors"]][, "background"]
    if (type == "density") {
        if (length(unique(bg)) > 1) {
            sm::sm.density.compare(as.numeric(scores), as.factor(bg), xlab = "Match Score", ...)
            graphics::abline(v = x[["details"]][[1]], col = "red")
            graphics::title(main = "Match Score Density by Discovery Method")
            graphics::legend("right", c("BSN", "Random"), fill = 2:3)
        } else {
            sm::sm.density(as.numeric(scores), ...)
            graphics::abline(v = x[["details"]][[1]], col = "red")
            graphics::title(main = "Match Score Density from BSN")
        }
    } else if (type == "box") {
        if (length(unique(bg)) > 1) {
            graphics::boxplot(as.numeric(scores) ~ as.factor(bg), xlab = "Discovery Method", ylab = "Match Score", main = "Match Score by Discovery Method",
                ...)
            graphics::abline(h = x[["details"]][[1]], col = "red")
        } else {
            graphics::boxplot(as.numeric(scores), ylab = "Match Score", main = "Match Score from BSN", ...)
            graphics::abline(h = x[["details"]][[1]], col = "red")
        }
    } else {
        print("Please provide a valid type argmument.")
    }

}

#' @rdname plot
#' @export
plot.Duplicates <- function(x, type, ...) {
    scores <- x[["Neighbors"]][, "matchScore"]
    thresh <- x[["details"]][[2]]
    bg <- as.numeric(scores > thresh)
    if (type == "density") {
        sm::sm.density.compare(as.numeric(scores), as.factor(bg), xlab = "Match Score", ...)
        graphics::abline(v = x[["details"]][[2]], col = "red")
        graphics::title(main = "Match Score Density by Duplicate Pair Status")
        graphics::legend("right", c("Not Dup", "Dup"), fill = 2:3)
    } else if (type == "box") {
        graphics::boxplot(as.numeric(scores) ~ as.factor(bg), xlab = "Duplicate Pair", ylab = "Match Score", main = "Match Score by Duplicate Pair Status",
            ...)
        graphics::abline(h = x[["details"]][[2]], col = "red")
    } else if (type == "family counts") {
      sizes <- matrix(data=NA,nrow=length(x$dupsList),ncol=1)
      for (i in 1:length(x$dupsList)) {
        sizes[i,1] <- nrow(x$dupsList[[i]])
      }
      graphics::hist(sizes, breaks=2:max(sizes), xlab = "Number of Occurrences", main = "Family Duplicate Counts",
                        ...)
    } else {
        print("Please provide a valid type argmument.")
    }
}

