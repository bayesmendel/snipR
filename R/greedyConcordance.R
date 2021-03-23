#' Calculate greedy match score.
#'
#' \code{greedyConcordance} scores pedigree members in a sequential greedy fashion
#' maximizing the match score and eliminating that family member from future
#' consideration at each iteration.
#'
#' @param fam1 Data frame containing family members 1
#' @param fam2 Data frame containing familiy members 2
#' @param keyVars key variables matarix
#' @return Score value, numeric.
#'
greedyConcordance <- function(fam1, fam2, keyVars) {
    numericVars <- keyVars[keyVars[, "keyType"] == "numeric", "keyVars"]
    numericWts <- as.numeric(keyVars[keyVars[, "keyType"] == "numeric", "keyWt"])
    matchVars <- keyVars[keyVars[, "keyType"] != "numeric", "keyVars"]
    matchWts <- as.numeric(keyVars[keyVars[, "keyType"] != "numeric", "keyWt"])
    if (nrow(fam1)>nrow(fam2)) {
        fam0 <- fam1
        fam1 <- fam2
        fam2 <- fam1
    }
    numericData1 <- fam1[, numericVars]
    numericData2 <- fam2[, numericVars]
    matchData1 <- fam1[, matchVars]
    matchData2 <- fam2[, matchVars]
    score <- 0
    for (i in 1:nrow(fam1)) {
        s1match <- data.frame(matrix(rep(matchData1[i, ], each = nrow(matchData2)), nrow = nrow(matchData2)))
        s2match <- data.frame(as.matrix(matchData2, nrow = nrow(matchData2)))
        s1numeric <- data.frame(matrix(rep(numericData1[i, ], each = nrow(numericData2)), nrow = nrow(numericData2)))
        s2numeric <- data.frame(as.matrix(numericData2, nrow = nrow(numericData2)))
        idx <- !is.na(s1match) & !is.na(s2match)  # score match variables
        if (any(rowSums(idx) > 0)) {
            diff <- data.matrix(s1match[idx]) - data.matrix(s2match[idx])
            matchScr <- ifelse(diff == 0, 1, 0)
            matchScr <- matchScr * matrix(rep(matchWts, each = nrow(matchScr)), nrow = nrow(matchScr))
        }
        idy <- !is.na(s1numeric) & !is.na(s2numeric)  # score numeric variables
        if (any(rowSums(idy) > 0)) {
            diff <- data.matrix(s1numeric[idy]) - data.matrix(s2numeric[idy])
            numericScr <- ifelse(diff == 0, 1, (1/(abs(diff)+1)))
            numericScr <- numericScr * matrix(rep(numericWts, each = nrow(numericScr)), nrow = nrow(numericScr))
        }
        if (nrow(s2match) == 1) {
            score <- (sum(matchScr, na.rm = TRUE) + sum(numericScr, na.rm = TRUE))/length(c(numericVars, matchVars))
            break
        } else {
            tempScr <- rowSums(matchScr, na.rm = TRUE) + rowSums(numericScr, na.rm = TRUE)
            id <- which.max(tempScr)
            score <- score + max(tempScr)/length(c(numericVars, matchVars))
            matchData2 <- matchData2[-id, ]
        }
    }
    return(score)
}
