#' Calculate match score generic.
#'
#' \code{scoreFamilies} returns a \code{\link{Scores}} object with the \code{Neighbors}
#' slot populated with a match score for each potential duplicate found from the BSN algorithm.
#' A match score is calculated by averaging the number of string/binary variable matches
#' consistent between the two records plus a similarity measure for the continuous variables.
#'
#' @param object object containing pre-computed Neighbors from the BSN algorithm.
#' @param familyDat data matrix containing family data so score.
#' @param selfID character referencing the column in familyDat which points to a self's ID.
#' @param motherID character referencing the column in familyDat which points to a self's mother's ID
#' @param fatherID character referencing the column in familyDat which points to a self's father's ID
#' @return An object of class \code{\link{Scores}} containing the scored neighbor pairs.
#'
#' @export
scoreFamilies <- function(object, familyDat, selfID, motherID, fatherID) {
    UseMethod("scoreFamilies", object)
}

#' @rdname scoreFamilies
#' @export
scoreFamilies.default <- function(object, familyDat, selfID, motherID, fatherID) {
    print("Scoring allowable on Blocks and Scores objects only.")
    print("Please first block the Neighbors object.")
    return(NULL)
}

#' @rdname scoreFamilies
#' @export
scoreFamilies.Blocks <- function(object, familyDat, selfID, motherID, fatherID) {
  # score selfs
    selfScore <- scoreNeighbors(object)$Neighbors[, "matchScore"]
    # score mothers
    mothersObj <- object
    newDat <- suppressWarnings(merge(object$rawData[, c(object[["ID"]], motherID)], familyDat, by.x = c(object[["ID"]], motherID),
        by.y = c(object[["ID"]], selfID), all.x = TRUE, suffixes = c("old", "new"))[, -2])
    mothersObj[["rawData"]] <- newDat
    motherScore <- scoreNeighbors(mothersObj)$Neighbors[, "matchScore"]
    # score fathers
    fathersObj <- object
    newDat <- suppressWarnings(merge(object$rawData[, c(object[["ID"]], fatherID)], familyDat, by.x = c(object[["ID"]], fatherID),
        by.y = c(object[["ID"]], selfID), all.x = TRUE, suffixes = c("old", "new"))[, -2])
    fathersObj[["rawData"]] <- newDat
    fatherScore <- scoreNeighbors(fathersObj)$Neighbors[, "matchScore"]
    # score fathers fathers
    ffathersObj <- fathersObj
    newDat <- suppressWarnings(merge(fathersObj$rawData[, c(object[["ID"]], fatherID)], familyDat, by.x = c(object[["ID"]], fatherID),
        by.y = c(object[["ID"]], selfID), all.x = TRUE, suffixes = c("old", "new"))[, -2])
    ffathersObj[["rawData"]] <- newDat
    ffatherScore <- scoreNeighbors(ffathersObj)$Neighbors[, "matchScore"]
    # score fathers mothers
    fmothersObj <- fathersObj
    newDat <- suppressWarnings(merge(fathersObj$rawData[, c(object[["ID"]], motherID)], familyDat, by.x = c(object[["ID"]], motherID),
        by.y = c(object[["ID"]], selfID), all.x = TRUE, suffixes = c("old", "new"))[, -2])
    fmothersObj[["rawData"]] <- newDat
    fmotherScore <- scoreNeighbors(fmothersObj)$Neighbors[, "matchScore"]
    # score mothers fathers
    mfathersObj <- mothersObj
    newDat <- suppressWarnings(merge(mothersObj$rawData[, c(object[["ID"]], fatherID)], familyDat, by.x = c(object[["ID"]], fatherID),
        by.y = c(object[["ID"]], selfID), all.x = TRUE, suffixes = c("old", "new"))[, -2])
    mfathersObj[["rawData"]] <- newDat
    mfatherScore <- scoreNeighbors(mfathersObj)$Neighbors[, "matchScore"]
    # score mothers mothers
    mmothersObj <- mothersObj
    newDat <- suppressWarnings(merge(mothersObj$rawData[, c(object[["ID"]], motherID)], familyDat, by.x = c(object[["ID"]], motherID),
        by.y = c(object[["ID"]], selfID), all.x = TRUE, suffixes = c("old", "new"))[, -2])
    mmothersObj[["rawData"]] <- newDat
    mmotherScore <- scoreNeighbors(mmothersObj)$Neighbors[, "matchScore"]
    # score everything else this is time consuming - any better ideas?
    allScore <- sapply(1:nrow(object$Neighbors), function(id) {
        tsc <- try({
            fam1 <- familyDat[familyDat[,object$ID] == object$Neighbors[id, 1], ]
            fam2 <- familyDat[familyDat[,object$ID] == object$Neighbors[id, 2], ]
            self1 <- fam1[fam1[, selfID] == 1, ]  # find self
            fam1 <- fam1[fam1[, selfID] != 1, ]  # remove self
            mother1 <- fam1[fam1[, selfID] == self1[, motherID], ]  # find mother
            fam1 <- fam1[fam1[, selfID] != self1[, motherID], ]  # remove mother
            fam1 <- fam1[fam1[, selfID] != mother1[, motherID], ]  # remove mothers mother
            fam1 <- fam1[fam1[, selfID] != mother1[, fatherID], ]  # remove mothers father
            father1 <- fam1[fam1[, selfID] == self1[, fatherID], ]  # find father
            fam1 <- fam1[fam1[, selfID] != self1[, fatherID], ]  # remove father
            fam1 <- fam1[fam1[, selfID] != father1[, motherID], ]  # remove fathers mother
            fam1 <- fam1[fam1[, selfID] != father1[, fatherID], ]  # remove fathers father
            sibs1 <- unique(rbind(fam1[fam1[, motherID] == mother1[, selfID], ],
                                  fam1[fam1[, fatherID] == father1[, selfID], ])) # find siblings
            fam1 <- fam1[!fam1[,selfID] %in% sibs1[,selfID]] # remove siblings
            aunts1 <- unique(rbind(fam1[fam1[, motherID] == mother1[, motherID], ],
                                  fam1[fam1[, motherID] == father1[, motherID], ],
                                  fam1[fam1[, fatherID] == mother1[, fatherID], ],
                                  fam1[fam1[, fatherID] == father1[, fatherID], ])) # find aunts/uncles
            fam1 <- fam1[!fam1[,selfID] %in% aunts1[,selfID]] # remove siblings
            kids1 <- unique(rbind(fam1[fam1[, motherID] == self1[, selfID], ],
                                  fam1[fam1[, fatherID] == self1[, selfID], ])) # find children
            fam1 <- fam1[!fam1[,selfID] %in% kids1[,selfID]] # remove children
            self2 <- fam2[fam2[, selfID] == 1, ]  # repeat for comparison family
            fam2 <- fam2[fam2[, selfID] != 1, ]  # remove self
            mother2 <- fam2[fam2[, selfID] == self2[, motherID], ]  # find mother
            fam2 <- fam2[fam2[, selfID] != self2[, motherID], ]  # remove mother
            fam2 <- fam2[fam2[, selfID] != mother2[, motherID], ]  # remove mothers mother
            fam2 <- fam2[fam2[, selfID] != mother2[, fatherID], ]  # remove mothers father
            father2 <- fam2[fam2[, selfID] == self2[, fatherID], ]  # find father
            fam2 <- fam2[fam2[, selfID] != self2[, fatherID], ]  # remove father
            fam2 <- fam2[fam2[, selfID] != father2[, motherID], ]  # remove fathers mother
            fam2 <- fam2[fam2[, selfID] != father2[, fatherID], ]  # remove fathers father
            sibs2 <- unique(rbind(fam2[fam2[, motherID] == mother2[, selfID], ],
                                  fam2[fam2[, fatherID] == father2[, selfID], ])) # find siblings
            fam2 <- fam2[!fam2[,selfID] %in% sibs2[,selfID]] # remove siblings
            aunts2 <- unique(rbind(fam2[fam2[, motherID] == mother2[, motherID], ],
                                   fam2[fam2[, motherID] == father2[, motherID], ],
                                   fam2[fam2[, fatherID] == mother2[, fatherID], ],
                                   fam2[fam2[, fatherID] == father2[, fatherID], ])) # find aunts/uncles
            fam2 <- fam2[!fam2[,selfID] %in% aunts2[,selfID]] # remove siblings
            kids2 <- unique(rbind(fam2[fam2[, motherID] == self2[, selfID], ],
                                  fam2[fam2[, fatherID] == self2[, selfID], ])) # find children
            fam2 <- fam2[!fam2[,selfID] %in% kids2[,selfID]] # remove children

            sibscore <- greedyConcordance(sibs1, sibs2, object$keyVars)  # try to concord siblings
            auntscore <- greedyConcordance(aunts1, aunts2, object$keyVars)  # try to concord aunts/uncles
            kidscore <- greedyConcordance(kids1, kids2, object$keyVars)  # try to concord children
            otherscore <- greedyConcordance(fam1, fam2, object$keyVars)  # try to concord all others
            as.numeric(sibscore) + as.numeric(auntscore) + as.numeric(kidscore) + as.numeric(otherscore)
        }, silent = TRUE)
        ifelse(class(tsc) == "try-error", 0, tsc)
    })
    scoreVec <- as.numeric(selfScore) + as.numeric(motherScore) +
                as.numeric(fatherScore) + as.numeric(ffatherScore) +
                as.numeric(fmotherScore) + as.numeric(mfatherScore) +
                as.numeric(mmotherScore) + as.numeric(allScore)
    # promote to Scores object
    object <- Scores(object, scoreVec)
    return(object)
}

#' @rdname scoreFamilies
#' @export
scoreFamilies.Scores <- function(object, familyDat, selfID, motherID, fatherID) {
  # score selfs
  selfScore <- scoreNeighbors(object)$Neighbors[, "matchScore"]
  # score mothers
  mothersObj <- object
  newDat <- suppressWarnings(merge(object$rawData[, c(object[["ID"]], motherID)], familyDat, by.x = c(object[["ID"]], motherID),
                                   by.y = c(object[["ID"]], selfID), all.x = TRUE, suffixes = c("old", "new"))[, -2])
  mothersObj[["rawData"]] <- newDat
  motherScore <- scoreNeighbors(mothersObj)$Neighbors[, "matchScore"]
  # score fathers
  fathersObj <- object
  newDat <- suppressWarnings(merge(object$rawData[, c(object[["ID"]], fatherID)], familyDat, by.x = c(object[["ID"]], fatherID),
                                   by.y = c(object[["ID"]], selfID), all.x = TRUE, suffixes = c("old", "new"))[, -2])
  fathersObj[["rawData"]] <- newDat
  fatherScore <- scoreNeighbors(fathersObj)$Neighbors[, "matchScore"]
  # score fathers fathers
  ffathersObj <- fathersObj
  newDat <- suppressWarnings(merge(fathersObj$rawData[, c(object[["ID"]], fatherID)], familyDat, by.x = c(object[["ID"]], fatherID),
                                   by.y = c(object[["ID"]], selfID), all.x = TRUE, suffixes = c("old", "new"))[, -2])
  ffathersObj[["rawData"]] <- newDat
  ffatherScore <- scoreNeighbors(ffathersObj)$Neighbors[, "matchScore"]
  # score fathers mothers
  fmothersObj <- fathersObj
  newDat <- suppressWarnings(merge(fathersObj$rawData[, c(object[["ID"]], motherID)], familyDat, by.x = c(object[["ID"]], motherID),
                                   by.y = c(object[["ID"]], selfID), all.x = TRUE, suffixes = c("old", "new"))[, -2])
  fmothersObj[["rawData"]] <- newDat
  fmotherScore <- scoreNeighbors(fmothersObj)$Neighbors[, "matchScore"]
  # score mothers fathers
  mfathersObj <- mothersObj
  newDat <- suppressWarnings(merge(mothersObj$rawData[, c(object[["ID"]], fatherID)], familyDat, by.x = c(object[["ID"]], fatherID),
                                   by.y = c(object[["ID"]], selfID), all.x = TRUE, suffixes = c("old", "new"))[, -2])
  mfathersObj[["rawData"]] <- newDat
  mfatherScore <- scoreNeighbors(mfathersObj)$Neighbors[, "matchScore"]
  # score mothers mothers
  mmothersObj <- mothersObj
  newDat <- suppressWarnings(merge(mothersObj$rawData[, c(object[["ID"]], motherID)], familyDat, by.x = c(object[["ID"]], motherID),
                                   by.y = c(object[["ID"]], selfID), all.x = TRUE, suffixes = c("old", "new"))[, -2])
  mmothersObj[["rawData"]] <- newDat
  mmotherScore <- scoreNeighbors(mmothersObj)$Neighbors[, "matchScore"]
  # score everything else this is time consuming - any better ideas?
  allScore <- sapply(1:nrow(object$Neighbors), function(id) {
    tsc <- try({
      if (is.na(Neighbors[id, "matchScore"])) {
        fam1 <- familyDat[familyDat[,object$ID] == object$Neighbors[id, 1], ]
        fam2 <- familyDat[familyDat[,object$ID] == object$Neighbors[id, 2], ]
        self1 <- fam1[fam1[, selfID] == 1, ]  # find self
        fam1 <- fam1[fam1[, selfID] != 1, ]  # remove self
        mother1 <- fam1[fam1[, selfID] == self1[, motherID], ]  # find mother
        fam1 <- fam1[fam1[, selfID] != self1[, motherID], ]  # remove mother
        fam1 <- fam1[fam1[, selfID] != mother1[, motherID], ]  # remove mothers mother
        fam1 <- fam1[fam1[, selfID] != mother1[, fatherID], ]  # remove mothers father
        father1 <- fam1[fam1[, selfID] == self1[, fatherID], ]  # find father
        fam1 <- fam1[fam1[, selfID] != self1[, fatherID], ]  # remove father
        fam1 <- fam1[fam1[, selfID] != father1[, motherID], ]  # remove fathers mother
        fam1 <- fam1[fam1[, selfID] != father1[, fatherID], ]  # remove fathers father
        sibs1 <- unique(rbind(fam1[fam1[, motherID] == mother1[, selfID], ],
                              fam1[fam1[, fatherID] == father1[, selfID], ])) # find siblings
        fam1 <- fam1[!fam1[,selfID] %in% sibs1[,selfID]] # remove siblings
        aunts1 <- unique(rbind(fam1[fam1[, motherID] == mother1[, motherID], ],
                               fam1[fam1[, motherID] == father1[, motherID], ],
                               fam1[fam1[, fatherID] == mother1[, fatherID], ],
                               fam1[fam1[, fatherID] == father1[, fatherID], ])) # find aunts/uncles
        fam1 <- fam1[!fam1[,selfID] %in% aunts1[,selfID]] # remove siblings
        kids1 <- unique(rbind(fam1[fam1[, motherID] == self1[, selfID], ],
                              fam1[fam1[, fatherID] == self1[, selfID], ])) # find children
        fam1 <- fam1[!fam1[,selfID] %in% kids1[,selfID]] # remove children
        self2 <- fam2[fam2[, selfID] == 1, ]  # repeat for comparison family
        fam2 <- fam2[fam2[, selfID] != 1, ]  # remove self
        mother2 <- fam2[fam2[, selfID] == self2[, motherID], ]  # find mother
        fam2 <- fam2[fam2[, selfID] != self2[, motherID], ]  # remove mother
        fam2 <- fam2[fam2[, selfID] != mother2[, motherID], ]  # remove mothers mother
        fam2 <- fam2[fam2[, selfID] != mother2[, fatherID], ]  # remove mothers father
        father2 <- fam2[fam2[, selfID] == self2[, fatherID], ]  # find father
        fam2 <- fam2[fam2[, selfID] != self2[, fatherID], ]  # remove father
        fam2 <- fam2[fam2[, selfID] != father2[, motherID], ]  # remove fathers mother
        fam2 <- fam2[fam2[, selfID] != father2[, fatherID], ]  # remove fathers father
        sibs2 <- unique(rbind(fam2[fam2[, motherID] == mother2[, selfID], ],
                              fam2[fam2[, fatherID] == father2[, selfID], ])) # find siblings
        fam2 <- fam2[!fam2[,selfID] %in% sibs2[,selfID]] # remove siblings
        aunts2 <- unique(rbind(fam2[fam2[, motherID] == mother2[, motherID], ],
                               fam2[fam2[, motherID] == father2[, motherID], ],
                               fam2[fam2[, fatherID] == mother2[, fatherID], ],
                               fam2[fam2[, fatherID] == father2[, fatherID], ])) # find aunts/uncles
        fam2 <- fam2[!fam2[,selfID] %in% aunts2[,selfID]] # remove siblings
        kids2 <- unique(rbind(fam2[fam2[, motherID] == self2[, selfID], ],
                              fam2[fam2[, fatherID] == self2[, selfID], ])) # find children
        fam2 <- fam2[!fam2[,selfID] %in% kids2[,selfID]] # remove children
        sibscore <- greedyConcordance(sibs1, sibs2, object$keyVars)  # try to concord siblings
        auntscore <- greedyConcordance(aunts1, aunts2, object$keyVars)  # try to concord aunts/uncles
        kidscore <- greedyConcordance(kids1, kids2, object$keyVars)  # try to concord children
        otherscore <- greedyConcordance(fam1, fam2, object$keyVars)  # try to concord all others
        as.numeric(sibscore) + as.numeric(auntscore) + as.numeric(kidscore) + as.numeric(otherscore)
      } else {as.numeric(Neighbors[id, "matchScore"])}
    }, silent = TRUE)
    ifelse(class(tsc) == "try-error", 0, tsc)
  })
  scoreVec <- as.numeric(selfScore) + as.numeric(motherScore) +
    as.numeric(fatherScore) + as.numeric(ffatherScore) +
    as.numeric(fmotherScore) + as.numeric(mfatherScore) +
    as.numeric(mmotherScore) + as.numeric(allScore)
  object <- set(object, "scoreVec", scoreVec)
  return(object)
}
