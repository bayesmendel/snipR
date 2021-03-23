#' Simulated family pedigrees. 
#'
#' A dataset containing (possibly duplicated) family pedigrees, 
#' with each family consisting of a reported 'self' individual
#' and various cancer risk indicators/status, as well as a 
#' pseudo-random collection of reported family members
#' and their associated indicators. Allele sharing between generations
#' follow basic Mendelian inheritance with parameters provided by
#' BayesMendel.
#'
#' @docType data
#'
#' @usage data(pedigrees)
#' @format An object of class \code{'cross'}; see \code{\link[qtl]{read.cross}}.
# 
#' @format A data frame with 1033 rows and 29 variables:
#' \describe{
#'   \item{ID}{Family member ID (self=1)}
#'   \item{Gender}{Reported gender}
#'   \item{FatherID}{ID for reported father}
#'   \item{MotherID}{ID for reported mother}
#'   \item{AffectedBreast}{Breast cancer affliction status (binary)}
#'   \item{AgeBreast}{Age at onset or current age if AffectedBreast=0}
#'   \item{Twins}{Does individual have a twin?}
#'   \item{BRCA1}{BRCA1 gene status}
#'   \item{Duplicate}{Is pedigree duplicated?}
#'   \item{nDuplicates}{Number of duplicate pedigrees present}
#'   \item{FamilyID}{Pedigree unique identifier (typically unknown)}
#'   \item{requestID}{Pedigree ID (unmatched identifier)}
#'   \item{senderIP}{Simulated IP address used}
#'   ...
#' }
#' @keywords datasets pedigrees data 
#' @source \url{http://bcb.dfci.harvard.edu/bayesmendel/software.php/}
#'
#' @examples
#' data(pedigrees)
#' self <- pedigrees[pedigrees$ID==1,]
#' \donttest{hist(self$nDuplicates,
#'    breaks=seq(0,max(sself$nDuplicates),by=1),
#'    main='Distribution of Duplications',xlab='Duplicate Count',
#'    ylab='Frequency',right=FALSE)}
"pedigrees"
