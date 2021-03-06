% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blockedSN.R
\name{blockedSN}
\alias{blockedSN}
\alias{blockedSN.default}
\alias{blockedSN.Neighbors}
\alias{blockedSN.Blocks}
\alias{blockedSN.Scores}
\title{Blocked sorted neighbors generic.}
\usage{
blockedSN(object, blockVar, repSN, windowSN, keyLength)

\method{blockedSN}{default}(object, blockVar, repSN, windowSN, keyLength)

\method{blockedSN}{Neighbors}(object, blockVar, repSN, windowSN, keyLength)

\method{blockedSN}{Blocks}(object, blockVar, repSN, windowSN, keyLength)

\method{blockedSN}{Scores}(object, blockVar, repSN, windowSN, keyLength)
}
\arguments{
\item{object}{\code{\link{Neighbors}} or \code{\link{Blocks}} object to perform SNIP algorithm on.}

\item{blockVar}{Character value.  Blocking variable.}

\item{repSN}{Integer value. How much iterations of sorted neighbors to perform?}

\item{windowSN}{Integer value. Size of sliding window to use during sorted neighbors.}

\item{keyLength}{Numeric value or vector. How many keyVars to concatenate per sort key.}
}
\value{
An object of class \code{\link{Blocks}} containing the neighbors found
and keys used during the blocked sorted neighbors iteration.
}
\description{
\code{blockedSN} returns a \code{\link{Blocks}} object created from a single
iteration of the SNIP algorithm.
}
