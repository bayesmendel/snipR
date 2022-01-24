<!-- README.md is generated from README.Rmd. Please edit that file -->

# snipR

snipR implements the Sorted NeIghborhoods for Pedigrees (SNIP) algorithm
for deduplicating pedigree data.

## Installation

You can install snipR from github with:

``` r
# install.packages("devtools")
devtools::install_github("bayesmendel/snipR")
```

If you encounter a bug, please file a minimal reproducible example on
[github](https://github.com/bayesmendel/snipR/issues).

## Learning snipR

To get started, please read the intro vignette:
`vignette("Introduction", package = "snipR")`.

## Key data structures

The base S3 object in bsnR is a *Neighbors*, a contained for the raw
data and user specified keys. At each step of the analysis this object
will be promoted, allowing subsequent steps to be performed.
