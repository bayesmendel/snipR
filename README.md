<!-- README.md is generated from README.Rmd. Please edit that file -->
bsnR
====

bsnR implements the blocked sorted neighbors algorithm for discovering potential duplicate pairs. The analysis is broken down into a series of steps:

-   Initialize Neighbors object

-   Perform BSN and collect possible duplicate pairs

-   Score duplicate pairs and retain those above specified threshold

-   Export, repeat, or verify accuracy of algorithm

Installation
------------

You can install bsnR from github with:

``` r
# install.packages("devtools")
devtools::install_github("mPloenzke/bsnR")
```

If you encounter a bug, please file a minimal reproducible example on [github](https://github.com/mPloenzke/bsnR/issues).

Learning bsnR
-------------

To get started, please read the intro vignette: `vignette("Introduction", package = "bsnR")`.

Key data structures
-------------------

The base S3 object in bsnR is a *Neighbors*, a contained for the raw data and user specified keys. At each step of the analysis this object will be promoted, allowing subsequent steps to be performed.
