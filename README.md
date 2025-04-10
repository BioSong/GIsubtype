
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GIsubtype

<!-- badges: start -->
<!-- badges: end -->

The goal of GIsubtype is to panGI cancer subtype.

## Installation

You can install the development version of GIsubtype from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("BioSong/GIsubtype")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(GIsubtype)
## basic example code
res <- GIclassifier(ExamExp,idType="SYMBOL")
#> Checking input dataset and parameters......
#> Transforming the input gene expression profile into dataframe format.
#> Making molecular subtype prediction......
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
