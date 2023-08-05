
# phyloraster

<!-- badges: start -->
  [![.travis](https://github.com/gabferreira/phyloraster/actions/workflows/.travis.yml/badge.svg)](https://github.com/gabferreira/phyloraster/actions/workflows/.travis.yml)
  <!-- badges: end -->

# phyloraster <a href="https://github.com/gabferreira/phyloraster"><img src="man/figures/logo.png" align="right" height="139" alt="phyloraster website" /></a>

[`phyloraster`](https://github.com/gabferreira/phyloraster) is an R package to calculate measures of endemism and evolutionary diversity using rasters of presence-absence as input, allowing to join the results derived from species distribution models (SDMs) with phylogenetic information.

## Installation

You can install the development version of *phyloraster* package from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("gabferreira/phyloraster")
```

## Steps to calculte phylogenetic diversity using ```phyloraster```

### First, load phylogenetic data and rasters of presence-absence for a set of species.

``` r 
library(phyloraster)
library(terra)
library(ape)
```

``` r
ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
```

### Now, let's prepare the dataset. The raster stack will be sorted according to the tree order and the branch lengths will be extracted from the tree for each species.

``` r
dataprep <- phyloraster::phylo.pres(ras, tree)
```

### Now, we are already able to calculate the phylogenetic diversity.

``` r
pd <- phyloraster::rast.pd(data$x, branch.length = data$branch.length)
pd
```

### The result can be visualized using the R `plot` function from the `terra` package.

``` r
terra::plot(pd, main = "Phylogenetic Endemism")
```

A vignette with other examples can be found loading:

``` r
browseVignettes("phyloraster")
```

