
# phyloraster

<!-- badges: start -->
[![Travis build status](https://travis-ci.com/gabferreira/phyloraster.svg?branch=master)](https://travis-ci.com/gabferreira/phyloraster)
[![.travis](https://github.com/gabferreira/phyloraster/actions/workflows/.travis.yml/badge.svg)](https://github.com/gabferreira/phyloraster/actions/workflows/.travis.yml)
<!-- badges: end -->

# phyloraster <a href="https://github.com/gabferreira/phyloraster"><img src="man/figures/logo.png" align="right" height="139" alt="phyloraster website" /></a>

The goal of *phyloraster* package is to calculate phylogenetic indices (such as phylogenetic diversity (PD. Faith 1992), evolutionary distinctiveness (Isaac et al. 2007; Laffan et al. 2016), phylogenetic endemism (PE. Rosauer et al. 2009, Laffan et al. 2016), weighted endemism (WE. Laffan et al. 2016)) for presence-absence rasters and return a raster object as an output. See more details on the package vignette.

## Installation

You can install the development version of *phyloraster* package from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("gabferreira/phyloraster")
```

## Steps to calculte PD, ED, PE and WE using ```phyloraster```

### First, load phylogenetic and spatial data

``` r 
library(phyloraster)
library(terra)
library(ape)
```

``` r
ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
```

### Now, let's prepare the dataset! The stack will be sorted according to the tree order and the branch lengths will be extracted for each species.

``` r
dataprep <- phyloraster::phylo.pres(ras, tree)
```

### Now, we are already able to calculate PD, ED, PE and WE!!

``` r
pd <- phyloraster::rast.pd(data$x, data$branch.length)
pe <- phyloraster::rast.pe(data$x, data$branch.length)
we <- phyloraster::rast.we(ras)
ed <- phyloraster::rast.ed(data$x, data$branch.length, data$n.descendants)
```

### The result can be visualized using the R `plot` function from the `terra` package.

``` r
terra::plot(pe, main = "Phylogenetic Endemism")
```

