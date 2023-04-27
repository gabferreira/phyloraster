
# phylogrid

<!-- badges: start -->
[![Travis build status](https://travis-ci.com/gabferreira/phylogrid.svg?branch=master)](https://travis-ci.com/gabferreira/phylogrid)
[![.travis](https://github.com/gabferreira/phylogrid/actions/workflows/.travis.yml/badge.svg)](https://github.com/gabferreira/phylogrid/actions/workflows/.travis.yml)
<!-- badges: end -->

The goal of *phylogrid* package is to calculate phylogenetic indices (such as phylogenetic diversity (PD. Faith 1992), evolutionary distinctiveness (Isaac et al. 2007; Laffan et al. 2016), phylogenetic endemism (PE. Rosauer et al. 2009, Laffan et al. 2016), weighted endemism (WE. Laffan et al. 2016)) for presence-absence rasters and return a raster object as an output. See more details on the package vignette.

## Installation

You can install the development version of *phylogrid* package from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("gabferreira/phylogrid")
```

## Steps to calculte PD, ED, PE and WE using ```phylogrid```

### First, load phylogenetic and spatial data

``` r 
library(phylogrid)
library(terra)
library(ape)
```

``` r
ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
```

### Now, let's prepare the dataset! The stack will be sorted according to the tree order and the branch lengths will be extracted for each species.

``` r
dataprep <- phylogrid::phylo.pres(ras, tree)
```

### Now, we are already able to calculate PD, ED, PE and WE!!

``` r
pd <- phylogrid::rast.pd(data$x, data$branch.length)
pe <- phylogrid::rast.pe(data$x, data$branch.length)
we <- phylogrid::rast.we(ras)
ed <- phylogrid::rast.ed(data$x, data$branch.length, data$n.descendants)
```

### The result can be visualized using the R `plot` function from the `terra` package.

``` r
terra::plot(pe, main = "Phylogenetic Endemism")
```

