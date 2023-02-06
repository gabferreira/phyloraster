
# phylogrid

<!-- badges: start -->
[![Travis build status](https://travis-ci.com/gabferreira/phylogrid.svg?branch=master)](https://travis-ci.com/gabferreira/phylogrid)
[![.travis](https://github.com/gabferreira/phylogrid/actions/workflows/.travis.yml/badge.svg)](https://github.com/gabferreira/phylogrid/actions/workflows/.travis.yml)
<!-- badges: end -->

The goal of *phylogrid* package is to calculate phylogenetic indices (such as phylogenetic diversity (PD. Faith 1992), phylogenetic endemism (PE. Rosauer et al. 2009), weighted endemism (WE)) for presence-absence rasters and return a raster object as an output.

## Installation

You can install the development version of *phylogrid* package from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("gabferreira/phylogrid")
```

## Steps to calculte PD, PE and WE using ```phylogrid```

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
dataprep <- phylogrid::phylo.pres(pres.stack = ras, tree = tree)
```

### After that, we will calculate the inverse of the range size and the inverse of the range size multiplied by the length of the branches. If you want to save these rasters directly to your computer, provide a path in the filename argument.

``` r
range <- phylogrid::inv.range(x = dataprep$x, branch.length = dataprep$branch.length, filename = NULL)
```

### Now, we are already able to calculate PD, PE and WE!!

``` r
pg <- phylogrid::geo.phylo(x = dataprep$x, area.inv = range$inv.R,
                           area.tips = range$LR, branch.length = dataprep$branch.length, filename = NULL)
```

