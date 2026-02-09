# Load an example dataset with presence-absence data of 33 tree frogs and a phylogenetic tree for this species

This function load a phylogenetic tree, a raster and a data.frame with
presence-absence of 33 Australian tree frogs from Rosauer (2017). We
also provide distribution shapefiles for ten species according to the
IUCN.

## Usage

``` r
load.data.rosauer()
```

## Source

Rosauer, 2017. Available on:
[Github](https://github.com/DanRosauer/phylospatial/tree/master/PhyloEndemism_in_R/Tree%20Frog%20Data/)

IUCN. 2022. The IUCN Red List of Threatened Species (spatial data).
Version 2022-1. [IUCN](https://www.iucnredlist.org)

## Value

data.frame, SpatRaster, SpatVector and phylo
