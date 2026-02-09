# Calculate phylogenetic diversity for each raster cell

Calculate the sum of the branch length for species present in raster
data.

## Usage

``` r
.rast.pd.B(x, edge.path, branch.length, filename = "", ...)
```

## Arguments

- x:

  SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for
  a set of species. The layers (species) will be sorted according to the
  tree order. See the phylo.pres function.

- edge.path:

  matrix. Matrix representing the paths through the tree from root to
  each tip. See
  [`phylo.pres`](https://gabferreira.github.io/phyloraster/reference/phylo.pres.md)

- branch.length:

  numeric. A Named numeric vector of branch length for each species. See
  [`phylo.pres`](https://gabferreira.github.io/phyloraster/reference/phylo.pres.md)

- filename:

  character. Output filename

- ...:

  additional arguments passed for terra::app

## Value

SpatRaster

## References

Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity.
Biological conservation, 61(1), 1-10.

## Author

Neander Marcel Heming and Gabriela Alves-Ferreira
