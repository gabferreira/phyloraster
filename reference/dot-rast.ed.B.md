# Calculate Evolutionary distinctiveness for each raster cell

This function calculates evolutionary distinctiveness according to the
fair-proportion index.

## Usage

``` r
.rast.ed.B(x, edge.path, branch.length, n.descen, filename = "", ...)
```

## Arguments

- x:

  SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for
  a set of species. The layers (species) must be sorted according to the
  tree order. See the phylo.pres function.

- edge.path:

  matrix. Matrix representing the paths through the tree from root to
  each tip. See
  [`phylo.pres`](https://gabferreira.github.io/phyloraster/reference/phylo.pres.md)

- branch.length:

  numeric. A Named numeric vector of branch length for each species. See
  [`phylo.pres`](https://gabferreira.github.io/phyloraster/reference/phylo.pres.md)

- n.descen:

  numeric. A Named numeric vector of number of descendants for each
  branch. See
  [`phylo.pres`](https://gabferreira.github.io/phyloraster/reference/phylo.pres.md)

- filename:

  character. Output filename

- ...:

  additional arguments passed for terra::app

## Value

SpatRaster

## References

Isaac, N. J., Turvey, S. T., Collen, B., Waterman, C. and Baillie, J. E.
(2007). Mammals on the EDGE: conservation priorities based on threat and
phylogeny. PLoS ONE 2, e296.

## Author

Neander Marcel Heming and Gabriela Alves-Ferreira
