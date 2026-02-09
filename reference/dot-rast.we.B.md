# Calculate weighted endemism for each raster cell

Calculate the sum of the inverse of the range size for species present
in each raster cell.

## Usage

``` r
.rast.we.B(x, inv.R, filename = "", ...)
```

## Arguments

- x:

  SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for
  a set of species. The layers (species) will be sorted according to the
  tree order. See the phylo.pres function.

- inv.R:

  SpatRaster. Inverse of range size. See
  [`inv.range`](https://gabferreira.github.io/phyloraster/reference/inv.range.md)

- filename:

  character. Output filename

- ...:

  additional arguments passed for terra::app

## Value

SpatRaster

## References

Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J.,
VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and
endemism in conservation. In: Systematics and Conservation Evaluation
(eds Forey PL, Humphries CJ, Vane-Wright RI), p. 438. Oxford University
Press, Oxford.

Crisp, M., Laffan, S., Linder, H., Monro, A. (2001). Endemism in the
Australian flora. Journal of Biogeography, 28, 183–198.

## Author

Neander Marcel Heming and Gabriela Alves-Ferreira
