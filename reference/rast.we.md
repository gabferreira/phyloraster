# Calculate weighted endemism for raster data

Calculate the weighted endemism for species present in raster data.

## Usage

``` r
rast.we(x, inv.R, filename = "", ...)
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

Laffan, S. W., Rosauer, D. F., Di Virgilio, G., Miller, J. T.,
González‐Orozco, C. E., Knerr, N., ... & Mishler, B. D. (2016).
Range‐weighted metrics of species and phylogenetic turnover can better
resolve biogeographic transition zones. Methods in Ecology and
Evolution, 7(5), 580-588.

Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J.,
VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and
endemism in conservation. In: Systematics and Conservation Evaluation
(eds Forey PL, Humphries CJ, Vane-Wright RI), p. 438. Oxford University
Press, Oxford.

Crisp, M., Laffan, S., Linder, H., Monro, A. (2001). Endemism in the
Australian flora. Journal of Biogeography, 28, 183–198.

## Author

Neander Marcel Heming and Gabriela Alves Ferreira

## Examples

``` r
# \donttest{
library(terra)
library(phyloraster)
x <- rast(system.file("extdata", "rast.presab.tif",
package="phyloraster"))
inv.R <- inv.range(x)
we <- rast.we(x, inv.R)
plot(we)

# }
```
