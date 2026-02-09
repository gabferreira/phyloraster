# Calculate range size for a set of species using a raster as input

This function calculate range size in square meters (by default) for all
cells that are not NA. The size of the cells is constant in degrees but
not in square meters, which was considered in the method applied to
calculate the area.

## Usage

``` r
range_size(x, cellSz, unit = "m", ...)
```

## Arguments

- x:

  SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for
  a set of species. The layers (species) will be sorted according to the
  tree order. See the phylo.pres function.

- cellSz:

  SpatRaster. A SpatRaster containing cellSize values. See
  [`cellSize`](https://rspatial.github.io/terra/reference/cellSize.html)

- unit:

  character. One of "m", "km", or "ha"

- ...:

  additional arguments to be passed passed down from a calling function.

## Value

vector

## Author

Gabriela Alves Ferreira and Neander Marcel Heming

## Examples

``` r
# \donttest{
x <- terra::rast(system.file("extdata", "rast.presab.tif",
package="phyloraster"))
range_size(x[[1:2]], cellSz <- terra::cellSize(x))
#> Litoria_revelata   Litoria_rothii 
#>       5021016045     196704341029 
# }
```
