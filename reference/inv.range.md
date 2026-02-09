# Calculate the inverse of range size

Get range size in square kilometers for all cells that are not NA, the
inverse of range size and the inverse of range size multiplied by branch
length for multiple species using a raster of presence-absence.

## Usage

``` r
inv.range(x, filename = "", overwrite = FALSE, ...)
```

## Arguments

- x:

  SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for
  a set of species. The layers (species) will be sorted according to the
  tree order. See the phylo.pres function.

- filename:

  character. Output filename

- overwrite:

  logical. If TRUE, filename is overwritten

- ...:

  additional arguments to be passed passed down from a calling function.

## Value

SpatRaster and numeric

## Author

Neander Marcel Heming and Gabriela Alves-Ferreira

## Examples

``` r
# \donttest{
# calculating the inverse of range size
x <- terra::rast(system.file("extdata", "rast.presab.tif",
package="phyloraster"))
inv.range(x[[5]])
#> class       : SpatRaster 
#> size        : 90, 68, 1  (nrow, ncol, nlyr)
#> resolution  : 0.1, 0.1  (x, y)
#> extent      : 144.0157, 150.8157, -23.044, -14.044  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#> source(s)   : memory
#> varname     : rast.presab 
#> name        : Litoria_rubella 
#> min value   :    0.0005237743 
#> max value   :    0.0005511653 
# }
```
