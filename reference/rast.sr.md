# Calculate species richness for raster data

Calculate the species richness for raster data.

## Usage

``` r
rast.sr(x, filename = "", ...)
```

## Arguments

- x:

  SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for
  a set of species.

- filename:

  character. Output filename.

- ...:

  additional arguments to be passed passed down from a calling function.

## Value

SpatRaster

## Author

Gabriela Alves Ferreira and Neander Marcel Heming

## Examples

``` r
# \donttest{
x <- terra::rast(system.file("extdata", "rast.presab.tif",
package="phyloraster"))
rse <- phyloraster::rast.sr(x)
terra::plot(rse)

# }
```
