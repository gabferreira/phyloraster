# Transform a data.frame to raster

The function transforms a data.frame or a matrix of presence- absence in
a raster of distribution.

## Usage

``` r
df2rast(x, CRS = "+proj=longlat +datum=WGS84", ...)
```

## Arguments

- x:

  data.frame. A data.frame or matrix with species names in columns and
  sites in rows. The first two columns must provide longitude and
  latitude, respectively.

- CRS:

  character. Description of the Coordinate Reference System (map
  projection) in PROJ.4.

- ...:

  additional arguments to be passed passed down from a calling function.

## Value

SpatRaster

## Examples

``` r
# \donttest{
dat <- phyloraster::load.data.rosauer()
df2rast(dat$presab, crs = "+proj=longlat +datum=WGS84 +ellps=WGS84
+towgs84=0,0,0")
#> class       : SpatRaster 
#> size        : 90, 68, 33  (nrow, ncol, nlyr)
#> resolution  : 0.1, 0.1  (x, y)
#> extent      : 144.0157, 150.8157, -23.044, -14.044  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=longlat +datum=WGS84 +no_defs 
#> source(s)   : memory
#> names       : Litor~elata, Litor~othii, Litor~stris, Litor~salis, Litor~bella, Litor~ermis, ... 
#> min values  :           0,           0,           0,           0,           0,           0, ... 
#> max values  :           1,           1,           1,           1,           1,           1, ... 
# }
```
