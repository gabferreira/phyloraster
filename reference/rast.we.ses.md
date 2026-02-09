# Calculate weighted endemism standardized for species richness

Calculates the standardized effect size for weighted endemism. See
Details for more information.

## Usage

``` r
rast.we.ses(
  x,
  inv.R,
  spat_alg = "bootspat_str",
  spat_alg_args = list(rprob = NULL, rich = NULL, fr_prob = NULL),
  aleats = 10,
  filename = "",
  ...
)
```

## Arguments

- x:

  SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for
  a set of species. The layers (species) will be sorted according to the
  tree order. See the phylo.pres function.

- inv.R:

  SpatRaster. Inverse of range size. See
  [`inv.range`](https://gabferreira.github.io/phyloraster/reference/inv.range.md)

- spat_alg:

  A function with the algorithm implementing the desired randomization
  method. It must work with SpatRaster objects. See examples. Example of
  functions that work are:
  [`bootspat_naive`](https://hemingnm.github.io/SESraster/reference/bootspat_naive.html),
  [`bootspat_str`](https://hemingnm.github.io/SESraster/reference/bootspat_str.html),
  [`bootspat_ff`](https://hemingnm.github.io/SESraster/reference/bootspat_ff.html).

- spat_alg_args:

  List of arguments passed to the randomization method chosen in
  'spat_alg'. See
  [`bootspat_naive`](https://hemingnm.github.io/SESraster/reference/bootspat_naive.html),
  [`bootspat_str`](https://hemingnm.github.io/SESraster/reference/bootspat_str.html),
  [`bootspat_ff`](https://hemingnm.github.io/SESraster/reference/bootspat_ff.html)

- aleats:

  positive integer. A positive integer indicating how many times the
  calculation should be repeated.

- filename:

  character. Output filename

- ...:

  additional arguments passed for terra::app

## Value

SpatRaster. The function returns the observed value of the metric, the
mean of the simulations calculated over n times, the standard deviation
of the simulations, the standardized effect size (SES) for the metric,
and the p-values.

## Details

The dependency ‘SESraster’ is used to calculate the null models. This
package currently implements six algorithms to randomize binary species
distribution with several levels of constraints: SIM1, SIM2, SIM3, SIM5,
SIM6 and SIM9 (sensu Gotelli 2000). The methods implemented in
‘SESraster’ are based on how species (originally rows) and sites
(originally columns) are treated (i.e. fixed, equiprobable, or
proportional sums) (Gotelli 2000). By default, the ‘phyloraster’ uses
the function bootspat\_ str() from the ‘SESraster’ package to conduct
the randomizations, but the user is free to choose any of the other
methods mentioned above through the spat_alg argument in the \*.ses()
functions of the ‘phyloraster’ package. The bootspat_str() is equivalent
to the SIM5 (proportional-fixed) method of Gotelli (2000), which
partially relaxes the spatial structure of species distributions, but
keeps the spatial structure of the observed richness pattern across
cells.

## References

Gotelli, N. J. 2000. Null model analysis of species co-occurrence
patterns. – Ecology 81: 2606–2621.

Heming, N. M., Mota, F. M. M. and Alves-Ferreira, G. 2023. SESraster:
raster randomization for null hypothesis testing.
https://CRAN.R-project.org/package=SESraster.

## See also

[`phylo.pres`](https://gabferreira.github.io/phyloraster/reference/phylo.pres.md),
[`inv.range`](https://gabferreira.github.io/phyloraster/reference/inv.range.md),
[`geo.phylo.ses`](https://gabferreira.github.io/phyloraster/reference/geo.phylo.ses.md),
[`rast.ed.ses`](https://gabferreira.github.io/phyloraster/reference/rast.ed.ses.md),
[`rast.pd.ses`](https://gabferreira.github.io/phyloraster/reference/rast.pd.ses.md),
`rast.we.ses`,
[`rast.pe.ses`](https://gabferreira.github.io/phyloraster/reference/rast.pe.ses.md),
[`bootspat_str`](https://hemingnm.github.io/SESraster/reference/bootspat_str.html),
[`bootspat_naive`](https://hemingnm.github.io/SESraster/reference/bootspat_naive.html),
[`bootspat_ff`](https://hemingnm.github.io/SESraster/reference/bootspat_ff.html),
[`SESraster`](https://hemingnm.github.io/SESraster/reference/SESraster.html)

## Author

Gabriela Alves-Ferreira and Neander Heming

## Examples

``` r
# \donttest{
library(terra)
library(SESraster)
x <- terra::rast(system.file("extdata", "rast.presab.tif",
package="phyloraster"))
t <- rast.we.ses(x[[1:10]], aleats = 3)
#> Please cite SESraster when using spatial null models.
#>           See: citation(SESraster)
plot(t)

# }
```
