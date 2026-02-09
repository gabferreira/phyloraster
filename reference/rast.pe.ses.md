# Standardized effect size for Phylogenetic endemism

Calculates the standardized effect size for phylogenetic endemism. See
Details for more information.

## Usage

``` r
rast.pe.ses(
  x,
  tree,
  branch.length,
  branch.length.alt,
  inv.R,
  full_tree_metr = TRUE,
  spat_alg = "bootspat_str",
  spat_alg_args = list(rprob = NULL, rich = NULL, fr_prob = NULL),
  metric = c("pe", "pe.alt", "rpe", "all")[4],
  aleats = 10,
  cores = 1,
  filename = "",
  overwrite = TRUE,
  ...
)
```

## Arguments

- x:

  SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for
  a set of species. The layers (species) will be sorted according to the
  tree order. See the phylo.pres function.

- tree:

  phylo. A dated tree.

- branch.length:

  numeric. A Named numeric vector of branch length for each species. See
  [`phylo.pres`](https://gabferreira.github.io/phyloraster/reference/phylo.pres.md)

- branch.length.alt:

  numeric. Branch length calculated by using an alternative phylogeny
  with non-zero branch lengths converted to a constant value (1) and
  rescaled so the sum of all branch lengths is 1.

- inv.R:

  SpatRaster. Inverse of range size. See
  [`inv.range`](https://gabferreira.github.io/phyloraster/reference/inv.range.md)

- full_tree_metr:

  logical. Whether edge.path, branch length and number of descendants
  should be calculated with the full (TRUE) or the prunned tree (FALSE).
  The default is TRUE.

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

- metric:

  character. Names of biodiversity metrics to calculate (pe, pe_alt,
  rpe, all). See details.

- aleats:

  positive integer. A positive integer indicating how many times the
  calculation should be repeated.

- cores:

  positive integer. If `cores > 1`, a 'parallel' package cluster with
  that many cores is created and used. You can also supply a cluster
  object. Ignored for functions that are implemented by terra in C++
  (see under fun)

- filename:

  character. Output filename

- overwrite:

  logical. If TRUE, filename is overwritten

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
cells. Biodiversity metrics available are:

- pe: Phylogenetic endemism (Rosauer et al., 2009)

- pe.alt: Alternate Phylogenetic endemism (Mishler et al., 2014)

- rpe: Relative Phylogenetic endemism (Mishler et al., 2014)

- all: Calculate all available metrics Alternate phylogenetic endemism
  (PE.alt, Mishler et al., 2014) is calculated using an alternate
  phylogeny with non-zero branch lengths converted to a constant value
  (here we use 1) and rescaled so the sum of all branch lengths is 1.
  Relative phylogenetic endemism (RPE, Mishler et al., 2014) is the
  ratio of phylogenetic endemism (PE, Rosauer et al., 2009) measured on
  the original tree versus PE measured on a alternate tree (PE.alt).

## References

Gotelli, N. J. 2000. Null model analysis of species co-occurrence
patterns. Ecology 81: 2606–2621.

Heming, N. M., Mota, F. M. M. and Alves-Ferreira, G. 2023. SESraster:
raster randomization for null hypothesis testing.
https://CRAN.R-project.org/package=SESraster.

Mishler, B. D., Knerr, N., González-Orozco, C. E., Thornhill, A. H.,
Laffan, S. W. and Miller, J. T. 2014. Phylogenetic measures of
biodiversity and neo- and paleo-endemism in Australian Acacia. – Nat.
Commun. 5: 4473.

Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C., &
Cook, L. G. (2009). Phylogenetic endemism: a new approach for
identifying geographical concentrations of evolutionary history.
Molecular ecology, 18(19), 4061-4072.

## See also

[`phylo.pres`](https://gabferreira.github.io/phyloraster/reference/phylo.pres.md),
[`inv.range`](https://gabferreira.github.io/phyloraster/reference/inv.range.md),
[`geo.phylo.ses`](https://gabferreira.github.io/phyloraster/reference/geo.phylo.ses.md),
[`rast.ed.ses`](https://gabferreira.github.io/phyloraster/reference/rast.ed.ses.md),
[`rast.pd.ses`](https://gabferreira.github.io/phyloraster/reference/rast.pd.ses.md),
[`rast.we.ses`](https://gabferreira.github.io/phyloraster/reference/rast.we.ses.md),
`rast.pe.ses`,
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
library(phyloraster)
library(SESraster)
x <- terra::rast(system.file("extdata", "rast.presab.tif",
package="phyloraster"))
tree <- ape::read.tree(system.file("extdata", "tree.nex",
package="phyloraster"))
data <- phylo.pres(x[[1:3]], tree)
#> Warning: Some species in the phylogeny 'tree' are missing from the
#>                   SpatRaster 'x' and were dropped: Litoria_dorsalis, Litoria_rubella, Litoria_nigrofrenata, Litoria_nasuta, Litoria_tornieri, Litoria_inermis, Litoria_pallida, Litoria_latopalmata, Litoria_bicolor, Litoria_fallax, Litoria_genimaculata, Litoria_andiirrmalin, Litoria_wilcoxii, Litoria_jungguy, Litoria_caerulea, Litoria_gracilenta, Litoria_chloris, Litoria_xanthomera, Cyclorana_brevipes, Cyclorana_novaehollandiae, Cyclorana_manya, Cyclorana_cultripes, Litoria_alboguttata, Cyclorana_longipes, Nyctimystes_dayi, Litoria_nannotis, Litoria_lorica, Litoria_rheocola, Litoria_nyakalensis, Litoria_infrafrenata
t <- rast.pe.ses(x = data$x, data$tree, aleats = 99, metric = "all")
#> Please cite SESraster when using spatial null models.
#>           See: citation(SESraster)
plot(t)

# }
```
