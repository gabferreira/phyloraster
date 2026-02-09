# Calculate phylogenetic endemism for raster data

Calculate the sum of the inverse of the range size multiplied by the
branch length for the species present in raster data.

## Usage

``` r
rast.pe(
  x,
  tree,
  inv.R,
  branch.length,
  full_tree_metr = TRUE,
  filename = "",
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

- inv.R:

  SpatRaster. Inverse of range size. See
  [`inv.range`](https://gabferreira.github.io/phyloraster/reference/inv.range.md)

- branch.length:

  numeric. A Named numeric vector of branch length for each species. See
  [`phylo.pres`](https://gabferreira.github.io/phyloraster/reference/phylo.pres.md)

- full_tree_metr:

  logical. Whether edge.path, branch length and number of descendants
  should be calculated with the full (TRUE) or the prunned tree (FALSE).
  The default is TRUE.

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

Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C. and
Cook, L. G. (2009). Phylogenetic endemism: a new approach for
identifying geographical concentrations of evolutionary history.
Molecular ecology, 18(19), 4061-4072.

## Author

Gabriela Alves-Ferreira and Neander Marcel Heming

## Examples

``` r
# \donttest{
library(terra)
library(phyloraster)
x <- rast(system.file("extdata", "rast.presab.tif",
package = "phyloraster"))
tree <- ape::read.tree(system.file("extdata", "tree.nex",
package = "phyloraster"))
pe <- rast.pe(x = x[[1:3]], tree)
#> Warning: Some species in the phylogeny 'tree' are missing from the
#>                   SpatRaster 'x' and were dropped: Litoria_dorsalis, Litoria_rubella, Litoria_nigrofrenata, Litoria_nasuta, Litoria_tornieri, Litoria_inermis, Litoria_pallida, Litoria_latopalmata, Litoria_bicolor, Litoria_fallax, Litoria_genimaculata, Litoria_andiirrmalin, Litoria_wilcoxii, Litoria_jungguy, Litoria_caerulea, Litoria_gracilenta, Litoria_chloris, Litoria_xanthomera, Cyclorana_brevipes, Cyclorana_novaehollandiae, Cyclorana_manya, Cyclorana_cultripes, Litoria_alboguttata, Cyclorana_longipes, Nyctimystes_dayi, Litoria_nannotis, Litoria_lorica, Litoria_rheocola, Litoria_nyakalensis, Litoria_infrafrenata
plot(pe)

# }
```
