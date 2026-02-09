# Calculate phylogenetic diversity for raster data

Calculate the sum of the branch length for species present in each cell
of the raster.

## Usage

``` r
rast.pd(
  x,
  tree,
  edge.path,
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

- edge.path:

  matrix. Matrix representing the paths through the tree from root to
  each tip. See
  [`phylo.pres`](https://gabferreira.github.io/phyloraster/reference/phylo.pres.md)

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

Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity.
Biological conservation, 61(1), 1-10.

## Author

Neander Marcel Heming and Gabriela Alves-Ferreira

## Examples

``` r
# \donttest{
library(terra)
library(phyloraster)
x <- rast(system.file("extdata", "rast.presab.tif",
package="phyloraster"))
tree <- ape::read.tree(system.file("extdata", "tree.nex",
package="phyloraster"))
data <- phylo.pres(x[[1:3]], tree)
#> Warning: Some species in the phylogeny 'tree' are missing from the
#>                   SpatRaster 'x' and were dropped: Litoria_dorsalis, Litoria_rubella, Litoria_nigrofrenata, Litoria_nasuta, Litoria_tornieri, Litoria_inermis, Litoria_pallida, Litoria_latopalmata, Litoria_bicolor, Litoria_fallax, Litoria_genimaculata, Litoria_andiirrmalin, Litoria_wilcoxii, Litoria_jungguy, Litoria_caerulea, Litoria_gracilenta, Litoria_chloris, Litoria_xanthomera, Cyclorana_brevipes, Cyclorana_novaehollandiae, Cyclorana_manya, Cyclorana_cultripes, Litoria_alboguttata, Cyclorana_longipes, Nyctimystes_dayi, Litoria_nannotis, Litoria_lorica, Litoria_rheocola, Litoria_nyakalensis, Litoria_infrafrenata
pd <- rast.pd(data$x, data$tree)
plot(pd)

# }
```
