# Prepare rasters and phylogenetic tree to run community metrics

Reorder a stack of rasters of species distribution to match the order of
the tips of the tree, and get branch length and number of descendants
for each species to calculate diversity metrics using
phyloraster::geo.phylo(). The branch length and the number of
descendants can be calculated based on the full tree or the raster based
tree subset. The names must be the same in the phylogenetic tree and in
the raster for the same species. For example, if you have the name
"Leptodactylus_latrans" in the raster and "Leptodactylus latrans" in the
tree, the function will not work. The same goes for uppercase and
lowercase letters.

## Usage

``` r
phylo.pres(x, tree, full_tree_metr = TRUE, ...)
```

## Arguments

- x:

  SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for
  a set of species.

- tree:

  phylo. A dated tree.

- full_tree_metr:

  logical. Whether edge.path, branch length and number of descendants
  should be calculated with the full (TRUE) or the prunned tree (FALSE).
  The default is TRUE.

- ...:

  additional arguments to be passed passed down from a calling function.

## Value

Returns a list containing a SpatRaster reordered according to the order
that the species appear in the phylogenetic tree, a subtree containing
only the species that are in the stack of rasters and finally two named
numerical vectors containing the branch length and the number of
descendants of each species.

## Author

Neander Marcel Heming and Gabriela Alves Ferreira

## Examples

``` r
# \donttest{
library(phyloraster)
x <- terra::rast(system.file("extdata", "rast.presab.tif",
package="phyloraster"))
tree <- ape::read.tree(system.file("extdata", "tree.nex",
package="phyloraster"))
phylo.pres(x[[1:3]], tree, full_tree_metr = TRUE)
#> Warning: Some species in the phylogeny 'tree' are missing from the
#>                   SpatRaster 'x' and were dropped: Litoria_dorsalis, Litoria_rubella, Litoria_nigrofrenata, Litoria_nasuta, Litoria_tornieri, Litoria_inermis, Litoria_pallida, Litoria_latopalmata, Litoria_bicolor, Litoria_fallax, Litoria_genimaculata, Litoria_andiirrmalin, Litoria_wilcoxii, Litoria_jungguy, Litoria_caerulea, Litoria_gracilenta, Litoria_chloris, Litoria_xanthomera, Cyclorana_brevipes, Cyclorana_novaehollandiae, Cyclorana_manya, Cyclorana_cultripes, Litoria_alboguttata, Cyclorana_longipes, Nyctimystes_dayi, Litoria_nannotis, Litoria_lorica, Litoria_rheocola, Litoria_nyakalensis, Litoria_infrafrenata
#> $x
#> class       : SpatRaster 
#> size        : 90, 68, 3  (nrow, ncol, nlyr)
#> resolution  : 0.1, 0.1  (x, y)
#> extent      : 144.0157, 150.8157, -23.044, -14.044  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#> source      : rast.presab.tif 
#> names       : Litoria_revelata, Litoria_rothii, Litoria_longirostris 
#> min values  :                0,              0,                    0 
#> max values  :                1,              1,                    1 
#> 
#> $tree
#> 
#> Phylogenetic tree with 3 tips and 1 internal node.
#> 
#> Tip labels:
#>   Litoria_revelata, Litoria_rothii, Litoria_longirostris
#> 
#> Unrooted; includes branch length(s).
#> 
#> $edge.path
#>                      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
#> Litoria_revelata        1    1    1    1    0    0    0    0    0     0     0
#> Litoria_rothii          1    1    1    0    1    0    0    0    0     0     0
#> Litoria_longirostris    1    1    1    0    0    1    1    0    0     0     0
#>                      [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20]
#> Litoria_revelata         0     0     0     0     0     0     0     0     0
#> Litoria_rothii           0     0     0     0     0     0     0     0     0
#> Litoria_longirostris     0     0     0     0     0     0     0     0     0
#>                      [,21] [,22] [,23] [,24] [,25] [,26] [,27] [,28] [,29]
#> Litoria_revelata         0     0     0     0     0     0     0     0     0
#> Litoria_rothii           0     0     0     0     0     0     0     0     0
#> Litoria_longirostris     0     0     0     0     0     0     0     0     0
#>                      [,30] [,31] [,32] [,33] [,34] [,35] [,36] [,37] [,38]
#> Litoria_revelata         0     0     0     0     0     0     0     0     0
#> Litoria_rothii           0     0     0     0     0     0     0     0     0
#> Litoria_longirostris     0     0     0     0     0     0     0     0     0
#>                      [,39] [,40] [,41] [,42] [,43] [,44] [,45] [,46] [,47]
#> Litoria_revelata         0     0     0     0     0     0     0     0     0
#> Litoria_rothii           0     0     0     0     0     0     0     0     0
#> Litoria_longirostris     0     0     0     0     0     0     0     0     0
#>                      [,48] [,49] [,50] [,51] [,52] [,53] [,54] [,55] [,56]
#> Litoria_revelata         0     0     0     0     0     0     0     0     0
#> Litoria_rothii           0     0     0     0     0     0     0     0     0
#> Litoria_longirostris     0     0     0     0     0     0     0     0     0
#>                      [,57] [,58]
#> Litoria_revelata         0     0
#> Litoria_rothii           0     0
#> Litoria_longirostris     0     0
#> 
#> $branch.length
#>  [1] 0.173157 0.072386 0.175442 0.589016 0.589015 0.193871 0.395144 0.395144
#>  [9] 0.589015 0.273621 0.490836 0.201929 0.288907 0.092767 0.196139 0.066625
#> [17] 0.129514 0.129514 0.196139 0.325841 0.511002 0.511002 0.011052 0.095197
#> [25] 0.398434 0.020969 0.017805 0.466543 0.043192 0.423350 0.322782 0.100568
#> [33] 0.100568 0.093810 0.390538 0.268077 0.122460 0.068385 0.054075 0.054075
#> [41] 0.251337 0.253980 0.253980 0.043002 0.210977 0.210978 0.032253 0.178724
#> [49] 0.178724 0.154130 0.749621 0.401051 0.348570 0.348570 0.081570 0.267000
#> [57] 0.267000 0.998948
#> 
#> $branch.length.alt
#>  [1] 0.01724138 0.01724138 0.01724138 0.01724138 0.01724138 0.01724138
#>  [7] 0.01724138 0.01724138 0.01724138 0.01724138 0.01724138 0.01724138
#> [13] 0.01724138 0.01724138 0.01724138 0.01724138 0.01724138 0.01724138
#> [19] 0.01724138 0.01724138 0.01724138 0.01724138 0.01724138 0.01724138
#> [25] 0.01724138 0.01724138 0.01724138 0.01724138 0.01724138 0.01724138
#> [31] 0.01724138 0.01724138 0.01724138 0.01724138 0.01724138 0.01724138
#> [37] 0.01724138 0.01724138 0.01724138 0.01724138 0.01724138 0.01724138
#> [43] 0.01724138 0.01724138 0.01724138 0.01724138 0.01724138 0.01724138
#> [49] 0.01724138 0.01724138 0.01724138 0.01724138 0.01724138 0.01724138
#> [55] 0.01724138 0.01724138 0.01724138 0.01724138
#> 
#> $n.descendants
#> [1] 0 1 1 1 3
#> 

# using the prunned tree
phylo.pres(x[[1:3]], tree, full_tree_metr = FALSE)
#> Warning: Some species in the phylogeny 'tree' are missing from the
#>                   SpatRaster 'x' and were dropped: Litoria_dorsalis, Litoria_rubella, Litoria_nigrofrenata, Litoria_nasuta, Litoria_tornieri, Litoria_inermis, Litoria_pallida, Litoria_latopalmata, Litoria_bicolor, Litoria_fallax, Litoria_genimaculata, Litoria_andiirrmalin, Litoria_wilcoxii, Litoria_jungguy, Litoria_caerulea, Litoria_gracilenta, Litoria_chloris, Litoria_xanthomera, Cyclorana_brevipes, Cyclorana_novaehollandiae, Cyclorana_manya, Cyclorana_cultripes, Litoria_alboguttata, Cyclorana_longipes, Nyctimystes_dayi, Litoria_nannotis, Litoria_lorica, Litoria_rheocola, Litoria_nyakalensis, Litoria_infrafrenata
#> $x
#> class       : SpatRaster 
#> size        : 90, 68, 3  (nrow, ncol, nlyr)
#> resolution  : 0.1, 0.1  (x, y)
#> extent      : 144.0157, 150.8157, -23.044, -14.044  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#> source      : rast.presab.tif 
#> names       : Litoria_revelata, Litoria_rothii, Litoria_longirostris 
#> min values  :                0,              0,                    0 
#> max values  :                1,              1,                    1 
#> 
#> $tree
#> 
#> Phylogenetic tree with 3 tips and 1 internal node.
#> 
#> Tip labels:
#>   Litoria_revelata, Litoria_rothii, Litoria_longirostris
#> 
#> Unrooted; includes branch length(s).
#> 
#> $edge.path
#>                      [,1] [,2] [,3]
#> Litoria_revelata        1    0    0
#> Litoria_rothii          0    1    0
#> Litoria_longirostris    0    0    1
#> 
#> $branch.length
#> [1] 0.3333337 0.3333331 0.3333331
#> 
#> $branch.length.alt
#> [1] 0.3333333 0.3333333 0.3333333
#> 
#> $n.descendants
#> [1] 1 1 1
#> 
# }
```
