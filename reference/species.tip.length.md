# Compute species tip length

Computation of species tip length using a phylogeny.

## Usage

``` r
species.tip.length(tree = NULL, edge.info = NULL, ...)
```

## Arguments

- tree:

  phylo. A dated tree.

- edge.info:

  Object returned by
  [`tip.root.path`](https://gabferreira.github.io/phyloraster/reference/tip.root.path.md)
  consisting of a list containing the edge matrix (H1) with the path
  from tip to root and and a numeric vector (edge.length) giving the
  length of each branch of the tree.

- ...:

  additional arguments to be passed passed down from a calling function.

## Value

returns a numeric vector giving the length of species branch.

## Details

Calculates tip lengths for all species in a phylogeny

## Author

Neander M. Heming

## Examples

``` r
library(phyloraster)
tree <- ape::read.tree(system.file("extdata", "tree.nex",
package="phyloraster"))

species.tip.length(tree)
#>          Litoria_revelata            Litoria_rothii      Litoria_longirostris 
#>                  0.589016                  0.589015                  0.395144 
#>          Litoria_dorsalis           Litoria_rubella      Litoria_nigrofrenata 
#>                  0.395144                  0.589015                  0.490836 
#>            Litoria_nasuta          Litoria_tornieri           Litoria_inermis 
#>                  0.288907                  0.196139                  0.129514 
#>           Litoria_pallida       Litoria_latopalmata           Litoria_bicolor 
#>                  0.129514                  0.196139                  0.511002 
#>            Litoria_fallax      Litoria_genimaculata      Litoria_andiirrmalin 
#>                  0.511002                  0.466543                  0.423350 
#>          Litoria_wilcoxii           Litoria_jungguy          Litoria_caerulea 
#>                  0.100568                  0.100568                  0.390538 
#>        Litoria_gracilenta           Litoria_chloris        Litoria_xanthomera 
#>                  0.122460                  0.054075                  0.054075 
#>        Cyclorana_brevipes Cyclorana_novaehollandiae           Cyclorana_manya 
#>                  0.253980                  0.253980                  0.210977 
#>       Cyclorana_cultripes       Litoria_alboguttata        Cyclorana_longipes 
#>                  0.210978                  0.178724                  0.178724 
#>          Nyctimystes_dayi          Litoria_nannotis            Litoria_lorica 
#>                  0.749621                  0.348570                  0.348570 
#>          Litoria_rheocola       Litoria_nyakalensis      Litoria_infrafrenata 
#>                  0.267000                  0.267000                  0.998948 

library(ape)
#> 
#> Attaching package: ‘ape’
#> The following objects are masked from ‘package:terra’:
#> 
#>     rotate, trans, zoom
set.seed(1)
tree <- rtree(n=40)

plot(tree)


species.tip.length(tree)
#>        t25        t15        t33        t20        t35         t6        t10 
#> 0.05893438 0.87626921 0.77891468 0.79730883 0.35319727 0.27026015 0.21320814 
#>        t37        t28        t38        t26        t12        t40        t23 
#> 0.47811803 0.92407447 0.97617069 0.35672691 0.14821156 0.10318424 0.44628435 
#>        t36        t32         t8        t29        t30         t7        t19 
#> 0.64010105 0.17344233 0.75482094 0.20754511 0.59571200 0.57487220 0.07706438 
#>        t34        t22        t14         t2        t13        t16        t18 
#> 0.92861520 0.56090075 0.98509522 0.68278808 0.60154122 0.23886868 0.45257083 
#>        t39         t1         t3        t31        t24        t27         t4 
#> 0.10498764 0.86454495 0.55715954 0.32877732 0.18086636 0.52963060 0.21269952 
#>         t9        t11        t21        t17         t5 
#> 0.28479048 0.44623532 0.77998489 0.41312421 0.06380848 

edge.info <- tip.root.path(tree)

species.tip.length(edge.info = edge.info)
#>        t25        t15        t33        t20        t35         t6        t10 
#> 0.05893438 0.87626921 0.77891468 0.79730883 0.35319727 0.27026015 0.21320814 
#>        t37        t28        t38        t26        t12        t40        t23 
#> 0.47811803 0.92407447 0.97617069 0.35672691 0.14821156 0.10318424 0.44628435 
#>        t36        t32         t8        t29        t30         t7        t19 
#> 0.64010105 0.17344233 0.75482094 0.20754511 0.59571200 0.57487220 0.07706438 
#>        t34        t22        t14         t2        t13        t16        t18 
#> 0.92861520 0.56090075 0.98509522 0.68278808 0.60154122 0.23886868 0.45257083 
#>        t39         t1         t3        t31        t24        t27         t4 
#> 0.10498764 0.86454495 0.55715954 0.32877732 0.18086636 0.52963060 0.21269952 
#>         t9        t11        t21        t17         t5 
#> 0.28479048 0.44623532 0.77998489 0.41312421 0.06380848 
```
