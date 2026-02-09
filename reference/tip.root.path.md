# Compute tree edge lengths and node paths from root to each tip

Computation of tree edge lengths and node paths from root to each tip to
calculate PD for a entire phylogeny (= sum of all edge or branch
lengths)

## Usage

``` r
tip.root.path(tree)
```

## Arguments

- tree:

  phylo. A dated tree.

## Value

returns a list with two components: matrix H1 representing the paths
through the tree from root to each tip, and edge.length a numeric vector
giving the length of each branch in the tree. Some matrix algebra and a
summation of the resulting vector gives the whole-tree PD value.

## Details

Based on the algorithm FastXtreePhylo of Peter D. Wilson

## Author

Peter Wilson

## Examples

``` r
library(phyloraster)
tree <- ape::read.tree(system.file("extdata", "tree.nex",
package="phyloraster"))

fxtp <- tip.root.path(tree)
H1 <- fxtp$H1
edge.length <- fxtp$edge.length
# PD for the whole community
pres <- rep(1, nrow(H1))
sum((crossprod(H1, pres)>0) * edge.length)
#> [1] 14.86832

# PD for a random subset of the community
pres <- sample(c(1, 0), nrow(H1), TRUE)
sum((crossprod(H1, pres)>0) * edge.length)
#> [1] 9.64867
```
