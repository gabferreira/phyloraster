#' Load an example dataset with presence-absence data of 33 tree frogs and a phylogenetic tree for this species
#'
#' This function load a phylogenetic tree and a dataset with presence-absence of 33 Australian tree frogs from Rosauer (2017). Also provides coordinates to each site.
#'
#' @param x a data.frame with coordinates and presence-absence data, and a phylogenetic tree.
#'
#' @return data.frame and phylo
# #' @export
#' @source Rosauer, 2017. Available on: <https://github.com/DanRosauer/phylospatial/tree/master/PhyloEndemism_in_R/Tree%20Frog%20Data>
#'
# #' @examples
#'
load.data.rosauer <- function(x){
  x <- phylogrid::dataR
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
  return(list(presab = x, tree = tree))
}
