#' Reorder stack according to tree order and get branch length
#'
#' Reorder stack according to tree order and get branch length
#' @param x SpatRaster. A raster of presence-absence
#' @param tree phylo. A dated tree
#'
#' @return SpatRaster, phylo and numeric vector
#' @export
#'
#' @examples
#' \dontrun{
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
#' phylo.pres(ras, tree)
#' }

phylo.pres <- function(x, tree) {
  spat.names <- as.character(names(x)) # to extract species names in the raster
  tree.phy4 <- phylobase::phylo4(tree) # transform phylo in phylo4 class
  labels <- as.character(phylobase::tipLabels(tree.phy4)) # extracting species names in the tree
  on.tree <- intersect(spat.names, labels) # species of the raster that are in the tree
  subtree <- ape::keep.tip(tree, on.tree) # to make a subset of the tree and keep only the species that are in the raster
  stack.reord <- x[[subtree[["tip.label"]]]] # to reorder the stack according to the tree
  if(!class(stack.reord) == "SpatRaster"){ # class "Raster" in "SpatRaster"
    x <- terra::rast(stack.reord)
  } else {
    x <- stack.reord
  }
  subtree <- phylobase::phylo4(subtree) # phylo in phylo4
  branch.length <- as.numeric(phylobase::edgeLength(subtree,
                                                    1:phylobase::nTips(subtree))) # extract branch lengths
  names(branch.length) <- names(x) # add names
  pp <- list(x = x, branch.length = branch.length, tree = subtree)
  return(pp)
}
