#' Reorder stack according to tree order and get branch length
#'
#' Reorder stack according to tree order and get branch length.
#' @param pres.stack a raster of presence-absence. It can be an object of class 'raster' or 'SpatRaster'
#' @param tree an object of class 'phylo'
#'
#' @return SpatRaster and numeric vector
#' @export
#'
#' @examples
phylo.pres <- function(pres.stack, tree){
  spat.names <- as.character(names(pres.stack))
  tree.phy4 <- phylobase::phylo4(tree) # to extract tipLabels
  labels <- as.character(phylobase::tipLabels(tree.phy4))
  on.tree <- intersect(spat.names,labels)   # raster species that are on the tree
  subtree <- ape::keep.tip(tree, on.tree)

  stack.reord <- pres.stack[[subtree[["tip.label"]]]] # reorder the species stack according to the species order in the tree
  pres.reord <- terra::rast(stack.reord) # SpatRaster
  #names(pres.reord) == subtree[["tip.label"]] # to check wether names match
  subtree <- phylobase::phylo4(subtree)
  branch.length <- as.numeric(phylobase::edgeLength(subtree, 1:phylobase::nTips(subtree)))
  names(branch.length) <- names(pres.reord)
  pp <- list(pres.reord = pres.reord, branch.length = branch.length)
  return(pp)
}
