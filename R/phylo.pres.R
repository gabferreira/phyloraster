#' Reorder stack according to tree order and get branch length
#'
#' Reorder stack according to tree order and get branch length.
#' @param pres.stack a raster of presence-absence of class 'SpatRaster'
#' @param tree an object of class 'phylo'
#'
#' @return SpatRaster, phylo and numeric vector
#' @export
#'
#' @examples
phylo.pres <- function(pres.stack, tree) {
  spat.names <- as.character(names(pres.stack))
  tree.phy4 <- phylobase::phylo4(tree)
  labels <- as.character(phylobase::tipLabels(tree.phy4))
  on.tree <- intersect(spat.names, labels)
  subtree <- ape::keep.tip(tree, on.tree)
  stack.reord <- pres.stack[[subtree[["tip.label"]]]]
  if(!class(stack.reord) == "SpatRaster"){
    pres.reord <- terra::rast(stack.reord)
  } else {
    pres.reord <- stack.reord
  }
  subtree <- phylobase::phylo4(subtree)
  branch.length <- as.numeric(phylobase::edgeLength(subtree,
                                                    1:phylobase::nTips(subtree)))
  names(branch.length) <- names(pres.reord)
  pp <- list(pres.reord = pres.reord, branch.length = branch.length, tree = subtree)
  return(pp)
}
