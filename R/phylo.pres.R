#' Prepare rasters and phylogenetic tree to run community metrics
#'
#' @description Reorder a stack of rasters of species distribution according to
#' tree order and get branch length for each species to calculate diversity
#' metrics using phyloraster::geo.phylo(). The names must be the same in the
#' phylogenetic tree and in the raster for the same species. For example, if you
#' have the name "Leptodactylus_latrans" in the raster and "Leptodactylus latrans"
#' in the tree, the function will not work. The same goes for uppercase and
#' lowercase letters.
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1)
#' for a set of species.
#' @param tree phylo. A dated tree.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @return Returns a list containing a SpatRaster reordered according to the order
#'  that the species appear in the phylogenetic tree, a subtree containing only
#'  the species that are in the stack of rasters and finally two named numerical
#'  vectors containing the branch length and the number of descendants of each species.
#' @author Neander Marcel Heming and Gabriela Alves Ferreira
#' @examples
#' library(phyloraster)
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
#' phylo.pres(x, tree)
#'
#' @export
phylo.pres <- function(x, tree, ...) {

  if(!inherits(x, "SpatRaster")){ # class "Raster" in "SpatRaster"
    # x <- terra::rast(x)
    # warning("Object 'x' has been converted to an object of class 'SpatRaster'")
    stop("Object 'x' must be of class 'SpatRaster'. See ?terra::rast")
  }
  if(!inherits(tree, "phylo")) {
    stop("The phylogeny 'tree' must be of class 'phylo'")
  }

  ## species name check
  spat.names <- names(x) # to extract species names in the raster
  tip.names <- as.character(phylobase::tipLabels(phylobase::phylo4(tree))) # extracting species names in the tree

  int.tip.spat <- intersect(spat.names, tip.names) # species in common for the raster and the tree
  tip.in.spat <- tip.names %in% spat.names # species of the tree that are in the raster
  spat.in.tip <- spat.names %in% tip.names # species of the raster that are in the tree

  if(identical(int.tip.spat, character(0))){
    stop("The SpatRaster 'x' and the phylogeny 'tree' have no species in common, or the species names do not match between them")
  }
  if (sum(!tip.in.spat)>0){
    warning(paste("Some species in the phylogeny 'tree' are missing from the SpatRaster 'x' and were dropped:", tip.names[!tip.in.spat]))
    tree <- ape::keep.tip(tree, int.tip.spat) # to make a subset of the tree and keep only the species that are in the raster
  }
  if (sum(!spat.in.tip)>0){
    warning(paste("Some species in the phylogeny 'tree' are missing from the SpatRaster 'x' and were dropped:", spat.names[!spat.in.tip]))
  }

  tree <- phylobase::phylo4(tree) # phylo in phylo4

  int.tip.spat[] <- phylobase::tipLabels(tree) # ensure tips names are in proper order

  # Get branch length
  branch.length <- stats::setNames(as.numeric(phylobase::edgeLength(tree, int.tip.spat)), int.tip.spat)

  # Get descendant node numbers
  n.descen <- stats::setNames(as.numeric(phylobase::ancestor(tree, int.tip.spat)), int.tip.spat)

  return(list(x = x[[int.tip.spat]], # reorder the stack according to the tree tips
              tree = tree,
              branch.length = branch.length,
              n.descendants = n.descen))

}
