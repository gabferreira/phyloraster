#'  Compute tree edge lengths and node paths from root to each tip
#'
#' Computation of tree edge lengths and node paths from root to each tip to
#' calculate PD for a entire phylogeny (= sum of all edge or branch lengths)
#'
#' @inheritParams phylo.pres
#'
#' @returns returns a list with two components: matrix H1
#' representing the paths through the tree from root to each tip,
#' and edge.length a
#' numeric vector giving the length of each branch in the tree. Some matrix
#' algebra and a summation of the resulting vector gives the whole-tree
#' PD value.
#'
#' @details
#' Based on the algorithm FastXtreePhylo of Peter D. Wilson
#'
#' @author Peter Wilson
#'
#' @examples
#' library(phyloraster)
#' tree <- ape::read.tree(system.file("extdata", "tree.nex",
#' package="phyloraster"))
#'
#' fxtp <- tip.root.path(tree)
#' H1 <- fxtp$H1
#' edge.length <- fxtp$edge.length
#' # PD for the whole community
#' pres <- rep(1, nrow(H1))
#' sum((crossprod(H1, pres)>0) * edge.length)
#'
#' # PD for a random subset of the community
#' pres <- sample(c(1, 0), nrow(H1), TRUE)
#' sum((crossprod(H1, pres)>0) * edge.length)
#'
#' @export
tip.root.path <- function(tree){
  # In edge matrix in a phylo object, OTUs are indexed 1:Ntip, the root node has
  # index Ntip + 1, and internal nodes are indexed from root.node to
  # total.nodes. This is clever because it lets us very efficiently
  # cross-reference between edge indicies and the left and right terminal node
  # indices of each edge. So, compute the key values:
  root.node <- ape::Ntip(tree) + 1
  total.nodes <- max(tree$edge)

  # Matrix H1 of Schumacher's Xtree method, which is exactly the same as
  # David Nipperess's phylomatrix:
  H1 <- matrix(0, ape::Ntip(tree), ape::Nedge(tree),
               dimnames=list(tree$tip.label, NULL))

  # A short-cut: we can instantly populate H1 with terminal edges:
  rc.ind <- cbind("row"=1:ape::Ntip(tree),
                  "col"=which(tree$edge[,2] < root.node))
  H1[rc.ind] <- 1

  # Make a vector of internal node indices in descending order
  # so we can traverse
  # tree from first nodes below tip nodes down to the root:
  internal.nodes <- seq(total.nodes,root.node,-1)

  child.edges <- next.edge <- logical(length(tree$edge[,2]))

  # Now, visit each internal node accumulating subtended child edge records:
  for(thisNode in internal.nodes){
    # Find the edge which has thisNode on its left:
    next.edge[] <- tree$edge[,2]==thisNode

    # Which edges are children of the node to the right:
    child.edges[] <- tree$edge[,1]==thisNode

    # Do the magic rowSums() trick to accumulate edges subtended by the
    # current edge:
    H1[,next.edge] <- rowSums(H1[,child.edges])
  }

  # All done...package results
  return(list("edge.length" = tree$edge.length,
              "H1" = H1))
  # return(H1)
}


#' Prepare rasters and phylogenetic tree to run community metrics
#'
#' @description Reorder a stack of rasters of species distribution to
#' match the order
#' of the tips of the tree, and get branch length and number of descendants for
#' each species to calculate diversity metrics using phyloraster::geo.phylo().
#' The
#'  branch length and the number of descendants can be calculated based on the
#'  full
#'  tree or the raster based tree subset.
#'  The names must be the same in the phylogenetic tree and in the raster for
#'  the
#'  same species. For example, if you have the name "Leptodactylus_latrans" in
#'  the raster and "Leptodactylus latrans" in the tree, the function will not
#'  work. The same goes for uppercase and lowercase letters.
#'
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1)
#' for a set of species.
#' @param tree phylo. A dated tree.
#' @param full_tree_metr logical. Whether edge.path, branch length and number
#' of descendants should be calculated with the full (TRUE) or the prunned tree
#' (FALSE).
#' @param ... additional arguments to be passed passed down from a calling
#' function.
#' @return Returns a list containing a SpatRaster reordered according to the
#' order
#'  that the species appear in the phylogenetic tree, a subtree containing only
#'  the species that are in the stack of rasters and finally two named numerical
#'  vectors containing the branch length and the number of descendants of
#'  each species.
#'
#' @author Neander Marcel Heming and Gabriela Alves Ferreira
#'
#' @examples
#' \donttest{
#' library(phyloraster)
#' x <- terra::rast(system.file("extdata", "rast.presab.tif",
#' package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex",
#' package="phyloraster"))
#' phylo.pres(x[[1:3]], tree, full_tree_metr = TRUE)
#'
#' # using the prunned tree
#' phylo.pres(x[[1:3]], tree, full_tree_metr = FALSE)
#' }
#' @export
phylo.pres <- function(x, tree, full_tree_metr = FALSE, ...) {

  if(!inherits(x, "SpatRaster")){ # class "Raster" in "SpatRaster"
    stop("Object 'x' must be of class 'SpatRaster'. See ?terra::rast")
  }

  if(!inherits(tree, "phylo")) {
    stop("The phylogeny 'tree' must be of class 'phylo'")
  }

  ## species name check
  spat.names <- names(x) # to extract species names in the raster
  tip.names <- tree$tip.label # as.character(phylobase::tipLabels
  # (phylobase::phylo4(tree))) # extracting species names in the tree

  int.tip.spat <- intersect(tip.names, spat.names) # species in common for the
  # raster and the tree
  tip.in.spat <- tip.names %in% spat.names # species of the tree that are in
  # the raster
  spat.in.tip <- spat.names %in% tip.names # species of the raster that are in
  # the tree

  if(identical(int.tip.spat, character(0))){
    stop("The SpatRaster 'x' and the phylogeny 'tree' have no species in common,
         or the species names do not match between them")
  }
  if (sum(!spat.in.tip)>0){
    warning(paste("Some species in the phylogeny 'tree' are missing from the
                  SpatRaster 'x' and were dropped:",
                  paste0(spat.names[!spat.in.tip], collapse = ", ")))
  }
  if (sum(!tip.in.spat)>0){
    warning(paste("Some species in the phylogeny 'tree' are missing from the
                  SpatRaster 'x' and were dropped:",
                  paste0(tip.names[!tip.in.spat], collapse = ", ")))

  }
  subtree <- ape::keep.tip(tree, int.tip.spat) # to make a subset of the tree
  # and keep only the species that are in the raster

  # int.tip.spat[] <- subtree$tip.label #phylobase::tipLabels(phylobase::
  #phylo4(subtree)) # ensure tips names are in proper order

  if(full_tree_metr){
    ## Compute node paths through the tree from root to each tip
    edge.info <- tip.root.path(tree)
    # Get descendant node numbers
    n.descen <- stats::na.exclude(
      as.numeric(phylobase::ancestor(phylobase::phylo4(tree))))

    ### sum common edges to reduce matrix dim
    # rownames(edge.info$H1) <- tree$tip.label
    H1agg <- stats::aggregate(stats::reformulate(int.tip.spat,
                                                 'edge.length'), FUN="sum",
                              data = cbind(t(edge.info$H1[int.tip.spat,]),
                                           edge.length = edge.info$edge.length))

    edge.info$H1 <- t(H1agg)[int.tip.spat,]
    edge.info$edge.length <- H1agg[,"edge.length"]

  } else {
    ## Compute node paths through the tree from root to each tip
    edge.info <- tip.root.path(subtree)
    # Get descendant node numbers
    n.descen <- stats::na.exclude(as.numeric(phylobase::ancestor(
      phylobase::phylo4(subtree))))

  }

  return(list(x = x[[int.tip.spat]], # reorder the stack according to the
              # tree tips
              tree = subtree,
              edge.path = edge.info$H1[int.tip.spat,],
              branch.length = edge.info$edge.length,
              n.descendants = n.descen))

}
