#' Compute species tip length
#'
#' Computation of species tip length using a phylogeny.
#'
#' @inheritParams phylo.pres
#' @param edge.info Object returned by \code{\link{tip.root.path}} consisting of
#' a list containing the edge matrix (H1) with the path from tip to root and and
#' a numeric vector (edge.length) giving the length of each branch of the tree.
#'
#' @returns returns a numeric vector giving the length of species branch.
#'
#' @details
#' Calculates tip lengths for all species in a phylogeny
#'
#' @author Neander M. Heming
#'
#' @examples
#' library(phyloraster)
#' tree <- ape::read.tree(system.file("extdata", "tree.nex",
#' package="phyloraster"))
#'
#' species.tip.length(tree)
#'
#' library(ape)
#' set.seed(1)
#' tree <- rtree(n=40)
#'
#' plot(tree)
#'
#' species.tip.length(tree)
#'
#' edge.info <- tip.root.path(tree)
#'
#' species.tip.length(edge.info = edge.info)
#'
#' @export
species.tip.length <- function(tree = NULL, edge.info = NULL, ...){

  if(is.null(edge.info)){
    edge.info <- phyloraster::tip.root.path(tree)
  }

  cp <- crossprod(edge.info$H1,
                  matrix(1, nrow = nrow(edge.info$H1))) == 1
  return(stats::setNames(edge.info$edge.length[cp],
                  rownames(edge.info$H1)))
}
