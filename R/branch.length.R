#' Compute species branch lengths
#'
#' Computation of species branch lengths using a phylogeny
#'
#' @inheritParams phylo.pres
#' @param edge.info Object returned by \code{\link{tip.root.path}} consisting of
#' a list containing the edge matrix (H1) with the path from tip to root and and
#' a numeric vector (edge.length) giving the length of each branch of the tree.
#'
#' @returns returns a numeric vector giving the length of species branch.
#'
#' @details
#' Calculates species branch lengths using the edge length and edge path
#'across the phylogeny
#'
#' @author Neander M. Heming
#'
#' @examples
#' library(phyloraster)
#' tree <- ape::read.tree(system.file("extdata", "tree.nex",
#' package="phyloraster"))
#'
#' species.branch.length(tree)
#'
#' library(ape)
#' set.seed(1)
#' tree <- rtree(n=40)
#'
#' plot(tree)
#'
#' species.branch.length(tree)
#'
#' edge.info <- tip.root.path(tree)
#'
#' species.branch.length(edge.info = edge.info)
#'
#' @export
species.branch.length <- function(tree, edge.info = NULL, ...){

  if(is.null(edge.info)){
    edge.info <- phyloraster::tip.root.path(tree)
  }

  return(crossprod(matrix(edge.info$edge.length, nrow = ncol(edge.info$H1)),
                   t(edge.info$H1))[1,])
  # colSums(t(edge.info$H1) * edge.info$edge.length)
}


#' Compute species tip length
#'
#' Computation of species tip length
#'
#' @inheritParams species.branch.length
#'
#' @returns returns a numeric vector giving the length of species branch.
#'
#' @details
#' Calculates tip lengths for all species
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
