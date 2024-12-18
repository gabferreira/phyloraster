#' Calculate Evolutionary distinctiveness for each species
#'
#' @description This function calculates evolutionary distinctiveness according
#' to the fair-proportion index for each species.
#'
#' @inheritParams geo.phylo
#'
#' @author Neander Marcel Heming and Gabriela Alves-Ferreira
#'
#' @references Isaac, N. J., Turvey, S. T., Collen, B.,
#' Waterman, C. and Baillie, J. E. (2007). Mammals on the EDGE:
#' conservation priorities based on threat and phylogeny. PLoS ONE 2, e296.
#'
#' @return data.frame
#' @examples
#' \donttest{
#' library(phyloraster)
#' tree <- ape::read.tree(system.file("extdata", "tree.nex",
#' package="phyloraster"))
#' plot(tree)
#' ed <- species.ed(tree)
#' ed
#' }
#' @export
species.ed <- function(tree){

  edge.info <- tip.root.path(tree)
  branch.length <- edge.info$edge.length
  n.descen <- colSums(edge.info$H1)
  spps <- tree$tip.label

  x <- rep(1, nrow(edge.info$H1))
  branch.ed <- (crossprod(edge.info$H1, x)>0) * (branch.length/n.descen)

  ed <- vapply(seq_along(x), function(i, H1, branch.ed) {
    sum(H1[i, ] * branch.ed)
  }, numeric(1), H1 = edge.info$H1, branch.ed = branch.ed)

  ed <- data.frame(ED = stats::setNames(ed, spps))
  return(ed)
}
