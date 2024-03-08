#' Calculate Evolutionary distinctiveness for each species
#'
#' @description This function calculates evolutionary distinctiveness according
#' to the fair-proportion index.
#'
#' @inheritParams geo.phylo
#'
#' @author Gabriela Alves-Ferreira and Neander Marcel Heming
#'
#' @references Isaac, N. J., Turvey, S. T., Collen, B.,
#' Waterman, C. and Baillie,
#'  J. E. (2007). Mammals on the EDGE: conservation priorities
#'  based on threat
#'  and phylogeny. PLoS ONE 2, e296.
#'
#' @return vector
#'
species.ed <- function(tree){

  edge.info <- tip.root.path(tree)
  branch.length <- edge.info$edge.length
  n.descen <- colSums(edge.info$H1)

  x <- rep(1, nrow(edge.info$H1))
  branch.ed <- (crossprod(edge.info$H1, x)>0) * (branch.length/n.descen)

  sapply(seq_along(x), function(i, H1, branch.ed){
    sum(H1[i,]*branch.ed)
    }, H1=edge.info$H1, branch.ed=branch.ed)
}
