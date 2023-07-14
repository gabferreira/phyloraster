#' Randomize branch lengths for n times
#'
#' @description This function randomize branch lengths or another vector for n times.
#' @param x numeric. A numeric vector containing the branch lengths.
#' @param aleats positive integer. A positive integer indicating how many times the calculation should be repeated.
#' @inheritParams base::sample
#'
#' @return list of numeric vectors
#'
#' @author Gabriela Alves-Ferreira
#'
#' @examples
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
#' data <- phylo.pres(x = x, tree = tree)
#' sa <- tip.rand(data$branch.length, aleats = 10)
#'
#'
#' @export
tip.rand <- function(x, aleats, replace = FALSE){
  # aleats <- aleats # number of null models

  bl.random <- list() # to store the branch length in the loop
  for(i in 1:aleats){
    bl.random[[i]] <- sample(x, replace = replace) # randomize branch lengths
  }
  return(bl.random)
}
