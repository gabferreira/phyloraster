#' Randomize branch lengths for n times
#'
#' @description This function randomize branch lengths or another vector for n times.
#' @param x numeric. A numeric vector containing the branch lengths.
#' @param aleats positive integer. A positive integer indicating how many times the calculation should be repeated.
#' @return list of numeric vectors
#' @author Gabriela Alves-Ferreira
#' @export
#' @examples
#' \dontrun{
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
#' data <- phylo.pres(x = x, tree = tree)
#' sa <- tip.rand(data$branch.length, aleats = 10)
#' }
#'
tip.rand <- function(x, aleats){
  aleats <- aleats # number of null models

  bl.random <- list() # to store the branch length in the loop
  for(i in 1:aleats){
    bl.random[[i]] <- sample(x) # randomize branch lengths
  }
  return(bl.random)
}
