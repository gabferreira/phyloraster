#' Randomize branch lengths for x times
#'
#' This function randomize branch lengths or another vector for a determined number of aleatorizations.
#' @param x numeric. A numeric vector.
#' @param aleats numeric. A number indicating how many times the calculation should be repeated.
#'
#' @return list
#' @author Gabriela Alves-Ferreira
#' @export
#'
# #' @examples
tip.rand <- function(x, aleats){
  aleats <- aleats # number of null models

  bl.random <- list() # to store the branch length in the loop
  for(i in 1:aleats){
    bl.random[[i]] <- sample(x) # randomize branch lengths
  }
  return(bl.random)
}
