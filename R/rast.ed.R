#' Calculate Evolutionary distinctiveness for a vector
#'
#' This function calculates evolutionary distinctiveness according to the fair-proportion index (Isaac et al., 2007).
#'
#' @param x numeric. A Named numeric vector of presence-absence
#' @param branch.length branch.length numeric. A Named numeric vector of branch length for each specie
#' @param n.descen numeric. A Named numeric vector of number of descendants for each branch
#' @references Isaac, N. J., Turvey, S. T., Collen, B., Waterman, C. & Baillie, J. E. (2007). Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE 2, e296.
#' @return numeric
#'
#' @return
#' @export
#'
#' @examples
.evol.distin <- function(x, branch.length, n.descen){

  x[is.na(x)] <- 0

  if(sum(x) == 0) {
    return(c(NA,NA))
  }

  if(sum(x) != 0){
    pres <- x == 1
    species <- names(x)
    ed <- sum(branch.length[pres]/n.descen[pres])
    ed1 <- ed # tive que adicionar pq a funcao nao funcionava se eu retornasse apenas um raster
  }
  return(c(ed = ed, ed1 = ed1))
}


#' Calculate Evolutionary distinctiveness for a raster
#'
#' This function calculates evolutionary distinctiveness according to the fair-proportion index (Isaac et al., 2007).
#'
#' @inheritParams rast.pd
#' @param x SpatRaster. A presence-absence SpatRaster with the layers ordered according to the tree order.
#' @param branch.length numeric. A Named numeric vector containing the branch length of each specie.
#' @param n.descen numeric. A Named numeric vector of number of descendants for each branch
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed for fun.
#' @references Isaac, N. J., Turvey, S. T., Collen, B., Waterman, C. & Baillie, J. E. (2007). Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE 2, e296.
#'
#' @return SpatRaster
#' @export
#'
#' @examples
#' \dontrun{
#' # raster
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' # phylogenetic tree
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
#' data <- phylogrid::phylo.pres(x, tree)
#' ed <- phylogrid::rast.ed(data$x, data$branch.length, data$n.descen, cores = 2)
#' }
#'

rast.ed <- function(x, branch.length, n.descen, filename = NULL, cores = 1, ...){

  if(!all.equal(names(x), names(branch.length))){
    stop("Species names are not in the same order on 'x' and 'branch.length' arguments! See 'phylogrid::phylo.pres' function.")
  }
    # phylogenetic diversity and richness-relative phylogenetic diversity
    if(cores>1){
      red <- terra::app(x, fun = .evol.distin,
                        branch.length, n.descen, cores = 1)
      red <- red[[1]] # separando so pro primeiro raster devido o problema da funcao app de nao calcular quando eh pra retornar apenas um raster
      names(red) <- c("ED")

    } else {
      red <- terra::app(x, fun = .evol.distin,
                        branch.length, n.descen)
      red <- red[[1]] # separando so pro primeiro raster devido o problema da funcao app de nao calcular quando eh pra retornar apenas um raster
      names(red) <- "ED"
    }

  if (!is.null(filename)){ # to save the rasters when the output filename is provide
    red <- terra::writeRaster(red, filename, ...)
  }

  return(red)
}
