.delta.vec <- function(x) {
  diff <- x[2] - x[1] # diference between metric in time 2 and 1
}

#' Calculates the delta of spatialized diversity metrics (richness, phylogenetic endemism, phylogenetic diversity, weighted endemism, evolutionary distinctiveness) between different time periods
#'
#' @description Calculates the delta of spatialized diversity metrics (richness, phylogenetic endemism, phylogenetic diversity, weighted endemism, evolutionary distinctiveness) between different time periods
#' @inheritParams geo.phylo
#' @param r1 SpatRaster Spatialized diversity metrics for time 1 (e.g phylogenetic diversity in present). To calculate some diversity metrics for rasters see phylogrid::geo.phylo function.
#' @param r2 SpatRaster Spatialized diversity metrics for time 2 (e.g phylogenetic diversity in future).
#'
#' @return SpatRaster
#' @export
#'
#' @examples
#' \dontrun{
#' # data
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
#'
#' # metric SE richness
#' riq.pres <- rast.se(x)
#' riq.fut <- rast.se(x[[c(1:15)]]) # imagine we lost some species in the future
#' dg <- delta.grid(riq.pres, riq.fut)
#' plot(dg)
#' }
delta.grid <- function(r1, r2, filename = NULL, cores = 1, ...){

  # delta
  rdiff <- terra::app(c(r1, r2), fun = .delta.vec, cores = cores, ...)
  names(rdiff) <- c("Delta")

  if (!is.null(filename)){ # to save the rasters when the output filename is provide
    rdiff <- terra::writeRaster(rdiff, filename)
  }

  return(rdiff)
}
