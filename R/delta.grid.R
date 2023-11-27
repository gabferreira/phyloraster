.delta.vec <- function(x) {
  diff <- x[2] - x[1] # diference between metric in time 2 and 1
}

#' Calculates the delta of spatialized diversity metrics (richness,
#' phylogenetic endemism, phylogenetic diversity, weighted endemism,
#' evolutionary distinctiveness) between different time periods
#'
#' @description Calculates the delta of spatialized diversity metrics
#' (richness, phylogenetic endemism, phylogenetic diversity, weighted endemism,
#' evolutionary distinctiveness) between different time periods.
#' @param r1 SpatRaster Spatialized diversity metrics for time 1 (e.g
#' phylogenetic diversity in present). To calculate some diversity metrics for
#' rasters see phyloraster::geo.phylo function.
#' @param r2 SpatRaster Spatialized diversity metrics for time 2 (e.g
#' phylogenetic diversity in future).
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed down from a calling
#' function.
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster
#' with that many cores is created and used.
#' @return SpatRaster
#' @details The two input rasters (r1 and r2) must have the same extent.
#' @export
#'
#' @examples
#' \donttest{
#' # data
#' x <- terra::rast(system.file("extdata", "rast.presab.tif",
#' package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex",
#' package="phyloraster"))
#'
#' # metric SR richness
#' riq.pres <- rast.sr(x)
#' # imagine we lost some species in the future
#' riq.fut <- rast.sr(x[[c(1:15)]])
#' dg <- delta.grid(riq.pres, riq.fut)
#' terra::plot(dg)
#' }
delta.grid <- function(r1, r2, filename = NULL, cores = 1, ...){

  # match extension for r1 and r2?
  if(!terra::ext(r1) == terra::ext(r2)){
   stop("Extents do not match. If you use multiple SpatRaster objects,
        they must have the same resolution and extent. See terra::crop and
        terra::intersect to resolve extension issues.")
  }

  # delta
  rdiff <- terra::app(c(r1, r2), fun = .delta.vec, cores = cores, ...)
  names(rdiff) <- c("Delta")

  if(!is.null(filename)){ # to save the rasters when the output filename is
    # provided
    rdiff <- terra::writeRaster(rdiff, filename)
  }

  return(rdiff)
}
