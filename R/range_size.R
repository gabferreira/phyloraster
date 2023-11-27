#' Calculate range size for a set of species using a raster as input
#'
#' @description This function calculate range size in square meters (by default)
#' for all cells that are not NA. The size of the cells is constant in degrees
#' but not in square meters, which was considered in the method applied to
#' calculate the area.
#'
#' @inheritParams geo.phylo.ses
#' @param cellSz SpatRaster. A SpatRaster containing cellSize values.
#' See \code{\link[terra]{cellSize}}
#' @inheritParams terra::cellSize
#' @param ... additional arguments to be passed passed down from a calling
#' function.
#'
#' @author Gabriela Alves Ferreira and Neander Marcel Heming
#'
#' @return vector
#'
#' @examples
#' \donttest{
#' x <- terra::rast(system.file("extdata", "rast.presab.tif",
#' package="phyloraster"))
#' range_size(x[[1:2]], cellSz <- terra::cellSize(x))
#'}
#' @export
range_size <- function(x, cellSz, unit = "m", ...){

  # colocar ifelse do mi
  temp <- vector("list", length = 2) # to create a temporary vector with
  # the raster number
  temp[[1]] <- paste0(tempfile(), ".tif")  # to store the first raster

  miss4 <- arg.check(match.call(), c("cellSz"))
  if(any(miss4)){
    cellSz <- terra::cellSize(terra::rast(x[[1]]), filename = temp[[1]],
                              unit = unit) # to calculate cell size
  }

  # The function bellow extracts the range size for each species and stores
  # it in a vector
  rs <- terra::global(x*cellSz, fun = "sum", na.rm = TRUE)[,1]

  names(rs) <- names(x) # to add names

    unlink(temp) # delete the archive that will not be used anymore
    return(rs)
}
