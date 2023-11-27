#' Calculate the inverse of range size
#'
#' @description Get range size in square kilometers for all cells that are not
#' NA, the inverse of range size and the inverse of range size multiplied by
#' branch length for multiple species using a raster of presence-absence.
#'
#' @inheritParams geo.phylo.ses
#'
#' @param overwrite logical. If TRUE, filename is overwritten
#' @param ... additional arguments to be passed passed down from a
#' calling function.
#'
#' @return SpatRaster and numeric
#'
#' @author Neander Marcel Heming and Gabriela Alves-Ferreira
#'
#' @examples
#' \donttest{
#' # calculating the inverse of range size
#' x <- terra::rast(system.file("extdata", "rast.presab.tif",
#' package="phyloraster"))
#' inv.range(x[[5]])
#' }
#' @export
inv.range <- function(x,
                      filename = "", overwrite = FALSE, ...){
  cores <- 1
  # Test if there is enough memory in the user's pc
  mi <- .fit.memory(x, 4) ## proc in memory = TRUE means that it fits in
  # the pc's memory, so you wouldn't have to use temporary files

  # temporary files
  temp.ir <- paste0(tempfile(), "ir.tif")  # to store the second raster
  temp1 <- paste0(tempfile(), "xe.tif")  # to store the xe raster

  # cell size
  cellSz <- terra::cellSize(x)

  # spp range size
  rs <- range_size(x, cellSz)
  rs[] <- ifelse(rs==0, Inf, rs)

  inv.R <- terra::ifel(is.na(x), NA, cellSz/rs, filename = filename, ...)

  terra::set.names(inv.R, names(x))

  unlink(temp1)

  return(inv.R)
}
