#' Calculate the inverse of range size
#'
#' @description Get range size in square kilometers for all cells that are not
#' NA, the inverse of range size and the inverse of range size multiplied by
#' branch length for multiple species using a raster of presence-absence.
#'
#' @inheritParams geo.phylo.ses
#' @inheritParams terra::app
#'
#' @param cores Not implemented yet.
#' @param ... additional arguments to be passed passed down from a calling function.
#'
#' @return SpatRaster and numeric
#'
#' @author Neander Marcel Heming and Gabriela Alves-Ferreira
#'
#' @examples
#' # calculating the inverse of range size
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' inv.range(x)
#'
#' @export
inv.range <- function(x, cores = 1,
                      filename = "", overwrite = FALSE, ...){
  cores = 1
  # Test if there is enough memory in the user's pc
  mi <- .fit.memory(x, 4) ## proc in memory = TRUE means that it fits in the pc's memory, so you wouldn't have to use temporary files

  # temporary files
  temp.ir <- paste0(tempfile(), "ir.tif")  # to store the second raster
  temp1 <- paste0(tempfile(), "xe.tif")  # to store the xe raster

  # cell size
  cellSz <- terra::cellSize(x)

  # spp range size
  rs <- range_size(x, cellSz)

  rs <- terra::app(x,
                   function(x, rs){
                     ifelse(is.na(x), NA, rs)
                   }, rs = rs,
                   filename = ifelse(mi, "", temp1),
                   overwrite = TRUE, cores = cores)

  # calculating inverse of area
  inv.R <- terra::sapp(rs,
                       function(x, cellSz){
                         cellSz/x
                       },
                       cellSz,
                       filename = filename, overwrite = overwrite)

  terra::set.names(inv.R, names(x))

  unlink(temp1)

  return(inv.R)
}
