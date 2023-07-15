#' Calculate the inverse of range size
#'
#' @description Get range size in square kilometers for all cells that are not NA, the inverse of range size and the inverse of range size multiplied by branch length for multiple species using a raster of presence-absence.
#' @param x SpatRaster. A presence-absence SpatRaster with the layers ordered according to the tree order.
#' @param LR logical. If LR = TRUE, the function will returns a raster with the clade range.
#' @param branch.length numeric. A vector containing the branch length of each specie.
#' @param filename character. Output filename.
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster with that many cores is created and used.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @return SpatRaster and numeric
#' @author Neander Marcel Heming and Gabriela Alves-Ferreira
#' @examples
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' inv.range(x)
#'
#' @export
inv.range <- function(x, branch.length = NULL, cores = 1, ...){

  # 2 rasters will be generated in this function, let's see if there is enough memory in the user's pc
  mi <- .fit.memory(x, 2) # proc in memory = T TRUE means that it fits in the pc's memory, so you wouldn't have to use temporary files

  # temporary files
  temp <- vector("list", length = 3) # to create a temporary vector with the raster number
  temp[[1]] <- paste0(tempfile(), "blr.tif")  # to store the second raster
  temp[[2]] <- paste0(tempfile(), "ir.tif")  # to store the second raster

  # calculating area
  rs <- phyloraster::range_size(x)

  # calculating inverse of area and inv area x branch length

  inv.R <- terra::app(c(terra::cellSize(x), x),
                      function(x, rs){
                        x[1]/(ifelse(is.na(x[-1]), NA, 1)*rs)
                      }, rs=rs,
                      filename = ifelse(mi, "", temp[[1]]), overwrite = T, cores = cores)

  if(!is.null(branch.length)){

    # branch.length[] <- branch.length/max(branch.length) # to reescale
    LR <- terra::app(inv.R,
                     function(x, bl){
                       x * bl
                     }, bl = branch.length,
                     filename = ifelse(mi, "", temp[[2]]), overwrite = T, cores = cores) # calculating the inverse of range size multiplied by branch length

    return(list(area.size = rs, inv.R = inv.R, LR = LR))

  } else {

    return(list(area.size = rs, inv.R = inv.R))

  }

}
