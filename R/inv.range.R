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
#' @export
#' @examples
#' \dontrun{
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' inv.range(x)
#' }
inv.range <- function(x, LR = F, branch.length, filename = NULL, cores = 1, ...){

  # 2 rasters will be generated in this function, let's see if there is enough memory in the user's pc
  sink(nullfile())    # suppress output
  mi <- terra::mem_info(x, 3)[5] != 0 # proc in memory = T TRUE means that it fits in the pc's memory, so you wouldn't have to use temporary files
  sink()

  # temporary files
  # if(!mi){
    temp <- vector("list", length = 3) # to create a temporary vector with the raster number
    temp[[1]] <- paste0(tempfile(), ".tif")  # to store the second raster
    temp[[2]] <- paste0(tempfile(), ".tif")  # to store the second raster
    temp[[3]] <- paste0(tempfile(), ".tif")  # to store the third raster
  # }

  # calculating area
  rs <- phyloraster::range_size(x)

  # calculating inverse of area and inv area x branch length

  # cell.area <- terra::app(c(terra::cellSize(terra::rast(x[[1]])),x),
  #                    function(x, rs){
  #                      x[1]/(x[-1]*rs)
  #                    }, rs = rs)

  cell.area <- terra::cellSize(terra::rast(x)) # to calculate cell size

  inv.R <- terra::ifel(x == 0, 0, cell.area/(x*rs),
                       filename = temp[[2]], overwrite = T) # calculating the inverse of range size
  names(inv.R) <- names(x)

  if(LR == TRUE){

    # branch.length[] <- branch.length/max(branch.length) # to reescale
    LR <- terra::app(inv.R, function(x, bl){
      x * bl
    }, bl = branch.length, filename = temp[[3]], overwrite = T, cores = cores) # calculating the inverse of range size multiplied by branch length

    return(list(area.size = rs, inv.R = inv.R, LR = LR))

  }

  return(list(area.size = rs, inv.R = inv.R))
}
