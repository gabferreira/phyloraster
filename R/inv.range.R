#' Get range size, the inverse of range size and the inverse of range size multiplied by branch length for multiple species using a presence-absence SpatRast
#'
#' @param pres.reord SpatRaster. A presence-absence SpatRaster with the layers ordered according to the tree order
#' @param branch.length numeric. A vector containing the branch length of each specie on the same order of the raster
#' @param filename character. A vector containing the path in your personal computer to save the rasters
#' @return SpatRaster and numeric
#' @export
#' @examples
#' \dontrun{
#' library(phylogrid)
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
#' data <- phylo.pres(ras, tree)
#' inv.range(data$pres.reord, data$branch.length)
#'}
inv.range <- function(pres.reord, branch.length, filename = NULL){

  # temporary files
  {
    temp <- vector("list", length = 3) # to create a temporary vector with the raster number
    temp[[1]] <- paste0(tempfile(), ".tif")  # to store the first raster
    if(!is.null(filename)){ # to save the rasters directly to your computer when the directory is given
      temp[[2]] <- paste0(filename, "inv.range.tif")
      temp[[3]] <- paste0(filename, "LR.tif")
    } else {
      temp[[2]] <- paste0(tempfile(), "inv.range.tif")
      temp[[3]] <- paste0(tempfile(), "LR.tif")
    }
  }

  # calculating area
  {
    area <- terra::cellSize(terra::rast(pres.reord[[1]]), filename = temp[[1]]) # to calculate cell size

    area.to <- terra::expanse(terra::ifel(any(!is.na(pres.reord)), 1, NA)) #  to calculate area total

    # The function bellow extracts the range size for each species and stores it in a vector
    rs <- sapply(1:terra::nlyr(pres.reord),
                 function(i, a, Z){
                   az <- terra::zonal(a, Z[[i]], sum)
                   az <- az[az[,1]==1,2]
                   ifelse(length(az)==0, 0, az) # avoids returning an error when there is no presence (1), that is, if any species had only 0 in the entire raster
                 }, a= area, Z=pres.reord)

    rs[] <- rs/area.to # to reescale
    names(rs) <- names(pres.reord) # to add names
    branch.length[] <- branch.length/max(branch.length) # to reescale
  }

  # calculating inverse of area and inv area x branch length
  {
    # message("Calculating the inverse of the range size") # show message while calculate

    inv.R <- terra::ifel(pres.reord == 0, 0, 1/(pres.reord*rs),
                         filename = temp[[2]], overwrite = TRUE) # calculating the inverse of range size
    # message("Calculating the inverse of the range size x branch lengths")

    LR <- terra::app(inv.R, function(x, bl){
      x*bl
    }, bl = branch.length, filename = temp[[3]]) # calculating the inverse of range size multiplied by branch length
  }

  unlink(temp[[1]]) # delete the archive

  return(list(area.size = rs, inv.R = inv.R, LR = LR))
}
