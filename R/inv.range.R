#' Calculate the inverse of range size
#'
#' @description Get range size, the inverse of range size and the inverse of range size multiplied by branch length for multiple species using a raster of presence-absence.
#' @param x SpatRaster. A presence-absence SpatRaster with the layers ordered according to the tree order.
#' @param branch.length numeric. A vector containing the branch length of each specie.
#' @param filename character. Output filename.
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster with that many cores is created and used.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @return SpatRaster and numeric
#' @author Neander Marcel Heming and Gabriela Alves-Ferreira
#' @export
#' @example
#' \dontrun{
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
#' data <- phylo.pres(ras, tree)
#' inv.range(data$x, data$branch.length)
#' }
inv.range <- function(x, branch.length, filename = NULL, cores = 1, ...){

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
    area <- terra::cellSize(terra::rast(x), filename = temp[[1]]) # to calculate cell size

    # area.to <- terra::expanse(terra::ifel(any(!is.na(x)), 1, NA)) #  to calculate total area

    # The function bellow extracts the range size for each species and stores it in a vector
    rs <- sapply(1:terra::nlyr(x), function(i, a, Z){
      az <- terra::zonal(a, Z[[i]], sum)
      az <- az[az[,1]==1,2]
      ifelse(length(az)==0, 0, az) # avoids returning an error when there is no presence (1), that is, if any species had only 0 in the entire raster
    }, a= area, Z = x)

    # rs[] <- rs/area.to # to reescale
    names(rs) <- names(x) # to add names
    # branch.length[] <- branch.length/max(branch.length) # to reescale
  }

  # calculating inverse of area and inv area x branch length
  {
    inv.R <- terra::ifel(x == 0, 0, 1/(x*rs),
                         filename = temp[[2]], overwrite = T) # calculating the inverse of range size

    LR <- terra::app(inv.R, function(x, bl){
      x*bl
    }, bl = branch.length, filename = temp[[3]], overwrite = T, cores = cores) # calculating the inverse of range size multiplied by branch length
  }

  unlink(temp[[1]]) # delete the files

  return(list(area.size = rs, inv.R = inv.R, LR = LR))
}
