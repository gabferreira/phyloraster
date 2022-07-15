
#' Phylogenetic diversity - Null model
#'
#' @param pres.reord SpatRaster. A presence-absence SpatRaster with the layers ordered according to the tree order
#' @param branch.length numeric. A Named numeric vector containing the branch length of each specie
#' @param aleats numeric. A vector containing the number of aleatorizations
#' @return SpatRaster
#' @export
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological conservation, 61(1), 1-10.
# #' @examples
pd.null.mod <- function(pres.reord, branch.length, aleats){
  {
    aleats <- 10
    temp <- vector("list", length = aleats) # to create a temporary vector with the raster number
    rast.pd <- list() # to store the rasters
  }

  ## null model (bootstrap structure)
  for(i in 1:aleats){
    bl.random <- sample(branch.length, replace = T) # aleatorize the branch lenght
    ## check if the values are differents
    # branch.length == bl.random
    for(j in 1:aleats){
      temp[[j]] <- paste0(tempfile(), j, ".tif") # directory to store the rasters
    }
    rast.pd[[i]] <- terra::app(pres.reord, fun = .vec.pd,
                               branch.length = bl.random, filename = temp[[i]])
  }

  rast.pd <- terra::rast(rast.pd) # to transform a list in raster

  {
    temp2 <- vector("list", length = 2)
    temp2[[1]] <- paste0(tempfile(), "pd.mean.bootstraped.tif")
    temp2[[2]] <- paste0(tempfile(), "pd.sd.bootstraped.tif")

    rast.pd.mean <- terra::app(rast.pd, fun = mean, na.rm = TRUE, filename = temp2[[1]]) # mean
    rast.pd.sd <- terra::app(rast.pd, fun = , na.rm = TRUE, filename = temp2[[2]]) # sd
  }

  unlink(temp) # delete the archive

  return(c(rast.pd.mean, rast.pd.sd))
}
