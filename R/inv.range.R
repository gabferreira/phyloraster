#' Get range size, the inverse of range size and the inverse of range size multiplied by branch length for multiple species using a presence-absence SpatRast
#'
#' @param pres.reord a presence-absence SpatRaster with the layers ordered according to the tree order
#' @param branch.length a numeric vector containing the branch length of each specie on the same order of the raster
#' @return SpatRaster and Named num vector
#' @export
#' @examples
inv.range <- function(pres.reord, branch.length){
  zna <- function(x){
    ifelse(x == 0, NA, x) # to transform 0 in NA
  }
  pres <- terra::app(pres.reord, zna) # 0 in NA raster

  srt <- terra::cellSize(pres) # to obtain cell size
  names(srt) <- names(pres.reord)

  # the function below added all the presences to know how many pixels the species is present
  ## That is, this function determines the size of the species distribution in number of pixels
  sz <- function(x, i){ # function to calculate range size
    terra::expanse(x[[i]])
  }
  rs <- sapply(1:terra::nlyr(srt), sz, x = pres) # range size
  names(rs) <- names(srt)

  ai <- function(x, a){
    x/a # to calculate the inverse of area size
  }
  area.inv <- terra::app(srt, ai, a = rs) # inverse of area size

  area.tips <- terra::app(area.inv, function(x, bl){
    x * bl
  }, bl = branch.length)

  resu <- list(area.size = rs, area.size.r = srt, area.inv = area.inv, area.tips = area.tips)
  return(resu)
}
