#' Get range size, the inverse of range size and the inverse of range size multiplied by branch length for multiple species using a presence-absence SpatRast
#'
#' @param pres.reord a presence-absence SpatRaster with the layers ordered according to the tree order
#' @param branch.length a numeric vector containing the branch length of each specie on the same order of the raster
#' @param filename a character vector containing the path in your personal computer to save the rasters
#' @return SpatRaster and Named num vector
#' @export
#' @examples
inv.range <- function(pres.reord, branch.length, filename = NULL){
  zna <- function(x){
    ifelse(x == 0, NA, x) # to transform 0 in NA
  }
  temp <- vector("list", length = 3)
  temp[[1]] <- paste0(tempfile(), ".tif")
  if(!is.null(filename)){
    temp[[2]] <- paste0(filename, "inv.range.tif")
    temp[[3]] <- paste0(filename, "LR.tif")
  } else {
    temp[[2]] <- paste0(tempfile(), "inv.range.tif")
    temp[[3]] <- paste0(tempfile(), "LR.tif")
  }
  message("Transforming 0 in NA") # para aparecer mensagem enquanto calcula
  pres <- terra::app(pres.reord, zna, filename = temp[[1]]) # 0 in NA raster
  srt <- terra::cellSize(pres) # to obtain cell size
  names(srt) <- names(pres.reord)

  # the function below added all the presences to know how many pixels the species is present
  ## That is, this function determines the size of the species distribution in number of pixels
  sz <- function(x, i){ # function to calculate range size
    terra::expanse(x[[i]])
  }
  rs <- sapply(1:terra::nlyr(pres), sz, x = pres) # range size
  names(rs) <- names(pres)

  # ai <- function(x, a){
  #   x/a # to calculate the inverse of area size
  # }
  # inv.R <- terra::app(srt, ai, a = rs) # inverse of area size
  message("Calculating the inverse of the range size") # para aparecer mensagem enquanto calcula
  inv.R <- terra::app(srt, function(x, a){
    x/a
  }, a = rs, filename = temp[[2]])

  message("Calculating the inverse of the range size x branch lengths") # para aparecer mensagem enquanto calcula
  LR <- terra::app(inv.R, function(x, bl){
    x*bl
  }, bl = branch.length, filename = temp[[3]])

  unlink(temp[[1]]) # vai apagar o arquivo que nao usaremos mais (raster de presenca)
  resu <- list(area.size = rs, inv.R = inv.R, LR = LR)
  return(resu)
}

