#' Get range size, the inverse of range size and the inverse of range size multiplied by branch length for multiple species using a presence-absence SpatRast
#'
#' @param pres.reord a presence-absence SpatRaster with the layers ordered according to the tree order
#' @param branch.length a numeric vector containing the branch length of each specie on the same order of the raster
#' @param filename a character vector containing the path in your personal computer to save the rasters
#' @return SpatRaster and Named num vector
#' @export
#' @examples
inv.range <- function(pres.reord, branch.length, filename = NULL){

  temp <- vector("list", length = 3)
  temp[[1]] <- paste0(tempfile(), ".tif")
  if(!is.null(filename)){
    temp[[2]] <- paste0(filename, "inv.range.tif")
    temp[[3]] <- paste0(filename, "LR.tif")
  } else {
    temp[[2]] <- paste0(tempfile(), "inv.range.tif")
    temp[[3]] <- paste0(tempfile(), "LR.tif")
  }
  # message("Transforming 0 in NA") # para aparecer mensagem enquanto calcula
  {
    # zna <- function(x){
    #   ifelse(x == 0, NA, x) # to transform 0 in NA
    # }
    #   pres <- terra::app(pres.reord, zna, filename = temp[[1]]) # 0 in NA raster
    #   srt <- terra::cellSize(pres) # to obtain cell size
    #   names(srt) <- names(pres.reord)
    #
    #   # the function below added all the presences to know how many pixels the species is present
    #   ## That is, this function determines the size of the species distribution in number of pixels
    #   sz <- function(x, i){ # function to calculate range size
    #     terra::expanse(x[[i]])
    #   }
    #   rs <- sapply(1:terra::nlyr(pres), sz, x = pres) # range size
    #   names(rs) <- names(pres)
  }
  area <- terra::cellSize(terra::rast(pres.reord[[1]]), filename = temp[[1]])

  area.to <- terra::expanse(terra::ifel(any(!is.na(pres.reord)), 1, NA))

  rs <- sapply(1:nlyr(pres.reord),
               function(i, a, Z){
                 az <- zonal(a, Z[[i]], sum)
                 az <- az[az[,1]==1,2]
                 ifelse(length(az)==0, 0, az) # evita voltar erro quando nao tem presenca (1), ou seja, se fosse tudo 0 pra alguma especie
               }, a= area, Z=pres.reord)

  rs[] <- rs/area.to # para reescalonar
  branch.length[] <- branch.length/max(branch.length) # para reescalonar

  message("Calculating the inverse of the range size") # para aparecer mensagem enquanto calcula
  # inv.R <- terra::app(pres.reord, function(x, a){
  #   ifelse(x == 0, 0, 1/(x*a))
  # }, a = rs, filename = temp[[2]])

  inv.R <- terra::ifel(pres.reord == 0, 0, 1/(pres.reord*rs),
                       filename = temp[[2]], overwrite = TRUE)
  message("Calculating the inverse of the range size x branch lengths") # para aparecer mensagem enquanto calcula
  LR <- terra::app(inv.R, function(x, bl){
    x*bl
  }, bl = branch.length, filename = temp[[3]])

  unlink(temp[[1]]) # vai apagar o arquivo que nao usaremos mais (raster de presenca)
  return(list(area.size = rs, inv.R = inv.R, LR = LR))
}
