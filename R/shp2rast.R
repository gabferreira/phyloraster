#' Rasterize shapefile
#'
#' The function transform a shapefile to a raster with the same extent.
#'
#' @param shp SpatialPolygonsDataFrame or a SpatVector. A shapefile representing species distribution
#' @param sps vector. A vector of characters with unique (not duplicated) species names presents in the shapefile
#' @param resolution numeric. A numeric vector of length 1 or 2 to set the resolution
#' @return SpatRaster
#' @export
#' @examples
#' \dontrun{
#' shp <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp", package="phylogrid"))
#' spps <- shp$BINOMIAL
#' shp2rast(dat$IUCN_shapefile, spps, resolution = 0.1)
#' }
#'
shp2rast <- function(shp, sps, resolution){

  if(!class(shp) == "SpatVector"){

    shp <- terra::vect(shp)
    exte <- terra::ext(shp) # extent
    rr <- terra::rast(exte, resolution = resolution) # extent raster

    r_list <- list() # list to store the objects created in a loop for

    for(i in 1:length(sps)){
      r_list[[i]] <- terra::rasterize(shp[i,], rr)
      }

    rt <- terra::rast(r_list) # raster stack
    names(rt) <- sps # names

    return(rt)
  } else {

    exte <- terra::ext(shp) # extent
    rr <- terra::rast(exte, resolution = resolution) # extent raster

    r_list <- list() # list to store the objects created in a loop for

    for(i in 1:length(sps)){
      r_list[[i]] <- terra::rasterize(shp[i,], rr)
    }

    rt <- terra::rast(r_list) # raster stack
    names(rt) <- sps # names

    return(rt)
  }
}
