#' Rasterize shapefile
#'
#' The function will rasterize the shapefile using the parameters of y, a spatraster. When the argument y is provided, the resolution parameter is ignored. When the argument ymask is TRUE, y is used as a mask for x.
#'
#' @inheritParams terra::rasterize
#' @param ymask SpatVector Mask used to delimit the region of interest, like the shapefile of a country for example
#' @param sps.col character. It should a variable name in x
#' @param resolution numeric. A numeric vector of length 1 or 2 to set the resolution
#' @return SpatRaster
#' @export
#' @examples
#' \dontrun{
#' shp <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp",
#'                               package="phylogrid"))
#' library(rnaturalearth)
#' library(terra)
#' countries <- terra::vect(ne_countries()) # mapa mundi
#' coun.crop <- terra::crop(countries, ext(shp)) # cut by the total extension of the polygons
#' coun.rast <- terra::rasterize(coun.crop,
#'                       terra::rast(ext(shp), resolution = 0.5))
#' plot(coun.rast, col = "green")
#'
#' rasterizing with a mask of a country for example
#' teste <- shp2rast(shp, y = coun.rast, sps.col = "BINOMIAL", ymask = TRUE, background = 0)
#' plot(teste[[1:3]], col = c("grey", "green"))
#'
# rasterizing based on extent and without using mask
#' teste2 <- shp2rast(shp, sps.col = "BINOMIAL", ymask = FALSE, background = NA, resolution = 0.5)
#' plot(teste2[[1:3]], col = c("grey", "green"))
#' }
shp2rast <- function(x, y = NULL, sps.col, ymask = FALSE, background = NA,
                     touches = TRUE, resolution, filename = NULL, ...){

  if(!class(x) == "SpatVector"){
    x <- terra::vect(x)
  }

  nm <- unique(data.frame(x)[,sps.col])

  ynull <- is.null(y) # y is null?

  if(ynull){
    exte <- terra::ext(x) # extent
    y <- terra::rast(exte, resolution = resolution) # extent raster
  }

  r_list <- list() # list to store the objects created in a loop for

  for(i in 1:length(nm)){
    r_list[[i]] <- terra::rasterize(x[data.frame(x)[,sps.col] == nm[i],], y,
                                    field = NULL,
                                    value = 1,
                                    background = background,
                                    touches = touches, filename = "", ...)
  }

  rt <- terra::rast(r_list) # raster stack
  names(rt) <- nm # names

  # applying a mask
  if(!ynull & ymask){
    rt <- terra::mask(rt, y)
  }

  if (!is.null(filename)){ # to save the rasters when the path is provide
    rt <- terra::writeRaster(rt, filename, ...)
  }

  return(rt)
}
