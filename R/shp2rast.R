#' Rasterize shapefile
#'
#' The function will rasterize the shapefile using the parameters of y, a
#' SpatRaster. When the argument y is provided, the resolution parameter is
#' ignored. When the argument ymask is TRUE, y is used as a mask for x.
#'
#' @inheritParams terra::rasterize
#' @inheritParams terra::rast
#' @param ymask logical. If TRUE, y will be used as a mask for x.
#' @param sps.col character. It should be a variable name in x.
#'
#' @return SpatRaster
#'
#'
#' @examples
#' \donttest{
#' library(terra)
#'
#' shp <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp",
#'                               package="phyloraster"))
#'
#' # create a polygon to use as mask with an extent
#' e <- terra::ext(113, 123, -43.64, -33.90)
#' p <- terra::as.polygons(e, crs="")
#' coun.crop <- terra::crop(p, terra::ext(shp))
#' coun.rast <- terra::rasterize(coun.crop,
#' terra::rast(terra::ext(shp), resolution = 0.5))
#'
#' plot(coun.rast, col = "green")
#'
#' # rasterizing with the mask of the polygon
#' shp.t <- shp2rast(shp, y = coun.rast, sps.col = "BINOMIAL",
#' ymask = TRUE, background = 0)
#' plot(shp.t, col = c("grey", "green"))
#'
#' # rasterizing without using mask
#' shp.t2 <- shp2rast(shp, sps.col = "BINOMIAL", ymask = FALSE,
#' background = NA, resolution = 0.1)
#' plot(shp.t2[[9]], col = c("grey", "green"))
#' }
#' @export
#'
shp2rast <- function(x, y = NULL, sps.col, ymask = FALSE, background = NA,
                     touches = TRUE, resolution, values = 1,
                     filename = NULL, ...){

  if(!inherits(x, "SpatVector")){
    stop("The object must be of the class 'SpatVector' from terra package.
         See terra::vect.")
  }

  nm <- unique(data.frame(x)[,sps.col]) # get the species names from shapefile x

  ynull <- is.null(y) # y is null?

  if(ynull){
    exte <- terra::ext(x) # extent
    y <- terra::rast(exte, resolution = resolution) # extent raster
  }

  r_list <- list() # list to store the objects created in a loop for

  for(i in seq_along(nm)){
    r_list[[i]] <- terra::rasterize(x[data.frame(x)[,sps.col] == nm[i],], y,
                                    # field = NULL,
                                    values = values,
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
