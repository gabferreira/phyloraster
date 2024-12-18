#' Transform a data.frame to raster
#'
#' @description The function transforms a data.frame or a matrix of presence-
#' absence in a raster of distribution.
#' @param x data.frame. A data.frame or matrix with species names in columns
#' and sites in rows. The first two columns must provide longitude and latitude,
#' respectively.
#' @param CRS character. Description of the Coordinate Reference System
#' (map projection) in PROJ.4.
#' @param ... additional arguments to be passed passed down from a calling
#' function.
#' @return SpatRaster
#' @export
#' @examples
#' \donttest{
#' dat <- phyloraster::load.data.rosauer()
#' df2rast(dat$presab, crs = "+proj=longlat +datum=WGS84 +ellps=WGS84
#' +towgs84=0,0,0")
#' }
df2rast <- function(x, CRS = "+proj=longlat +datum=WGS84", ...){

  if(!inherits(x, "matrix")){
    xm <- as.matrix(x) # to matrix
  }

  r <- terra::rast(xm, type = "xyz") # transforming in raster

  if(!CRS == "+proj=longlat +datum=WGS84") {
    terra::crs(r) <- CRS
  } else {
    terra::crs(r) <- "+proj=longlat +datum=WGS84"
  }

  return(r)
}
