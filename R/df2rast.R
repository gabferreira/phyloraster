#' Transform a data.frame to raster
#'
#' @description The function transform a data.frame or a matrix of presence-absence in a raster of distribution.
#' @param x data.frame. A data.frame or matrix with species names in columns and sites in rows. The first two columns must provide longitude and latitude, respectively.
#' @param CRS character. A coordinate reference system (projection) for a Raster object.
#' @return SpatRaster
#' @export
#' @examples
#' \dontrun{
#' dat <- phylogrid::load.data.rosauer()
#' df2rast(dat$presab, resolution = 0.1, CRS = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
#' }

df2rast <- function(x, CRS){
  if(!class(x) == "matrix"){
    xm <- as.matrix(x) # to matrix
  }
    r <- terra::rast(xm, type = "xyz") # transforming in raster
    terra::crs(r) <- CRS
    return(r)
}
