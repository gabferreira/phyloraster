#' Transform data.frame to raster
#'
#' @param x data.frame. A data.frame or matrix with species names in columns and sites in rows. The first two columns must provide x and y coordinates.
#' @inheritParams terra
#' @param CRS character. A coordinate reference system (projection) for a Raster object.
#' @return SpatRaster
#' @export
#'
#' @examples
#' \dontrun{
#' dat <- phylogrid::load.data.rosauer()
#' df2rast(dat$presab, CRS = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
#' }
#'
#'
df2rast <- function(x, CRS){
  if(!class(x) == "matrix"){
    xm <- as.matrix(x) # to matrix
  }
    r <- terra::rast(xm, type = "xyz") # transforming in raster
    terra::crs(r) <- CRS
    return(r)
}

df2rast <- function(x, resolution = NULL){
  if(!class(x) == "matrix"){
    xm <- as.matrix(x) # to matrix
  }
  if(!is.null(resolution))
  {
  r <- terra::rast(xm, type = "xyz") # transforming in raster
  res(r) <- resolution
  return(r)
  } else {
  r <- terra::rast(xm, type = "xyz") # transforming in raster
  }
}
