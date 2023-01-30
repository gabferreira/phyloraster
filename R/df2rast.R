#' Transform a data.frame to raster
#'
#' @description The function transform a data.frame or a matrix of presence-absence in a raster of distribution.
#' @param x data.frame. A data.frame or matrix with species names in columns and sites in rows. The first two columns must provide longitude and latitude, respectively.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @return SpatRaster
#' @export
#' @examples
#' \dontrun{
#' dat <- phylogrid::load.data.rosauer()
#' df2rast(dat$presab)
#' }

df2rast <- function(x, ...){

  if(!class(x) == "matrix"){
    xm <- as.matrix(x) # to matrix
  }
    r <- terra::rast(xm, type = "xyz") # transforming in raster
    return(r)
}
