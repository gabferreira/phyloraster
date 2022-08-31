#' Transform data.frame to raster
#'
#' @param x data.frame. A data.frame or matrix with species names in columns and sites in rows. The first two columns must provide x and y coordinates.
#'
#' @return SpatRaster
#' @export
#'
#' @examples
#' \dontrun{
#' dat <- phylogrid::dataR
#' df.to.rast(dat)
#' }
df2rast <- function(x){
  if(!class(x) == "matrix"){
    xm <- as.matrix(x) # to matrix
  }
  terra::rast(xm, type = "xyz") # transforming in raster
}
