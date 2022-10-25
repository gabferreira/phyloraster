#' Intern function to resampling a vector according to the observed frequency
#'
#' @param x numeric. A vector containing values to resampling.
#' @param fr data.frame A data.frame with 3 columns (layer, value, count).
#'
#' @return vector
# #' @export
#'
#' @examples
.lyr.sample <- function(x, fr){
  sapply(x, function(x, fr){
    if(is.na(x)){
      return(NA)
    } else {
      return(sample(fr$value, 1, prob = fr$count))
    }
  }, fr = fr)
}

#' Randomize a set of rasters according to the observed frequency.
#'
#' Randomize a set of rasters according to the observed frequency using the sites (by cells), species (by layer) or both (layers and cells).
#'
#' @param x SpatRaster. A presence-absence SpatRaster.
#' @param aleats numeric. A number indicating how many times the calculation should be repeated.
#' @param random character. Character indicating the type of randomization to be used. The available types are by "site", "specie" or "fullspat".
#' @return SpatRaster
#' @author Neander Marcel Heming and Gabriela Alves-Ferreira
#' @export
#'
#' @examples
#' \dontrun{
#' aleats = 10
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' sr <- spat.rand(ras, aleats = 10, random = "site")
#' plot(sr)
#' }
spat.rand <- function(x, aleats, random = c("site", "specie", "fullspat")){

  aleats <- aleats # number of null models
  fr <- terra::freq(x)

  if (random == "site"){

    resu <- list()

    for(i in 1:aleats){
      # randomize by layers- sites
      resu[[i]] <- terra::rast(lapply(1:terra::nlyr(x),
                                      function(i, r, fr){
                                        terra::app(r[[i]], fun=.lyr.sample, fr=fr[fr$layer==i,])
                                      }, r = x, fr = fr))
    }

    resu <- terra::rast(resu) # to transform a list in raster

  } else if (random == "specie") {

    resu <- list()

    for(i in 1:aleats){
      ### randomize by cells- species in each site
      resu[[i]] <- terra::app(x, sample)
    }

    resu <- terra::rast(resu) # to transform a list in raster

  } else if (random == "fullspat") {

    resu <- list()

    for(i in 1:aleats){
      ### randomize by sites and species!
      resu[[i]] <- terra::app(x, fun = .lyr.sample, fr = fr)
    }

    resu <- terra::rast(resu) # to transform a list in raster

  } else {

    stop("Choose a valid randomization method! The methods currently available are: 'site', 'specie', 'fullspat'.")

  }

  return(resu)

}
