#' Intern function to resampling a vector according to the observed frequency
#'
#' @param x numeric. A vector containing values to resampling.
#' @param fr data.frame A data.frame with 3 columns (layer, value, count).
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

#' Intern function to sample vectors with non-NA values
#'
#' @param x  numeric. A vector containing values to resampling.
#'
#' @return vector
# #' @export
#'
#' @examples
.sample.not.NA <- function(x){
  s <- !is.na(x)
  x[s] <- sample(x[s])
  x
}

#' Randomize a set of rasters according to the observed frequency.
#'
#' Randomize a set of rasters according to the observed frequency using the methods: sites (by cells), species (by layer) or both (layers and cells).
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster with that many cores is created and used.
#' @param x SpatRaster. A presence-absence SpatRaster.
#' @param random character. Character indicating the type of randomization to be used. The available types are by "site", "specie" or "fullspat".
#' @return SpatRaster
#' @author Neander Marcel Heming and Gabriela Alves-Ferreira
#' @export
#'
#' @examples
#' \dontrun{
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' sr <- spat.rand(ras, random = "site")
#' plot(sr)
#' }
spat.rand <- function(x, random = c("site", "species", "fullspat"),
                      filename = "", memory = FALSE, cores = 1){

  fr <- terra::freq(x)

  if (random == "species"){
    if(memory){ # quando cabe na memoria

      # randomize by layers- sites
      resu <- terra::rast(lapply(1:terra::nlyr(x),
                                 function(i, r){
                                   terra::app(r[[i]],
                                              fun = .sample.not.NA)
                                 }, r = x))
    } else { # nao cabe na memoria

      # randomize by layers- sites
      resu <- terra::writeRaster(terra::rast(lapply(1:terra::nlyr(x),
                                                    function(i, r, fr){
                                                      terra::app(r[[i]], fun = .lyr.sample,
                                                                 fr = fr[fr$layer==i,])
                                                    }, r = x, fr = fr)),
                                 filename = filename)


    }

  } else if (random == "site") {

    if(memory){
      ### randomize by cells- species in each site
      resu <- terra::app(x, sample, cores = cores, filename = filename)

    } else {
      ### randomize by cells- species in each site
      resu <- terra::app(x, sample, cores = cores, filename = filename)

    }

  } else if (random == "fullspat") {

    if(memory){
      ### randomize by sites and species!
      resu <- x
      resu[] <- .sample.not.NA(x[])

    } else {
      ### randomize by sites and species!
      resu <- terra::app(x, fun = .lyr.sample, fr = fr, cores = cores,
                         filename = filename)

    }

  } else {
    stop("Choose a valid randomization method! The methods currently available are: 'site', 'species', 'fullspat'.")
  }

  return(resu)

}

