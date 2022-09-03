#' Calculate Weighted Endemism (WE. Williams et al. 1994, Crisp et al. 2001) standardized for species richness
#'
#' The function calculates the WE corrected for species richness. In each null model run, richness is kept constant and range size of each species are randomized. The function provides the mean, standard deviation of all null models and also calculates the standardized effect size (SES).
#'
#' @param x SpatRaster. A SpatRaster of presence-absence.
#' @param aleats numeric. A number indicating how many times the calculation should be repeated.
#' @param filename character. Output filename.
#' @param random character. A character indicating what type of randomization. Could be by tip, site, specie or full spat(site and specie).
#' @return SpatRaster
#' @export
#' @references Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J., VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and endemism in conservation. In: Systematics and Conservation Evaluation (eds Forey PL, Humphries CJ, Vane-Wright RI), p. 438. Oxford University Press, Oxford.
#' @references Crisp, M., Laffan, S., Linder, H., Monro, A. (2001). Endemism in theAustralian flora. Journal of Biogeography, 28, 183–198.
#' @examples
#' \dontrun{
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' t <- rast.we.ses(ras, aleats = 10)
#' }

rast.we.ses <- function(x, aleats,
                        random = c("area.size", "site", "specie", "fullspat"),
                        filename = NULL){

    aleats <- aleats # number of randomizations
    temp <- vector("list", length = aleats) # "temp" is to store each round of the null model
    rs <- range.size(x) # Calculating area


  ## Null model (bootstrap structure)
    if(random == "area.size"){

      we.rand <- list() # to store the rasters in the loop

      for(i in 1:aleats){
        rs.random <- sample(rs, replace = TRUE) # randomize the range size
        ## check if the values are differents
        # rs == rs.random

        temp[[i]] <- paste0(tempfile(), i, ".tif") # directory to store the rasters

        we.rand[[i]] <- rast.we(pres.rast = x, range.size = rs.random,
                                filename = temp[[i]])
      }

      we.rand <- terra::rast(we.rand) # to transform a list in raster

    } else if(random == "site"){

      we.rand <- list() # to store the rasters in the loop

      for(i in 1:aleats){
        # temporary names to rasters
        temp[[i]] <- paste0(tempfile(), i, ".tif")

        ### embaralha por lyr - ordem dos sítios de cada espécie separada
        site.rand <- spat.rand(x, aleats = 1, random = "site")

        we.rand[[i]] <- rast.we(site.rand, range.size = rs,
                                filename = temp[[i]])
      }

      we.rand <- terra::rast(we.rand) # to transform a list in raster

    } else if(random == "specie"){

      ### randomize by cells - species in each site
      we.rand <- list() # to store the rasters in the loop

      for(i in 1:aleats){
        temp[[i]] <- paste0(tempfile(), i, ".tif") # temporary names to rasters

        sp.rand <- spat.rand(x, aleats = 1, random = "specie")

        we.rand[[i]] <- rast.we(sp.rand, range.size = rs,
                                filename = temp[[i]])
      }

      we.rand <- terra::rast(we.rand) # to transform a list in raster

    } else if (random == "fullspat") {

      we.rand <- list() # to store the rasters in the loop
      fr <- terra::freq(x)

      for(i in 1:aleats){

        temp[[i]] <- paste0(tempfile(), i, ".tif") # temporary names to rasters

        ### randomize sites and species
        full.rand <- spat.rand(x, aleats = 1, random = "fullspat")

        we.rand[[i]] <- rast.we(full.rand, range.size = rs, filename = temp[[i]])

      }

      we.rand <- terra::rast(we.rand) # to transform a list in raster

    } else {
      stop("Choose a valid randomization method! The methods currently available are: 'site', 'specie', 'fullspat'.")
    }

  ## WE observed
  we.obs <- rast.we(pres.rast = x, range.size = rs, filename = filename)

  ## WE rand mean and WE rand SD
  we.rand.mean <- terra::mean(we.rand, na.rm = TRUE, filename = filename) # mean pd
  we.rand.sd <- terra::stdev(we.rand, na.rm = TRUE, filename = filename) # sd pd

  unlink(temp) # delete the archive that will not be used

  ## Calculating the standard effect size (SES)
  {
    ses <- function(x){
      (x[1] - x[2])/sqrt(x[3])
    }
    we.ses <- terra::app(c(we.obs, we.rand.mean, we.rand.sd), fun = ses)
    names(we.ses) <- "SES"
  }
  out <- c(we.obs, we.rand.mean, we.rand.sd, we.ses)
  names(out) <- c("WE Observed", "Mean", "SD", "SES" )
  return(out)
}
