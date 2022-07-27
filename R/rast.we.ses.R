#' Calculate Weighted Endemism (WE. Williams et al. 1994, Crisp et al. 2001) standardized for species richness
#'
#' The function calculates the WE corrected for species richness. In each null model run, richness is kept constant and range size of each species are randomized. The function provides the mean, standard deviation of all null models and also calculates the standardized effect size (SES).
#'
#' @param pres.rast SpatRaster. A SpatRaster of presence-absence.
#' @param aleats numeric. A number indicating how many times the calculation should be repeated.
#' @return SpatRaster
#' @export
#' @references Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J., VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and endemism in conservation. In: Systematics and Conservation Evaluation (eds Forey PL, Humphries CJ, Vane-Wright RI), p. 438. Oxford University Press, Oxford.
#' @references Crisp, M., Laffan, S., Linder, H., Monro, A. (2001). Endemism in theAustralian flora. Journal of Biogeography, 28, 183–198.
#' @examples
#' \dontrun{
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' t <- rast.we.ses(ras, aleats = 10)
#' }
rast.we.ses <- function(pres.rast, aleats){
  {
    aleats <- aleats # number of null models

    # "temp" is to store each round of the null model
    temp <- vector("list", length = aleats) # to create a temporary vector with the raster number

    # "temp2" is to store the exportable results of the null model
    temp2 <- vector("list", length = 4) # to create a temporary vector
    temp2[[1]] <- paste0(tempfile(), "we.obs.tif") # to save the pd obs
    temp2[[2]] <- paste0(tempfile(), "we.rand.mean.tif") # to save the mean pd
    temp2[[3]] <- paste0(tempfile(), "we.rand.sd.tif") # to save the sd pd
    temp2[[4]] <- paste0(tempfile(), "we.ses.tif") # to save the pd ses
  }

  # Calculating area
  rs <- range.size(pres.rast)

  ## Null model (bootstrap structure)
  {
    we.rand <- list() # to store the rasters in the loop

    for(i in 1:aleats){
      rs.random <- sample(rs, replace = TRUE) # aleatorize the branch lenght
      ## check if the values are differents
      # rs == rs.random

      for(j in 1:aleats){
        temp[[j]] <- paste0(tempfile(), j, ".tif") # directory to store the rasters
      }
      we.rand[[i]] <- rast.we(pres.rast = pres.rast, range.size = rs.random, filename = temp[[i]])
    }

    we.rand <- terra::rast(we.rand) # to transform a list in raster
  }

  ## WE observed
  we.obs <- rast.we(pres.rast = pres.rast, range.size = rs, overwrite=TRUE, filename = temp2[[1]])

  ## WE rand mean and WE rand SD
  we.rand.mean <- terra::mean(we.rand, na.rm = TRUE, overwrite = TRUE, filename = temp2[[2]]) # mean pd
  we.rand.sd <- terra::stdev(we.rand, na.rm = TRUE, overwrite = TRUE, filename = temp2[[3]]) # sd pd

  unlink(temp) # delete the archive that will not be used
  # unlink(temp2) # delete the archive that will not be used

  ## Calculating the standard effect size (SES)
  {
    ses <- function(x, y, z){
      (x - y)/z
    }
    we.ses <- ses(x = we.obs, y = we.rand.mean, z = we.rand.sd)
    names(we.ses) <- "SES"
  }
  out <- c(we.obs, we.rand.mean, we.rand.sd, we.ses)
  names(out) <- c("WE Observed", "Mean", "SD", "SES" )
  return(out)
}

