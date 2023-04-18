#' Calculate weighted endemism for each raster cell
#'
#' @description Calculate the sum of the inverse of the range size for species present in each raster cell.
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for a set of species.
#' @param range_size numeric. A numerical vector containing the range size for each species. See the function range_size.
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster with that many cores is created and used.
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @author Neander Marcel Heming and Gabriela Alves-Ferreira
#' @return SpatRaster
# #' @export
#' @references Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J., VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and endemism in conservation. In: Systematics and Conservation Evaluation (eds Forey PL, Humphries CJ, Vane-Wright RI), p. 438. Oxford University Press, Oxford.
#' @references Crisp, M., Laffan, S., Linder, H., Monro, A. (2001). Endemism in theAustralian flora. Journal of Biogeography, 28, 183–198.
#' @examples
#' \dontrun{
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' rs <- range_size(ras)
#' .rast.we.B(ras, rs)
#' }
#'
.rast.we.B <- function(x, range_size, filename = NULL, cores = 1, ...){

  temp <- vector("list", length = 2) # to create a temporary vector with the raster number
  temp[[1]] <- paste0(tempfile(), ".tif")  # to store the first raster

  # inverse of range size
  inv.R <- terra::ifel(x == 0, 0, 1/(x * range_size),
                       filename = temp[[1]],
                       overwrite = TRUE) # calculating the inverse of range size
  # weighted endemism
  { # calculating we
    rend <- terra::app(inv.R,
                       function(x){
                         if(all(is.na(x))){
                           return(NA)}
                         sum(x, na.rm = T)
                       }, cores = cores)
  }
  names(rend) <- "WE" # layer name

  if (!is.null(filename)){ # to save the rasters when the path is provide
    rend <- terra::writeRaster(rend, filename, ...)
  }
  unlink(temp[[1]]) # delete the archive

  return(rend)
}


#' Calculate weighted endemism standardized for species richness
#'
#' @description Calculates the standardized effect size for weighted endemism. The function has three different methods for spatial randomization. See Details for more information.
#' @param x  A SpatRaster containing presence-absence data (0 or 1) for a set of species.
#' @param aleats positive integer. A positive integer indicating how many times the calculation should be repeated.
#' @param random character. A character indicating what type of randomization. Could be by "site", "species" or "both" (site and species).
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster with that many cores is created and used.
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @return SpatRaster
#' @details The "site" method keeps species richness constant in each pixel, but randomizes the position of species in the stack. In the "species" method, the order of species in the stack is kept constant, but the pixels where each species is present are randomized in space. The third method, "full.spat", combines site and species randomization at the same time.
#' @export
#' @references Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J., VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and endemism in conservation. In: Systematics and Conservation Evaluation (eds Forey PL, Humphries CJ, Vane-Wright RI), p. 438. Oxford University Press, Oxford.
#' @references Crisp, M., Laffan, S., Linder, H., Monro, A. (2001). Endemism in theAustralian flora. Journal of Biogeography, 28, 183–198.
#' @examples
#' \dontrun{
# ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
# t <- rast.we.ses(ras, aleats = 10, random = "site")
#' }
rast.we.ses <- function(x, aleats,
                        random = c("site", "species", "both"),
                        cores = 1, filename = NULL, ...){

  aleats <- aleats # number of randomizations
  temp <- vector("list", length = aleats) # "temp" is to store each round of the null model
  rs <- range(x) # Calculating area

  # x rasters will be generated in this function, let's see if there is enough memory in the user's pc
  mi <- .fit.memory(x)

  temp.raster <- paste0(tempfile(), ".tif") # temporary names to rasters

  ## Null model (bootstrap structure)

  # if(random == "area.size"){
  #
  #   we.rand <- list() # to store the rasters in the loop
  #
  #   for(i in 1:aleats){
  #     rs.random <- sample(rs, replace = TRUE) # randomize the range size
  #     ## check if the values are differents
  #     # rs == rs.random
  #     temp[[i]] <- paste0(tempfile(), i, ".tif") # directory to store the rasters
  #     we.rand[[i]] <- .rast.we.B(x, range_size = rs.random,
  #                             filename = temp[[i]], cores = cores)
  #   }
  #
  #     we.rand <- terra::rast(we.rand) # to transform a list in raster
  #
  #   } else if(random == "site"){

  if(random == "site"){

    we.rand <- list() # to store the rasters in the loop

    for(i in 1:aleats){
      # temporary names to rasters
      temp[[i]] <- paste0(tempfile(), i, ".tif")
      ### embaralha por lyr - ordem dos sítios de cada espécie separada
      site.rand <- spat.rand(x, random = "site", filename = temp.raster, memory = mi)
      we.rand[[i]] <- .rast.we.B(site.rand, range_size = rs,
                                 filename = temp[[i]], cores = cores)
    }

    we.rand <- terra::rast(we.rand) # to transform a list in raster

  } else if(random == "species"){

    ### randomize by cells - species in each site
    we.rand <- list() # to store the rasters in the loop

    for(i in 1:aleats){
      temp[[i]] <- paste0(tempfile(), i, ".tif") # temporary names to rasters
      sp.rand <- spat.rand(x, random = "species", filename = temp.raster, memory = mi)
      we.rand[[i]] <- .rast.we.B(sp.rand, range_size = rs,
                                 filename = temp[[i]], cores = cores)
    }

    we.rand <- terra::rast(we.rand) # to transform a list in raster

  } else if (random == "both") {

    we.rand <- list() # to store the rasters in the loop
    fr <- terra::freq(x)

    for(i in 1:aleats){
      temp[[i]] <- paste0(tempfile(), i, ".tif") # temporary names to rasters
      ### randomize sites and species
      full.rand <- spat.rand(x, random = "both", filename = temp.raster, memory = mi)
      we.rand[[i]] <- .rast.we.B(full.rand, range_size = rs,
                                 filename = temp[[i]], cores = cores)
    }

    we.rand <- terra::rast(we.rand) # to transform a list in raster

  } else {
    stop("Choose a valid randomization method! The methods currently available are: 'site', 'species', 'both'.")
  }

  ## WE observed
  we.obs <- rast.we(x = x, filename = filename)

  ## WE rand mean and WE rand SD
  we.rand.mean <- terra::mean(we.rand, na.rm = TRUE, filename = filename) # mean pd
  we.rand.sd <- terra::stdev(we.rand, na.rm = TRUE, filename = filename) # sd pd

  unlink(temp) # delete the file that will not be used

  ## Calculating the standard effect size (SES)
  {
    ses <- function(x){
      (x[1] - x[2])/x[3]
    }
    we.ses <- terra::app(c(we.obs, we.rand.mean, we.rand.sd), fun = ses,
                         cores = cores)
  }

  names(we.ses) <- "SES"
  out <- c(we.rand.mean, we.rand.sd, we.obs, we.ses)
  names(out) <- c("Mean", "SD", "WE Observed", "SES")

  if(!is.null(filename)){ # to save the rasters when the path is provide
    out <- terra::writeRaster(out, filename)
  }

  return(out)
}
