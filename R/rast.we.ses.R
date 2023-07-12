#' Calculate weighted endemism for each raster cell
#'
#' @description Calculate the sum of the inverse of the range size for species present in each raster cell.
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for a set of species.
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
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' rs <- range_size(x)
#' .rast.we.B(x, rs)
#' }
#'
.rast.we.B <- function(x, spp_seq, spp_seqINV, cores = 1, filename = "", ...){

  # weighted endemism
  rend <- terra::app(x,
                     .vec.wpe,
                     spp_seq, spp_seqINV,
                     cores = cores, filename = filename, ...)

  terra::set.names(rend, "WE") # layer name

  return(rend)
}


#' Calculate weighted endemism for raster data
#'
#' @description Calculate the sum of the inverse of the range size for species present in raster data.
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for a set of species.
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster with that many cores is created and used.
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @author Neander Marcel Heming and Gabriela Alves Ferreira
#' @return SpatRaster
#' @references Laffan, S. W., Rosauer, D. F., Di Virgilio, G., Miller, J. T., González‐Orozco, C. E., Knerr, N., ... & Mishler, B. D. (2016). Range‐weighted metrics of species and phylogenetic turnover can better resolve biogeographic transition zones. Methods in Ecology and Evolution, 7(5), 580-588.
#' @references Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J., VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and endemism in conservation. In: Systematics and Conservation Evaluation (eds Forey PL, Humphries CJ, Vane-Wright RI), p. 438. Oxford University Press, Oxford.
#' @references Crisp, M., Laffan, S., Linder, H., Monro, A. (2001). Endemism in theAustralian flora. Journal of Biogeography, 28, 183–198.
#' @examples
#' library(terra)
#' library(phyloraster)
#' x <- rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' we <- rast.we(x)
#' plot(we)
#'
#' @export
rast.we <- function(x, cores = 1, filename = "", ...){

  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  # 2 rasters will be generated in this function, let's see if there is enough memory in the user's pc
  sink(nullfile())    # suppress output
  mi <- terra::mem_info(x, 2)[5] != 0 # proc in memory = T TRUE means that it fits in the pc's memory, so you wouldn't have to use temporary files
  sink()

  # if(!mi){
  temp <- vector("list", length = 2) # to create a temporary vector with the raster number
  temp[[1]] <- paste0(tempfile(), ".tif")  # to store the first raster
  # }

  # inverse of range size
  inv.R <- inv.range(x, filename = ifelse(mi, "", temp[[1]]))$inv.R # calculate the inverse of range size multiplied by branch length of each species

  nspp <- terra::nlyr(x)
  spp_seq <- seq_len(nspp)
  spp_seqINV <- spp_seq + nspp

  # weighted endemism
  rend <- .rast.we.B(c(x, inv.R), spp_seq, spp_seqINV, cores, filename, ...)
  # rend <- .rast.we.B(c(x, inv.R), spp_seq, spp_seqINV, cores)

  # if(rescale == TRUE){
  #   rend <- terra::app(rend, function(x, m){ # rescale the values
  #     (x/m)
  #   }, m = terra::minmax(rend)[2,], cores = cores, filename = ifelse(mi, "", temp[[3]])) # 2 is the max
  # }

  unlink(temp[[1]]) # delete the files

  return(rend)
}


#' Calculate weighted endemism standardized for species richness
#'
#' @description Calculates the standardized effect size for weighted endemism. See Details for more information.
#' @param x  A SpatRaster containing presence-absence data (0 or 1) for a set of species.
#' @param aleats positive integer. A positive integer indicating how many times the calculation should be repeated.
#' @param random character. A character indicating what type of randomization. Could be by "site", "species" or "both" (site and species).
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster with that many cores is created and used.
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @return SpatRaster
#' @details The spatial randomization (spat) keeps the richness exact and samples species presences proportionally to their observed frequency (i.e. number of occupied pixels). The randomization will not assign values to cells with nodata.
#' @export
#' @references Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J., VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and endemism in conservation. In: Systematics and Conservation Evaluation (eds Forey PL, Humphries CJ, Vane-Wright RI), p. 438. Oxford University Press, Oxford.
#' @references Crisp, M., Laffan, S., Linder, H., Monro, A. (2001). Endemism in theAustralian flora. Journal of Biogeography, 28, 183–198.
#' @examples
#' \dontrun{
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' t <- rast.we.ses(data$x, data$branch.length, aleats = 10, random = "spat")
#' plot(t)
#' }
rast.we.ses <- function(x, aleats,
                        random = c("spat"),
                        cores = 1, filename = NULL, ...){

  aleats <- aleats # number of randomizations
  temp <- vector("list", length = aleats) # "temp" is to store each round of the null model

  # x rasters will be generated in this function, let's see if there is enough memory in the user's pc
  mi <- .fit.memory(x)

  temp.raster <- paste0(tempfile(), ".tif") # temporary names to rasters

  ## Null model (bootstrap structure)

  if(random == "spat"){

    we.rand <- list() # to store the rasters in the loop
    rich <- rast.se(x)
    prob <- terra::app(x,
                       function(x){
                         ifelse(is.na(x), 0, 1)
                       })

    fr_prob <- SESraster::fr2prob(x)

    for(i in 1:aleats){
      # temporary names to rasters
      temp[[i]] <- paste0(tempfile(), i, ".tif")

      ### embaralha por lyr - ordem dos sítios de cada espécie separada
      ### shuffle
      site.rand <- SESraster::bootspat_str(x = x, rich = rich, prob = prob,
                                                fr_prob = fr_prob)
      we.rand[[i]] <- .rast.we.B(site.rand,
                                 filename = temp[[i]], cores = cores)
    }

    we.rand <- terra::rast(we.rand) # to transform a list in raster

  } else {
    stop("Choose a valid randomization method! The method currently available is: 'spat'.")
  }

  ## WE observed
  we.obs <- rast.we(x = x, filename = filename)

  ## WE rand mean and WE rand SD
  we.rand.mean <- terra::mean(we.rand, na.rm = TRUE, filename = filename) # mean pd
  we.rand.sd <- terra::stdev(we.rand, na.rm = TRUE, filename = filename) # sd pd

  unlink(temp) # delete the file that will not be used
  unlink(temp.raster) # delete the file that will not be used

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
