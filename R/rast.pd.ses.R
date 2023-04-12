#' Phylogenetic diversity standardized for species richness
#'
#' @description Calculates the standardized effect size for phylogenetic diversity. The function has four different methods for spatial and phylogenetic randomization. See Details for more information.
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for a set of species.
#' @param branch.length numeric. A named numerical vector containing the branch length of each species. See phylo.pres function.
#' @param aleats positive integer. A positive integer indicating how many times the calculation should be repeated.
#' @param random character. A character indicating the type of randomization. The currently available randomization methods are "tip", "site", "species" or "both" (site and species).
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster with that many cores is created and used.
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @return SpatRaster. The function returns the observed phylogenetic diversity, the mean of the simulations calculated over n times, the standard deviation of the simulations, and the standardized effect size (SES) for the phylogenetic diversity.
#' @details The "tip" method shuffles the taxon names among all those in the phylogeny. The "site" method keeps species richness constant in each pixel, but randomizes the position of species in the stack. In the "species" method, the order of species in the stack is kept constant, but the pixels where each species is present are randomized in space. The third method, "full.spat", combines site and species randomization at the same time.
#' @export
#' @author Gabriela Alves-Ferreira and Neander Heming
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological conservation, 61(1), 1-10.
#' @examples
#' \dontrun{
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
#' data <- phylogrid::phylo.pres(ras, tree)
#' t <- phylogrid::rast.pd.ses(data$x, data$branch.length, aleats = 10, random = "species")
#' plot(t)
#' }
rast.pd.ses <- function(x, branch.length, aleats,
                        random = c("tip", "site", "species", "both"),
                        cores = 1, filename = "", ...){

  aleats <- aleats # number of null models
  temp <- vector("list", length = aleats) # to create a temporary vector with the raster number

  # x rasters will be generated in this function, let's see if there is enough memory in the user's pc
  mi <- .fit.memory(x)

  temp.raster <- paste0(tempfile(), ".tif") # temporary names to rasters

  ## Null model (bootstrap structure)
  if(random == "tip"){

    pd.rand <- list() # to store the rasters in the loop
    pd.rand2 <- list() # list to store the rasters
    bl.random <- branch.length # to store the branch length in the loop

    for(i in 1:aleats){

      bl.random[] <- sample(branch.length) # randomize branch lengths
      temp[[i]] <- paste0(tempfile(), i, ".tif") # temporary names to rasters

      pd.rand[[i]] <- terra::app(x, fun = .vec.pd,
                                 branch.length = bl.random,
                                 filename = temp[[i]], cores = cores, ..., overwrite = T)
      pd.rand2[[i]] <- pd.rand[[i]][[1]] # only the first layer for each specie
    }

    pd.rand2 <- terra::rast(pd.rand2) # to transform a list in raster

  } else if(random == "site"){

    pd.rand <- list() # to store the rasters in the loop
    pd.rand2 <- list() # list to store the rasters

    for(i in 1:aleats){

      # temporary names to rasters
      temp[[i]] <- paste0(tempfile(), i, ".tif")

      ### shuffle by layer - order of sites for each separate species
      pres.site.null <- spat.rand(x, random = "site", cores = cores,
                                  filename = temp.raster, memory = mi)

      # calculate pd
      pd.rand[[i]] <- terra::app(pres.site.null, fun = .vec.pd,
                                 branch.length = branch.length,
                                 filename = temp[[i]], cores = cores, ..., overwrite = T)

      pd.rand2[[i]] <- pd.rand[[i]][[1]] # only the first layer for each specie
    }

    pd.rand2 <- terra::rast(pd.rand2) # to transform a list in raster

  } else if(random == "species") {

    ### shuffle by cells - species in each site
    pd.rand <- list() # to store the rasters in the loop
    pd.rand2 <- list() # list to store the rasters

    for(i in 1:aleats){

      temp[[i]] <- paste0(tempfile(), i, ".tif") # temporary names to rasters
      sp.rand <- spat.rand(x, random = "species", cores = cores,
                           filename = temp.raster, memory = mi)
      pd.rand[[i]] <- terra::app(sp.rand,
                                 fun = .vec.pd,
                                 branch.length = branch.length,
                                 filename = temp[[i]], cores = cores, ..., overwrite = T)

      pd.rand2[[i]] <- pd.rand[[i]][[1]] # only the first layer for each specie
    }

    pd.rand2 <- terra::rast(pd.rand2) # to transform a list in raster

  } else if(random == "both") {

    pd.rand <- list() # to store the rasters in the loop
    pd.rand2 <- list() # list to store the rasters
    fr <- terra::freq(x)

    for(i in 1:aleats){

      temp[[i]] <- paste0(tempfile(), i, ".tif") # temporary names to rasters

      ### shuffle sites and species - "full.spat"
      pres.null <- terra::app(x, fun = .lyr.sample, fr = fr, cores = cores,
                              filename = temp.raster, ..., overwrite = T)

      pd.rand[[i]] <- terra::app(pres.null, fun = .vec.pd,
                                 branch.length = branch.length,
                                 filename = temp[[i]], cores = cores, ..., overwrite = T)

      pd.rand2[[i]] <- pd.rand[[i]][[1]] # only the first layer for each specie
    }

    pd.rand2 <- terra::rast(pd.rand2) # to transform a list in raster

  } else {
    stop("Choose a valid randomization method! The methods currently available are: 'tip', 'site', 'species', 'both'.")
  }

  ### PD rand mean
  pd.rand.mean <- terra::mean(pd.rand2, na.rm = TRUE) # mean pd
  ### PD rand SD
  pd.rand.sd <- terra::stdev(pd.rand2, na.rm = TRUE) # sd pd

  ### PD observed
  {
    x.reord <- x[[names(branch.length)]] # to reorder the stack according to the tree

    pd.obs <- phylogrid::rast.pd(x.reord, branch.length, cores = cores)
    pd.obs <- pd.obs$PD # selecting only PD
  }

  ### Concatenate rasters
  pd <- c(pd.rand.mean, pd.rand.sd, pd.obs)

  ## Calculating the standard effect size (SES)
  {
    ses <- function(x){
      (x[1] - x[2])/x[3]
    }
    pd.ses <- terra::app(c(pd.obs, pd.rand.mean, pd.rand.sd),
                         ses, cores = cores, ...)
    names(pd.ses) <- "SES"
  }

  out <- c(pd, pd.ses)
  names(out) <- c("Mean", "SD", "PD Observed", "SES")

  if (!is.null(filename)){ # to save the rasters when the path is provide
    out <- terra::writeRaster(out, filename)
  }

  unlink(temp) # delete the archive that will not be used
  unlink(temp.raster) # delete the archive that will not be used

  return(out)
}

