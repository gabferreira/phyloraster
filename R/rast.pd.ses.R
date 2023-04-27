#' Calculate phylogenetic diversity for each raster cell
#'
#' @description Calculate the sum of the branch length for species present in each cell of the raster.
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for a set of species. The layers (species) must be sorted according to the tree order. See the phylo.pres function.
#' @param branch.length numeric. A Named numerical vector containing the branch length for a set of species.
#' @param filename character. Output filename.
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster with that many cores is created and used.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @author Neander Marcel Heming and Gabriela Alves-Ferreira
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological conservation, 61(1), 1-10.
#' @return SpatRaster
# #' @export
#' @examples
#' \dontrun{
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
#' data <- phylo.pres(x, tree)
#' rast.pd(data$x, data$branch.length)
#' }
.rast.pd.B <- function(x, branch.length, filename = NULL, cores = 1, ...){

  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  # if(!all.equal(names(x), names(branch.length))){
  #
  #   stop("Species names are not in the same order on 'x' and 'branch.length' arguments! See 'phyloraster::phylo.pres' function.")
  #
  # } else {

    # 1 rasters will be generated in this function, let's see if there is enough memory in the user's pc
    sink(nullfile())    # suppress output
    mi <- terra::mem_info(x, 1)[5] != 0 # proc in memory = T TRUE means that it fits in the pc's memory, so you wouldn't have to use temporary files
    sink()

    temp <- vector("list", length = 1) # to create a temporary vector with the raster number
    temp[[1]] <- paste0(tempfile(), ".tif")  # to store the first raster

    # phylogenetic diversity
    rpd <- terra::app(x, fun = .vec.pd,
                      branch.length = branch.length, cores = cores,
                      filename = ifelse(mi, "", temp[[1]]))
    rpd <- rpd[[1]] # select only the first raster
    names(rpd) <- c("PD")
  # }

  if(!is.null(filename)){ # to save the rasters when the output filename is provide
    rpd <- terra::writeRaster(rpd, filename)
  }

  return(rpd)
}

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
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
#' data <- phyloraster::phylo.pres(x, tree)
#' t <- phyloraster::rast.pd.ses(data$x, data$branch.length, aleats = 10, random = "species")
#' plot(t)
#' }
rast.pd.ses <- function(x, branch.length, aleats,
                        random = c("tip", "spat"),
                        cores = 1, filename = NULL, ...){

  aleats <- aleats # number of null models
  temp <- vector("list", length = aleats) # to create a temporary vector with the raster number

  # x rasters will be generated in this function, let's see if there is enough memory in the user's pc
  mi <- .fit.memory(x)

  temp.raster <- paste0(tempfile(), ".tif") # temporary names to rasters

  ## Null model (bootstrap structure)
  if(random == "tip"){

    pd.rand <- list() # to store the rasters in the loop
    bl.random <- branch.length # to store the branch length in the loop

    for(i in 1:aleats){

      bl.random[] <- sample(branch.length) # randomize branch lengths
      temp[[i]] <- paste0(tempfile(), i, ".tif") # temporary names to rasters

      pd.rand[[i]] <- .rast.pd.B(x, branch.length = bl.random,
                                 filename = temp[[i]],
                                 cores = cores, ..., overwrite = T)
    }

    pd.rand <- terra::rast(pd.rand) # to transform a list in raster

  } else if(random == "spat"){

    pd.rand <- list() # to store the rasters in the loop
    rich <- rast.se(x)
    prob <- terra::app(x,
                       function(x){
                         ifelse(is.na(x), 0, 1)
                       })

    for(i in 1:aleats){

      # temporary names to rasters
      temp[[i]] <- paste0(tempfile(), i, ".tif")

      ### shuffle
      pres.site.null <- SESraster::bootspat_str(x = x, rich = rich,
                                                prob = prob)

      # calculate pd
      pd.rand[[i]] <- .rast.pd.B(pres.site.null, branch.length = branch.length,
                                 filename = temp[[i]],
                                 cores = cores, ..., overwrite = T)
    }

    pd.rand <- terra::rast(pd.rand) # to transform a list in raster

  } else {
    stop("Choose a valid randomization method! The methods currently available are: 'tip', 'site', 'species', 'both'.")
  }

  ### PD rand mean
  pd.rand.mean <- terra::mean(pd.rand, na.rm = TRUE) # mean pd
  ### PD rand SD
  pd.rand.sd <- terra::stdev(pd.rand, na.rm = TRUE) # sd pd

  ### PD observed
  {
    x.reord <- x[[names(branch.length)]] # to reorder the stack according to the tree

    pd.obs <- phyloraster::rast.pd(x.reord, branch.length, cores = cores)
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
