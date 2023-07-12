#' Calculate phylogenetic diversity (Faith 1992) for a vector
#'
#' @description This function calculates the sum of the branch length for a set of species for one sample.
#' @param x numeric. A named numerical vector of presence-absence for one sample.
#' @param branch.length numeric. A named numerical vector containing the branch length for each species.
#' @author Neander Marcel Heming and Gabriela Alves-Ferreira
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological conservation, 61(1), 1-10.
#' @return numeric
# #' @export
.vec.pd <- function(x, branch.length, pd = c(PD=NA)){

  x[is.na(x)] <- 0 # 0 for all value = NA

  if(sum(x)== 0) { # return NA if x = 0
    return(pd)
  }

  if(sum(x) != 0) { # if the sum of x is non-zero then do this:
    # x == 1, only species present in the vector
    pd[] <- sum(branch.length[x == 1]) # pd Faith 1992
  }

  return(pd)
}


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
.rast.pd.B <- function(x, branch.length, cores = 1, filename = "", ...){

  # phylogenetic diversity
  rpd <- terra::app(x,
                    .vec.pd,
                    branch.length = branch.length,
                    cores = cores, filename = filename, ...)

  terra::set.names(rpd, "PD")

  return(rpd)
}


#' Calculate phylogenetic diversity for raster data
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
#' @examples
#' library(terra)
#' library(phyloraster)
#' x <- rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
#' data <- phylo.pres(x, tree)
#' pd <- rast.pd(data$x, data$branch.length)
#' plot(pd)
#' @export
rast.pd <- function(x, branch.length, cores = 1, filename = "", ...){

  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  if(!all.equal(names(x), names(branch.length))){
    stop("Species names are not in the same order on 'x' and 'branch.length' arguments! See 'phyloraster::phylo.pres' function.")
  }

  rpd <- .rast.pd.B(x, branch.length, cores, filename, ...)

  return(rpd)

}


#' Phylogenetic diversity standardized for species richness
#'
#' @description Calculates the standardized effect size for phylogenetic diversity. See Details for more information.
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for a set of species.
#' @param branch.length numeric. A named numerical vector containing the branch length of each species. See phylo.pres function.
#' @param aleats positive integer. A positive integer indicating how many times the calculation should be repeated.
#' @param random character. A character indicating the type of randomization. The currently available randomization methods are "tip", "site", "species" or "both" (site and species).
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster with that many cores is created and used.
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @return SpatRaster. The function returns the observed phylogenetic diversity, the mean of the simulations calculated over n times, the standard deviation of the simulations, and the standardized effect size (SES) for the phylogenetic diversity.
#' @details The spatial randomization (spat) keeps the richness exact and samples species presences proportionally to their observed frequency (i.e. number of occupied pixels). The randomization will not assign values to cells with nodata. The phylogenetic randomization shuffles taxa names across all taxa included in phylogeny.
#' @export
#' @author Gabriela Alves-Ferreira and Neander Heming
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological conservation, 61(1), 1-10.
#' @examples
#' \dontrun{
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
#' data <- phyloraster::phylo.pres(x, tree)
#' t <- phyloraster::rast.pd.ses(data$x, data$branch.length, aleats = 10, random = "spat")
#' plot(t)
#' }
rast.pd.ses <- function(x, branch.length, aleats,
                        random = c("tip", "spat"),
                        cores = 1, filename = "", ...){

  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  if(!all.equal(names(x), names(branch.length))){
    stop("Species names are not in the same order on 'x' and 'branch.length' arguments! See 'phyloraster::phylo.pres' function.")
  }

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

    fr_prob <- SESraster::fr2prob(x)

    for(i in 1:aleats){
      # temporary names to rasters
      temp[[i]] <- paste0(tempfile(), i, ".tif")

      ### shuffle
      pres.site.null <- SESraster::bootspat_str(x = x, rich = rich, prob = prob,
                                                fr_prob = fr_prob)
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

  unlink(temp) # delete the archive that will not be used
  unlink(temp.raster) # delete the archive that will not be used

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

  return(out)
}
