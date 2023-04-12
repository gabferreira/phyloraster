#' Calculate phylogenetic endemism for a raster
#'
#' Calculate phylogenetic endemism using rasters as input and output.
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for a set of species. The layers (species) must be sorted according to the tree order. See the phylo.pres function.
#' @param branch.length numeric. A named numerical vector containing the branch length of each specie.
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster with that many cores is created and used.
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @return SpatRaster
# #' @export
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C., & Cook, L. G. (2009). Phylogenetic endemism: a new approach for identifying geographical concentrations of evolutionary history. Molecular ecology, 18(19), 4061-4072.
#' @examples
#' \dontrun{
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
#' data <- phylo.pres(ras, tree)
#' .rast.pe.B(data$x, data$branch.length, cores = 1)
#' }
.rast.pe.B <- function(x, branch.length, cores = 1, filename = NULL, ...){

  # if(!all.equal(names(x), names(branch.length))){
  #
  #   stop("Species names are not in the same order on 'x' and 'branch.length' arguments! See 'phylogrid::phylo.pres' function.")
  #
  # } else {

    area.branch <- phylogrid::inv.range(x, branch.length)

    rpe <- terra::app(area.branch$LR,
                      function(x){
                        if(all(is.na(x))){
                          return(NA)
                        }
                        sum(x, na.rm = T)
                      }, cores = cores)
    rpe <- terra::app(rpe, function(x, m){ # to reescale values from 0 to 1
      (x/m)
    }, m = terra::minmax(rpe)[2,], cores = cores)
  # }

  names(rpe) <- c("PE")

  if (!is.null(filename)){ # to save the rasters when the path is provide
    rpe <- terra::writeRaster(rpe, filename, ...)
  }
  return(rpe)
}

#' Phylogenetic endemism standardized for specie richness
#'
#' @description Calculates the standardized effect size for phylogenetic endemism. The function has four different methods for spatial and phylogenetic randomization. See Details for more information.
#' @inheritParams rast.pd.ses
#' @return SpatRaster
#' @author Gabriela Alves-Ferreira and Neander Heming
#' @export
#' @details The "tip" method shuffles the taxon names among all those in the phylogeny. The "site" method keeps species richness constant in each pixel, but randomizes the position of species in the stack. In the "species" method, the order of species in the stack is kept constant, but the pixels where each species is present are randomized in space. The third method, "full.spat", combines site and species randomization at the same time.
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C., & Cook, L. G. (2009). Phylogenetic endemism: a new approach for identifying geographical concentrations of evolutionary history. Molecular ecology, 18(19), 4061-4072.
#' @examples
#' \dontrun{
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
#' data <- phylo.pres(ras, tree)
#' t <- rast.pe.ses(data$x, data$branch.length, aleats = 10, random = "both")
#' plot(t)
#' }
rast.pe.ses <- function(x, branch.length, aleats,
                        random = c("area.size", "site", "species", "both"),
                        cores = 1, filename = "", ...){

  aleats <- aleats # number of null models
  temp <- vector("list", length = aleats) # to create a temporary vector with the raster number


  # x rasters will be generated in this function, let's see if there is enough
  # memory in the user's pc
  mi <- .fit.memory(x)
  temp.raster <- paste0(tempfile(), ".tif") # temporary names to rasters

  ## Null model (bootstrap structure)
  if(random == "tip"){

    pe.rand <- list() # to store the rasters in the loop
    bl.random <- branch.length

    for(i in 1:aleats){
      bl.random[] <- sample(branch.length, replace = T) # aleatorize the branch lenght
      ## check if the values are differents
      # branch.length == bl.random
      temp[[i]] <- paste0(tempfile(), i, ".tif") # directory to store the rasters
      pe.rand[[i]] <- .rast.pe.B(x, branch.length = bl.random,
                                 filename = temp[[i]], cores = cores)
    }

    pe.rand <- terra::rast(pe.rand) # to transform a list in raster

  } else if (random == "site"){

    pe.rand <- list() # to store the rasters in the loop

    for(i in 1:aleats){
      # temporary names to rasters
      temp[[i]] <- paste0(tempfile(), i, ".tif")
      ### shuffle by layer - order of sites for each separate species
      pres.site.null <- spat.rand(x, random = "site", cores = cores,
                                  filename = temp.raster, memory = mi, overwrite = T)
      # calculate pe
      pe.rand[[i]] <- .rast.pe.B(pres.site.null, branch.length = branch.length,
                                 filename = temp[[i]], cores = cores)
    }

    pe.rand <- terra::rast(pe.rand) # to transform a list in raster

  } else if (random == "species") {

    ### randomize by cells - species in each site
    pe.rand <- list() # to store the rasters in the loop

    for(i in 1:aleats){
      temp[[i]] <- paste0(tempfile(), i, ".tif") # temporary names to rasters
      sp.rand <- spat.rand(x, random = "species", cores = cores,
                           filename = temp.raster, memory = mi)
      pe.rand[[i]] <- .rast.pe.B(sp.rand, branch.length = branch.length,
                                 filename = temp[[i]], cores = cores)
    }

    pe.rand <- terra::rast(pe.rand) # to transform a list in raster

  } else if (random == "both") {

    pe.rand <- list() # to store the rasters in the loop
    fr <- terra::freq(x)

    for(i in 1:aleats){
      temp[[i]] <- paste0(tempfile(), i, ".tif") # temporary names to rasters
      ### randomize sites and species
      pres.null <- spat.rand(x, random = "both", cores = cores,
                             filename = temp.raster, memory = mi)
      pe.rand[[i]] <- .rast.pe.B(pres.null, branch.length = branch.length,
                                 filename = temp[[i]], cores = cores)
    }

    pe.rand <- terra::rast(pe.rand) # to transform a list in raster

  } else {
    stop("Choose a valid randomization method! The methods currently available are: 'tip', 'site', 'species', 'both'.")
  }

  ## PE observed
  x.reord <- x[[names(branch.length)]] # to reorder the stack according to the tree

  pe.obs <- phylogrid::rast.pe(x.reord, branch.length = branch.length,
                               filename = filename, cores = cores)

  ## PD rand mean and PD rand SD
  pe.rand.mean <- terra::mean(pe.rand, na.rm = TRUE, filename = filename) # mean pd
  pe.rand.sd <- terra::stdev(pe.rand, na.rm = TRUE, filename = filename) # sd pd

  unlink(temp) # delete the file that will not be used
  unlink(temp.raster) # delete the file that will not be used

  ## Calculating the standard effect size (SES)
  {
    ses <- function(x){
      (x[1] - x[2])/x[3]
    }
    pe.ses <- terra::app(c(pe.obs, pe.rand.mean, pe.rand.sd),
                         fun = ses, cores = cores)
  }

  names(pe.ses) <- "SES"
  out <- c(pe.rand.mean, pe.rand.sd, pe.obs, pe.ses)
  names(out) <- c("Mean", "SD", "PE Observed", "SES")

  if (!is.null(filename)){ # to save the rasters when the path is provide
    out <- terra::writeRaster(out, filename)
  }

  return(out)
}
