#' Calculate Evolutionary Distinctiveness for a vector
#'
#' @description This function calculates evolutionary distinctiveness for a set of species using the fair-proportion index (Isaac et al., 2007).
#' @usage .evol.distin(x, branch.length, n.descen)
#' @param x numeric. A Named numeric vector of presence-absence
#' @param branch.length numeric. A Named numeric vector of branch length for each specie
#' @param n.descen numeric. A Named numeric vector of number of descendants for each branch
#' @author Gabriela Alves-Ferreira and Neander Marcel Heming
#' @references Isaac, N. J., Turvey, S. T., Collen, B., Waterman, C. and Baillie, J. E. (2007). Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE 2, e296.
#' @return numeric
# #' @export
.evol.distin <- function(x, branch.length, n.descen){

  x[is.na(x)] <- 0 # 0 for all value = NA

  if(sum(x) == 0) { # return NA if x = 0
    return(c(NA,NA))
  }

  if(sum(x) != 0){ # if the sum of x is non-zero then do this:
    pres <- x == 1 # only species present in the vector
    # species <- names(x)
    ed <- sum(branch.length[pres]/n.descen[pres]) # evolutionary distinctiveness
    ed1 <- ed # terra::app function does not work when this intern function returns only one raster
  }

  return(c(ed = ed, ed1 = ed1))

}

#' Calculate Evolutionary distinctiveness for each raster cell
#'
#' @description This function calculates evolutionary distinctiveness according to the fair-proportion index.
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for a set of species. The layers (species) must be sorted according to the tree order. See the phylo.pres function.
#' @param branch.length numeric. A Named numeric vector containing the branch length of each specie.
#' @param n.descen numeric. A Named numeric vector of number of descendants for each branch
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster with that many cores is created and used.
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed for fun.
#' @author Gabriela Alves-Ferreira and Neander Marcel Heming
#' @references Isaac, N. J., Turvey, S. T., Collen, B., Waterman, C. and Baillie, J. E. (2007). Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE 2, e296.
#' @return SpatRaster
# #' @export
#' @examples
.rast.ed.B <- function(x, branch.length, n.descen, cores = 1, filename = NULL, ...){

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

    # evolutionary distinctiveness
    red <- terra::app(x, fun = .evol.distin,
                      branch.length, n.descen, cores = cores,
                      filename = ifelse(mi, "", temp[[1]]))
    red <- red[[1]] # only the first raster
    names(red) <- "ED" # layer name
  # }

  if(!is.null(filename)){ # to save the rasters when the output filename is provide
    red <- terra::writeRaster(red, filename)
  }

  return(red)

}

#' Evolutionary distinctiveness standardized for specie richness
#'
#' @description Calculates the standardized effect size for evolutionary distinctiveness. See Details for more information.
#' @inheritParams rast.pd.ses
#' @inheritParams rast.ed
#' @return SpatRaster
#' @author Gabriela Alves-Ferreira and Neander Heming
#' @export
#' @details The spatial randomization (spat) keeps the richness exact and samples species presences proportionally to their observed frequency (i.e. number of occupied pixels). The randomization will not assign values to cells with nodata. The phylogenetic randomization shuffles taxa names across all taxa included in phylogeny.
#' @references Isaac, N. J., Turvey, S. T., Collen, B., Waterman, C. and Baillie, J. E. (2007). Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE 2, e296.
#' @references Laffan, S. W., Rosauer, D. F., Di Virgilio, G., Miller, J. T., González‐Orozco, C. E., Knerr, N., ... & Mishler, B. D. (2016). Range‐weighted metrics of species and phylogenetic turnover can better resolve biogeographic transition zones. Methods in Ecology and Evolution, 7(5), 580-588.
#' @examples
#' \dontrun{
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
#' data <- phyloraster::phylo.pres(x, tree)
#' t <- rast.ed.ses(data$x, data$branch.length, data$n.descendants, aleats = 5, random = "spat")
#' terra::plot(t)
#' }
rast.ed.ses <- function(x, branch.length, n.descen, aleats, random =
                          c("tip", "spat"),
                        cores = 1, filename = NULL, ...){
  aleats <- aleats
  temp <- vector("list", length = aleats) # to create a temporary vector with the raster number

  # x rasters will be generated in this function, let's see if there is enough memory in the user's pc
  mi <- .fit.memory(x)

  temp.raster <- paste0(tempfile(), ".tif") # temporary names to rasters

  ## Null model (bootstrap structure)
  if(random == "tip"){

    ed.rand <- list() # to store the rasters in the loop
    bl.random <- branch.length # to store the branch length in the loop

    for(i in 1:aleats){

      bl.random[] <- sample(branch.length) # randomize branch lengths
      temp[[i]] <- paste0(tempfile(), i, ".tif") # temporary names to rasters

      ed.rand[[i]] <- .rast.ed.B(x, branch.length = bl.random, n.descen,
                                 filename = temp[[i]], cores = cores)
    }

    ed.rand <- terra::rast(ed.rand) # to transform a list in raster

  } else if(random == "spat"){

    ed.rand <- list() # to store the rasters in the loop
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
      pres.site.null <- SESraster::bootspat_str(x = x, rich = rich,
                                                fr_prob = fr_prob, prob = prob)

      # calculate ed
      ed.rand[[i]] <- .rast.ed.B(pres.site.null, branch.length = branch.length,
                                 n.descen = n.descen,
                               filename = temp[[i]], cores = cores)
    }

    ed.rand <- terra::rast(ed.rand) # to transform a list in raster

  }  else {
    stop("Choose a valid randomization method! The methods currently available are: 'tip', 'spat'.")
  }

  ## ed observed
  x.reord <- x[[names(branch.length)]] # to reorder the stack according to the tree

  ed.obs <- phyloraster::rast.ed(x.reord, branch.length = branch.length,
                               n.descen = n.descen,
                               filename = filename, cores = cores)

  ## ED rand mean and ED rand SD
  ed.rand.mean <- terra::mean(ed.rand, na.rm = TRUE, filename = filename) # mean pd
  ed.rand.sd <- terra::stdev(ed.rand, na.rm = TRUE, filename = filename) # sd pd

  unlink(temp) # delete the file that will not be used
  unlink(temp.raster) # delete the file that will not be used

  ## Calculating the standard effect size (SES)
  {
    ses <- function(x){
      (x[1] - x[2])/x[3]
    }
    ed.ses <- terra::app(c(ed.obs, ed.rand.mean, ed.rand.sd),
                         fun = ses, cores = cores)
  }

  names(ed.ses) <- "SES"
  out <- c(ed.rand.mean, ed.rand.sd, ed.obs, ed.ses)
  names(out) <- c("Mean", "SD", "ED Observed", "SES")

  if (!is.null(filename)){ # to save the rasters when the path is provide
    out <- terra::writeRaster(out, filename)
  }

  return(out)
}
