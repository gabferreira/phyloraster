#' Calculate phylogenetic endemism for a vector
#'
#' @description This function calculates phylogenetic endemism for a vector
#' @param x numeric. A named numerical vector of presence-absence for one sample.
#' @param branch.length numeric. A named numerical vector containing the branch length for each species.
#' @author Neander Marcel Heming and Gabriela Alves-Ferreira
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological conservation, 61(1), 1-10.
#' @return numeric
# #' @export
.vec.wpe <- function(x, spp_seq, spp_seq2, wpe = c(PE=NA)){

  if(all(is.na(x[spp_seq]))){
    return(wpe)
  }

  wpe[] <- sum(x[spp_seq]*x[spp_seq2], na.rm = T)

  return(wpe)
}

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
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
#' data <- phylo.pres(x, tree)
#' ### TODO correct example
#' .rast.pe.B(data$x, data$branch.length, cores = 1)
#' }
.rast.pe.B <- function(x, spp_seq, spp_seqLR, cores = 1, filename = "", ...){

  # phylogenetic endemism
  rpe <-terra::app(x,
                   .vec.wpe,
                   spp_seq, spp_seqLR,
                   cores = cores, filename = filename, ...)

  terra::set.names(rpe, "PE") # layer name

  # if(rescale == TRUE){
  #   rpe <- terra::app(rpe, function(x, m){ # rescale the values from 0 to 1
  #     (x/m)
  #   }, m = terra::minmax(rpe)[2,], filename = ifelse(mi, "", temp[[3]]))
  #
  # }

  return(rpe)

}

#' Calculate phylogenetic endemism for raster data
#'
#' @description Calculate the sum of the inverse of the range size multiplied by the branch length for the species present in raster data.
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for a set of species. The layers (species) must be sorted according to the tree order. See the phylo.pres function.
#' @param branch.length numeric. A Named numeric vector containing the branch length of each specie
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster with that many cores is created and used.
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @author Gabriela Alves-Ferreira and Neander Marcel Heming
#' @references Laffan, S. W., Rosauer, D. F., Di Virgilio, G., Miller, J. T., González‐Orozco, C. E., Knerr, N., ... & Mishler, B. D. (2016). Range‐weighted metrics of species and phylogenetic turnover can better resolve biogeographic transition zones. Methods in Ecology and Evolution, 7(5), 580-588.
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C. and Cook, L. G. (2009). Phylogenetic endemism: a new approach for identifying geographical concentrations of evolutionary history. Molecular ecology, 18(19), 4061-4072.
#' @return SpatRaster
#' @examples
#' library(terra)
#' library(phyloraster)
#' x <- rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
#' data <- phylo.pres(x, tree)
#' pe <- rast.pe(data$x, data$branch.length)
#' plot(pe)
#'
#' @export
rast.pe <- function(x, branch.length, cores = 1, filename = "", ...){

  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  # 3 rasters will be generated in this function, let's see if there is enough memory in the user's pc
  if(!all.equal(names(x), names(branch.length))){
    stop("Species names are not in the same order on 'x' and 'branch.length' arguments! See 'phyloraster::phylo.pres' function.")
  }

  mi <- .fit.memory(x)

  temp <- vector("list", length = 3) # to create a temporary vector with the raster number
  temp[[1]] <- paste0(tempfile(), ".tif")  # to store the first raster

  LR <- inv.range(x, branch.length, # LR = T,
                  filename = ifelse(mi, "", temp[[1]]))$LR # calculate the inverse of range size multiplied by branch length of each species


  nspp <- terra::nlyr(x)
  spp_seq <- seq_len(nspp)
  spp_seqLR <- spp_seq + nspp

  # phylogenetic endemism
  rpe <- .rast.pe.B(c(x, LR), spp_seq, spp_seqLR, cores, filename, ...)

  unlink(temp[[1]])

  return(rpe)

}


#' Phylogenetic endemism standardized for specie richness
#'
#' @description Calculates the standardized effect size for phylogenetic endemism. See Details for more information.
#' @inheritParams rast.pd.ses
#' @return SpatRaster
#' @author Gabriela Alves-Ferreira and Neander Heming
#' @export
#' @details The spatial randomization (spat) keeps the richness exact and samples species presences proportionally to their observed frequency (i.e. number of occupied pixels). The randomization will not assign values to cells with nodata. The phylogenetic randomization shuffles taxa names across all taxa included in phylogeny.
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C., & Cook, L. G. (2009). Phylogenetic endemism: a new approach for identifying geographical concentrations of evolutionary history. Molecular ecology, 18(19), 4061-4072.
#' @examples
#' \dontrun{
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
#' data <- phylo.pres(x, tree)
#' t <- rast.pe.ses(data$x, data$branch.length, aleats = 10, random = "spat")
#' plot(t)
#' }
rast.pe.ses <- function(x, branch.length, aleats,
                        random = c("tip", "spat"),
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

  } else if (random == "spat"){

    pe.rand <- list() # to store the rasters in the loop
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

      # calculate pe
      pe.rand[[i]] <- .rast.pe.B(pres.site.null, branch.length = branch.length,
                                 filename = temp[[i]], cores = cores)
    }

    pe.rand <- terra::rast(pe.rand) # to transform a list in raster
  } else {
    stop("Choose a valid randomization method! The methods currently available are: 'tip', 'spat'.")
  }

  ## PE observed
  x.reord <- x[[names(branch.length)]] # to reorder the stack according to the tree

  pe.obs <- rast.pe(x.reord, branch.length = branch.length,
                    filename = filename, cores = cores)

  ## PE rand mean and PE rand SD
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
