#' Calculate Evolutionary Distinctiveness for a vector
#'
#' @description This function calculates evolutionary distinctiveness for a set
#' of species using the fair-proportion index (Isaac et al., 2007).
#'
#' @usage .vec.ed(x, branch.length, n.descen, ed=c(ED=NA))
#'
#' @param x numeric. A Named numeric vector of presence-absence
#' @param branch.length numeric. A Named numeric vector of branch length for
#' each species
#' @param n.descen numeric. A Named numeric vector of number of descendants for
#' each branch
#' @param ed numeric. numeric vector with the values of evolutionary distinctiveness
#'
#' @author Gabriela Alves-Ferreira and Neander Marcel Heming
#'
#' @references Isaac, N. J., Turvey, S. T., Collen, B., Waterman, C. and
#' Baillie, J. E. (2007). Mammals on the EDGE: conservation priorities based on
#' threat and phylogeny. PLoS ONE 2, e296.
#'
#' @return numeric
# #' @export
.vec.ed <- function(x, branch.length, n.descen, ed=c(ED=NA)){

  x[is.na(x)] <- 0 # 0 for all value = NA

  if(sum(x) == 0) { # return NA if x = 0

    return(ed)

  } else { # if the sum of x is non-zero then do this:

    pres <- x == 1 # only species present in the vector
    # species <- names(x)
    ed[] <- sum(branch.length[pres]/n.descen[pres]) # evolutionary distinctiveness

  }

  return(ed)

}

#' Calculate Evolutionary distinctiveness for each raster cell
#'
#' @description This function calculates evolutionary distinctiveness according
#' to the fair-proportion index.
#'
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1)
#' for a set of species. The layers (species) must be sorted according to the
#' tree order. See the phylo.pres function.
#' @param branch.length numeric. A Named numeric vector containing the branch
#' length of each specie.
#' @param n.descen numeric. A Named numeric vector of number of descendants for
#' each branch
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster
#' with that many cores is created and used.
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed for fun.
#'
#' @author Gabriela Alves-Ferreira and Neander Marcel Heming
#'
#' @references Isaac, N. J., Turvey, S. T., Collen, B., Waterman, C. and Baillie,
#'  J. E. (2007). Mammals on the EDGE: conservation priorities based on threat
#'  and phylogeny. PLoS ONE 2, e296.
#'
#' @return SpatRaster
#'
# #' @export
# #' @examples
.rast.ed.B <- function(x, branch.length, n.descen, cores = 1, filename = "", ...){

  # evolutionary distinctiveness
  red <- terra::app(x, fun = .vec.ed,
                    branch.length, n.descen, cores = cores,
                    filename = filename, ...)

  terra::set.names(red, "ED") # layer name

  return(red)

}


#' Calculate Evolutionary distinctiveness for raster data
#'
#' @description This function calculates evolutionary distinctiveness according
#' to the fair-proportion index.
#'
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1)
#' for a set of species. The layers (species) must be sorted according to the
#' tree order. See the phylo.pres function.
#' @param branch.length numeric. A Named numeric vector containing the branch
#' length of each specie.
#' @param n.descen numeric. A Named numeric vector of number of descendants for
#' each branch
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster
#' with that many cores is created and used.
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed for fun.
#' @inheritParams geo.phylo
#'
#' @author Gabriela Alves-Ferreira and Neander Marcel Heming
#'
#' @references Isaac, N. J., Turvey, S. T., Collen, B., Waterman, C. and Baillie,
#' J. E. (2007). Mammals on the EDGE: conservation priorities based on threat
#' and phylogeny. PLoS ONE 2, e296.
#'
#' @return SpatRaster
#'
#' @examples
#' library(terra)
#' library(phyloraster)
#' x <- rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' # phylogenetic tree
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
#' data <- phylo.pres(x, tree)
#' ed <- rast.ed(data$x, branch.length = data$branch.length, n.descen = data$n.descen)
#' plot(ed)
#'
#' @export
rast.ed <- function(x, tree,
                    branch.length, n.descen,
                    cores = 1, filename = "", ...){

  ## object checks
  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  ### initial argument check
  {
    miss4 <- arg.check(match.call(), c("LR", "inv.R", "branch.length", "n.descen"))
    miss.tree <- arg.check(match.call(), "tree")

    if(any(miss4) & miss.tree){

      stop("Either argument 'tree' or all 'LR', 'inv.R', 'branch.length', and 'n.descen' need to be supplied")

    } else if(any(miss4)){

      data <- phylo.pres(x, tree)
      # area.branch <- inv.range(data$x, data$branch.length)

      x <- data$x
      # LR <- area.branch$LR
      # inv.R <- area.branch$inv.R
      branch.length <- data$branch.length
      n.descen <- data$n.descendants

    } else if(any(#isFALSE(identical(names(x), names(LR))),
      #isFALSE(identical(names(x), names(inv.R))),
      isFALSE(identical(names(x), names(branch.length))),
      isFALSE(identical(names(x), names(n.descen))))) {

      data <- phylo.pres(x, tree)
      # area.branch <- inv.range(data$x, data$branch.length)

      x <- data$x
      # LR <- area.branch$LR
      # inv.R <- area.branch$inv.R
      branch.length <- data$branch.length
      n.descen <- data$n.descendants
    }

  }


  ## evolutionary distinctiveness
  red <- .rast.ed.B(x,
                    branch.length, n.descen,
                    cores, filename, ...)

  return(red)

}


#' Evolutionary distinctiveness standardized for specie richness
#'
#' @description Calculates the standardized effect size for evolutionary
#' distinctiveness. See Details for more information.
#'
#' @param random character. A character indicating the type of randomization.
#' The currently available randomization methods are "tip", "site", "species" or
#' "both" (site and species).
#' @inheritParams geo.phylo
#' @inheritParams SESraster::SESraster
#'
#' @return SpatRaster
#'
#' @author Neander M. Heming and Gabriela Alves-Ferreira
#'
#' @details The spatial randomization (spat) keeps the richness exact and samples
#'  species presences proportionally to their observed frequency (i.e. number
#'  of occupied pixels). The randomization will not assign values to cells with
#'  nodata. The phylogenetic randomization shuffles taxa names across all taxa
#'  included in phylogeny.
#'
#'
#' @seealso \code{\link{phylo.pres}}, \code{\link{inv.range}},
#' \code{\link{geo.phylo.ses}},
#' \code{\link{rast.ed.ses}}, \code{\link{rast.pd.ses}},
#' \code{\link{rast.we.ses}}, \code{\link{rast.pe.ses}},
#' \code{\link[SESraster]{bootspat_str}}, \code{\link[SESraster]{bootspat_naive}},
#' \code{\link[SESraster]{bootspat_ff}}, \code{\link[SESraster]{SESraster}}
#'
#' @references Isaac, N. J., Turvey, S. T., Collen, B., Waterman, C. and Baillie,
#' J. E. (2007). Mammals on the EDGE: conservation priorities based on threat
#' and phylogeny. PLoS ONE 2, e296.
#' @references Laffan, S. W., Rosauer, D. F., Di Virgilio, G., Miller, J. T.,
#' González‐Orozco, C. E., Knerr, N., ... & Mishler, B. D. (2016). Range‐weighted
#' metrics of species and phylogenetic turnover can better resolve biogeographic
#' transition zones. Methods in Ecology and Evolution, 7(5), 580-588.
#'
#' @examples
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
#' t <- rast.ed.ses(x, tree, aleats = 5, random = "spat")
#' terra::plot(t)
#'
#' @export
rast.ed.ses <- function(x, tree,
                        branch.length, n.descen,
                        spat_alg = "bootspat_str",
                        spat_alg_args = list(rprob = NULL,
                                             rich = NULL,
                                             fr_prob = NULL),
                        random = c("tip", "spat")[2],
                        aleats = 10,
                        cores = 1, filename = "", ...){

  ## object checks
  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  ### initial argument check
  {
    miss4 <- arg.check(match.call(), c("LR", "inv.R", "branch.length", "n.descen"))
    miss.tree <- arg.check(match.call(), "tree")

    if(any(miss4) & miss.tree){

      stop("Either argument 'tree' or all 'LR', 'inv.R', 'branch.length', and 'n.descen' need to be supplied")

    } else if(any(miss4)){

      data <- phylo.pres(x, tree)
      # area.branch <- inv.range(data$x, data$branch.length)

      x <- data$x
      # LR <- area.branch$LR
      # inv.R <- area.branch$inv.R
      branch.length <- data$branch.length
      n.descen <- data$n.descendants

    } else if(any(#isFALSE(identical(names(x), names(LR))),
      #isFALSE(identical(names(x), names(inv.R))),
      isFALSE(identical(names(x), names(branch.length))),
      isFALSE(identical(names(x), names(n.descen))))) {

      data <- phylo.pres(x, tree)
      # area.branch <- inv.range(data$x, data$branch.length)

      x <- data$x
      # LR <- area.branch$LR
      # inv.R <- area.branch$inv.R
      branch.length <- data$branch.length
      n.descen <- data$n.descendants
    }

  }

  require(SESraster)

  ## function arguments
  #    .rast.ed.B(x, branch.length = bl.random, n.descen,
  #                              filename = temp[[i]], cores = cores)
  FUN_args = list(branch.length = branch.length,
                  n.descen = n.descen,
                  # spp_seq = spp_seq,
                  # spp_seqLR = spp_seqLR,
                  # spp_seqINV = spp_seqINV,
                  # resu = resu,
                  cores = cores)


  ## Null model (bootstrap structure)
  if(random == "tip"){

    ed.ses <- SESraster::SESraster(x,
                                   FUN = ".rast.ed.B", FUN_args = FUN_args,
                                   Fa_sample = "branch.length",
                                   Fa_alg = "sample", Fa_alg_args = list(replace=FALSE),
                                   spat_alg = NULL, spat_alg_args = list(),
                                   # spat_alg = spat_alg, spat_alg_args = spat_alg_args,
                                   aleats = aleats,
                                   cores = cores, filename = filename, ...)

  } else if(random == "spat"){

    ed.ses <- SESraster::SESraster(x,
                                   FUN = ".rast.ed.B", FUN_args = FUN_args,
                                   # Fa_sample = "branch.length",
                                   # Fa_alg = "sample", Fa_alg_args = list(replace=FALSE),
                                   # spat_alg = NULL, spat_alg_args = list(),
                                   spat_alg = spat_alg, spat_alg_args = spat_alg_args,
                                   aleats = aleats,
                                   cores = cores, filename = filename, ...)

  }  else {

    stop("Choose a valid randomization method! The methods currently available are: 'tip', 'spat'.")

  }

  return(ed.ses)

}
