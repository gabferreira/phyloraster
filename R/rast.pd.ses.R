#' Calculate phylogenetic diversity (Faith 1992) for a vector
#'
#' @description This function calculates the sum of the branch length for a set
#'  of species for one sample.
#'
#' @param x numeric. A named numerical vector of presence-absence for one sample.
#' @param branch.length numeric. A named numerical vector containing the branch
#' length for each species.
#' @param pd numeric. A numeric vector with the values of phylogenetic diversity.
#'
#' @return numeric
#'
#'
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic
#' diversity. Biological conservation, 61(1), 1-10.
#'
#' @author Neander Marcel Heming and Gabriela Alves-Ferreira
#'
# #' @export
# #' @examples
.vec.pd <- function(x, branch.length, pd = c(PD=NA)){

  if(all(is.na(x))){
    return(pd)
  }

  pd[] <- sum(x*branch.length, na.rm = TRUE) # pd Faith 1992

  return(pd)
}


#' Calculate phylogenetic diversity for each raster cell
#'
#' @description Calculate the sum of the branch length for species present in
#' raster data.
#'
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1)
#' for a set of species. The layers (species) must be sorted according to the tree
#' order. See the phylo.pres function.
#' @param branch.length numeric. A Named numerical vector containing the branch
#' length for a set of species.
#' @param filename character. Output filename.
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster with
#' that many cores is created and used.
#' @param ... additional arguments to be passed passed down from a calling function.
#'
#' @return SpatRaster
#'
#'
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic
#' diversity. Biological conservation, 61(1), 1-10.
#'
#' @author Neander Marcel Heming and Gabriela Alves-Ferreira
#'
# #' @export
# #' @examples
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
#' @description Calculate the sum of the branch length for species present in
#' each cell of the raster.
#'
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1)
#' for a set of species. The layers (species) must be sorted according to the
#' tree order. See the phylo.pres function.
#' @param branch.length numeric. A Named numerical vector containing the branch
#'  length for a set of species.
#' @param filename character. Output filename.
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster
#' with that many cores is created and used.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @inheritParams geo.phylo
#'
#' @return SpatRaster
#'
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological conservation, 61(1), 1-10.
#'
#' @author Neander Marcel Heming and Gabriela Alves-Ferreira
#'
#' @examples
#' library(terra)
#' library(phyloraster)
#' x <- rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
#' data <- phylo.pres(x, tree)
#' pd <- rast.pd(data$x, branch.length = data$branch.length)
#' plot(pd)
#'
#' @export
rast.pd <- function(x, tree,
                    branch.length,
                    cores = 1, filename = "", ...){

  ## object checks
  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  ### initial argument check
  {
    miss4 <- arg.check(match.call(), c("LR", "inv.R",
                                       "branch.length", "n.descen"))
    miss.tree <- arg.check(match.call(), "tree")

    if(any(miss4) & miss.tree){

      stop("Either argument 'tree' or 'branch.length' need to be supplied")

    } else if(any(miss4)){

      data <- phylo.pres(x, tree)
      # area.branch <- inv.range(data$x, data$branch.length)

      x <- data$x
      # LR <- area.branch$LR
      # inv.R <- area.branch$inv.R
      branch.length <- data$branch.length
      # n.descen <- data$n.descendants

    } else if(any(isFALSE(identical(names(x), names(branch.length)))
                  # isFALSE(identical(names(x), names(inv.R))),
                  # isFALSE(identical(names(x), names(LR))),
                  # isFALSE(identical(names(x), names(n.descen)))
    )) {

      data <- phylo.pres(x, tree)
      # area.branch <- inv.range(data$x, data$branch.length)

      x <- data$x
      # LR <- area.branch$LR
      # inv.R <- area.branch$inv.R
      branch.length <- data$branch.length
      # n.descen <- data$n.descendants
    }

  }


  ## phylogenetic diversity
  rpd <- .rast.pd.B(x, branch.length = branch.length,
                    cores = cores, filename = filename, ...)

  return(rpd)

}


#' Phylogenetic diversity standardized for species richness
#'
#' @description Calculates the standardized effect size for phylogenetic diversity.
#'  See Details for more information.
#'
#' @param random character. A character indicating the type of randomization.
#' The currently available randomization methods are "tip", "site", "species" or
#' "both" (site and species).
#' @inheritParams geo.phylo
#' @inheritParams SESraster::SESraster
#'
#' @return SpatRaster
#'
#' @details The spatial randomization (spat) keeps the richness exact and samples
#'  species presences proportionally to their observed frequency (i.e. number
#'  of occupied pixels). The randomization will not assign values to cells with
#'  nodata. The phylogenetic randomization shuffles taxa names across all taxa
#'  included in phylogeny.
#'
#' @seealso \code{\link{phylo.pres}}, \code{\link{inv.range}},
#' \code{\link{geo.phylo.ses}},
#' \code{\link{rast.ed.ses}}, \code{\link{rast.pd.ses}},
#' \code{\link{rast.we.ses}}, \code{\link{rast.pe.ses}},
#' \code{\link[SESraster]{bootspat_str}}, \code{\link[SESraster]{bootspat_naive}},
#' \code{\link[SESraster]{bootspat_ff}}, \code{\link[SESraster]{SESraster}}
#'
#' @return SpatRaster. The function returns the observed phylogenetic diversity,
#' the mean of the simulations calculated over n times, the standard deviation of
#' the simulations, and the standardized effect size (SES) for the phylogenetic
#' diversity.
#'
#' @details The spatial randomization (spat) keeps the richness exact and
#' samples species presences proportionally to their observed frequency
#' (i.e. number of occupied pixels). The randomization will not assign values
#' to cells with nodata. The phylogenetic randomization shuffles taxa names
#' across all taxa included in phylogeny.
#'
#' @author Gabriela Alves-Ferreira and Neander Heming
#'
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic
#' diversity. Biological conservation, 61(1), 1-10.
#'
#' @examples
#' library(terra)
#' library(phyloraster)
#' x <- rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
#' data <- phylo.pres(x, tree)
#' t <- rast.pd.ses(data$x, branch.length = data$branch.length, aleats = 10, random = "spat")
#' plot(t)
#'
#' @export
rast.pd.ses <- function(x, tree,
                        branch.length,
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

      # stop("Either argument 'tree' or all 'LR', 'inv.R', 'branch.length', and 'n.descen' need to be supplied")
      stop("Either argument 'tree' or 'branch.length' need to be supplied")

    } else if(any(miss4)){

      data <- phylo.pres(x, tree)
      # area.branch <- inv.range(data$x, data$branch.length)

      x <- data$x
      # LR <- area.branch$LR
      # inv.R <- area.branch$inv.R
      branch.length <- data$branch.length
      # n.descen <- data$n.descendants

    } else if(any(#isFALSE(identical(names(x), names(LR))),
      #isFALSE(identical(names(x), names(inv.R))),
      isFALSE(identical(names(x), names(branch.length)))
      # isFALSE(identical(names(x), names(n.descen)))
    )) {

      data <- phylo.pres(x, tree)
      # area.branch <- inv.range(data$x, data$branch.length)

      x <- data$x
      # LR <- area.branch$LR
      # inv.R <- area.branch$inv.R
      branch.length <- data$branch.length
      # n.descen <- data$n.descendants
    }

  }

  require(SESraster)

  ## function arguments
  #    .rast.ed.B(x, branch.length = bl.random, n.descen,
  #                              filename = temp[[i]], cores = cores)
  FUN_args = list(branch.length = branch.length,
                  # n.descen = n.descen,
                  # spp_seq = spp_seq,
                  # spp_seqLR = spp_seqLR,
                  # spp_seqINV = spp_seqINV,
                  # resu = resu,
                  cores = cores)


  ## Null model (bootstrap structure)
  if(random == "tip"){

    pd.ses <- SESraster::SESraster(x,
                                   FUN = ".rast.pd.B", FUN_args = FUN_args,
                                   Fa_sample = "branch.length",
                                   Fa_alg = "sample", Fa_alg_args = list(replace=FALSE),
                                   spat_alg = NULL, spat_alg_args = list(),
                                   # spat_alg = spat_alg, spat_alg_args = spat_alg_args,
                                   aleats = aleats,
                                   cores = cores, filename = filename, ...)

  } else if(random == "spat"){

    pd.ses <- SESraster::SESraster(x,
                                   FUN = ".rast.pd.B", FUN_args = FUN_args,
                                   # Fa_sample = "branch.length",
                                   # Fa_alg = "sample", Fa_alg_args = list(replace=FALSE),
                                   # spat_alg = NULL, spat_alg_args = list(),
                                   spat_alg = spat_alg, spat_alg_args = spat_alg_args,
                                   aleats = aleats,
                                   cores = cores, filename = filename, ...)

  } else {

    stop("Choose a valid randomization method! The methods currently available are: 'tip', 'spat'.")

  }

  return(pd.ses)

}
