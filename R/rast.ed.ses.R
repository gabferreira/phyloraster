#' Calculate Evolutionary distinctiveness for each raster cell
#'
#' @description This function calculates evolutionary distinctiveness according
#' to the fair-proportion index.
#'
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1)
#' for a set of species. The layers (species) must be sorted according to the
#' tree order. See the phylo.pres function.
#' @inheritParams geo.phylo
#'
#' @author Neander Marcel Heming and Gabriela Alves-Ferreira
#'
#' @references Isaac, N. J., Turvey, S. T., Collen, B.,
#' Waterman, C. and Baillie,
#'  J. E. (2007). Mammals on the EDGE: conservation priorities
#'  based on threat
#'  and phylogeny. PLoS ONE 2, e296.
#'
#' @return SpatRaster
#'
#' @keywords internal
.rast.ed.B <- function(x, edge.path, branch.length, n.descen,
                       filename = "", ...){

  # evolutionary distinctiveness
  red <- terra::app(x,
                    function(x, H1, branch.length, n.descen){
                      if(all(is.na(x))) return(NA)
                      if(sum(x, na.rm = T)==0) return(0)

                    sum((crossprod(H1, x)>0) * (branch.length/n.descen))/sum(x)
                    }, H1 = edge.path, branch.length = branch.length,
                    n.descen = n.descen,
                    filename = filename, ...)
  # red <- sum(x*(branch.length/n.descen),
  #            filename = filename, ...)

  terra::set.names(red, "ED") # layer name

  return(red)

}


#' Calculate Evolutionary distinctiveness for raster data
#'
#' @description This function calculates evolutionary distinctiveness according
#' to the fair-proportion index. The values represents the mean ED for species
#' presents in each raster cell.
#'
#' @inheritParams geo.phylo.ses
#' @inheritParams phylo.pres
#'
#' @author Neander Marcel Heming and Gabriela Alves-Ferreira
#'
#' @references Isaac, N. J., Turvey, S. T., Collen, B., Waterman, C. and
#' Baillie,
#' J. E. (2007). Mammals on the EDGE: conservation priorities based on threat
#' and phylogeny. PLoS ONE 2, e296.
#'
#' @return SpatRaster
#'
#' @examples
#' \donttest{
#' library(terra)
#' library(phyloraster)
#' x <- rast(system.file("extdata", "rast.presab.tif",
#' package="phyloraster"))
#' # phylogenetic tree
#' tree <- ape::read.tree(system.file("extdata", "tree.nex",
#' package="phyloraster"))
#' data <- phylo.pres(x[[1:3]], tree)
#' ed <- rast.ed(data$x, edge.path = data$edge.path,
#'               branch.length = data$branch.length,
#'               n.descen = data$n.descen)
#' plot(ed)
#' }
#' @export
rast.ed <- function(x, tree,
                    edge.path, branch.length, n.descen,
                    full_tree_metr = TRUE,
                    filename = "", ...){

  ## object checks
  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  ### initial argument check
  {
    miss4 <- arg.check(match.call(), c("inv.R", "edge.path",
                                       "branch.length", "n.descen"))
    miss.tree <- arg.check(match.call(), "tree")

    if(any(miss4) & miss.tree){

      stop("Either argument 'tree' or all 'inv.R', 'edge.path',
           'branch.length', and 'n.descen' need to be supplied")

    } else if(any(miss4)){

      data <- phylo.pres(x, tree, full_tree_metr = full_tree_metr)
      # area.branch <- inv.range(data$x, data$branch.length)

      x <- data$x
      # LR <- area.branch$LR
      # inv.R <- area.branch$inv.R
      edge.path <- data$edge.path
      branch.length <- data$branch.length
      n.descen <- data$n.descendants

    } else if(any(isFALSE(identical(names(x), rownames(edge.path))) #,
                  #isFALSE(identical(names(x), names(LR))),
                  # isFALSE(identical(names(x), names(branch.length))),
                  # isFALSE(identical(names(x), names(n.descen)))
    )) {

      data <- phylo.pres(x, tree, full_tree_metr = full_tree_metr)
      # area.branch <- inv.range(data$x, data$branch.length)

      x <- data$x
      # LR <- area.branch$LR
      # inv.R <- area.branch$inv.R
      edge.path <- data$edge.path
      branch.length <- data$branch.length
      n.descen <- data$n.descendants
    }

  }
  ## evolutionary distinctiveness
  red <- .rast.ed.B(x,
                    edge.path, branch.length, n.descen,
                    filename = filename, ...)

  return(red)

}


#' Standardized effect size for Evolutionary distinctiveness
#'
#' @description Calculates the standardized effect size for evolutionary
#' distinctiveness. See Details for more information.
#'
#' @inheritParams geo.phylo.ses
#' @inheritParams phylo.pres
#'
#' @seealso \code{\link{phylo.pres}},
#' \code{\link{inv.range}},
#' \code{\link{geo.phylo.ses}},
#' \code{\link{rast.ed.ses}},
#' \code{\link{rast.pd.ses}},
#' \code{\link{rast.we.ses}},
#' \code{\link{rast.pe.ses}},
#' \code{\link[SESraster]{bootspat_str}},
#' \code{\link[SESraster]{bootspat_naive}},
#' \code{\link[SESraster]{bootspat_ff}},
#' \code{\link[SESraster]{SESraster}}
#'
#' @return SpatRaster. The function returns the observed value of the metric,
#' the mean of the simulations calculated over n times, the standard deviation
#' of the simulations, the standardized effect size (SES) for the metric,
#' and the p-values.
#'
#' @details The dependency ‘SESraster’ is used to calculate the null models.
#' This package currently implements six algorithms to randomize binary species
#'  distribution with several levels of constraints:
#'  SIM1, SIM2, SIM3, SIM5, SIM6 and SIM9 (sensu Gotelli 2000).
#'  The methods implemented in ‘SESraster’ are based on how species
#'  (originally rows) and sites (originally columns) are treated
#'  (i.e. fixed, equiprobable, or proportional sums) (Gotelli 2000).
#'  By default, the ‘phyloraster’ uses the function bootspat_ str() from the
#'  ‘SESraster’ package to conduct the randomizations, but the user is free
#'  to choose any of the other methods mentioned above through the spat_alg
#'  argument in the *.ses() functions of the ‘phyloraster’ package.
#'  The bootspat_str() is equivalent to the SIM5 (proportional-fixed) method of
#'  Gotelli (2000), which partially relaxes the spatial structure of species
#'  distributions, but keeps the spatial structure of the observed richness
#'  pattern across cells.
#'
#' @references Gotelli, N. J. 2000.
#' Null model analysis of species co-occurrence patterns.
#' – Ecology 81: 2606–2621.
#' @references Heming, N. M., Mota, F. M. M. and Alves-Ferreira, G. 2023.
#' SESraster: raster randomization for null hypothesis testing.
#' https://CRAN.R-project.org/package=SESraster.
#'
#' @author Neander M. Heming and Gabriela Alves-Ferreira
#'
#' @examples
#' \donttest{
#' library(phyloraster)
#' library(SESraster)
#' x <- terra::rast(system.file("extdata", "rast.presab.tif",
#' package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex",
#' package="phyloraster"))
#' t <- rast.ed.ses(x[[1:10]], tree, aleats = 3)
#' terra::plot(t)
#' }
#' @export
rast.ed.ses <- function(x, tree,
                        edge.path, branch.length, n.descen,
                        full_tree_metr = TRUE,
                        spat_alg = "bootspat_str",
                        spat_alg_args = list(rprob = NULL,
                                             rich = NULL,
                                             fr_prob = NULL),
                        aleats = 10,
                        filename = "", ...){

  requireNamespace("SESraster")
  message("Please cite SESraster when using spatial null models.
          See: citation(SESraster)")


  ## object checks
  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  ### initial argument check
  {
    miss4 <- arg.check(match.call(), c("inv.R", "edge.path", "branch.length",
                                       "n.descen"))
    miss.tree <- arg.check(match.call(), "tree")

    if(any(miss4) & miss.tree){

      stop("Either argument 'tree' or all 'inv.R',
      'edge.path', 'branch.length',
           and 'n.descen' need to be supplied")

    } else if(any(miss4)){

      data <- phylo.pres(x, tree, full_tree_metr = full_tree_metr)
      # area.branch <- inv.range(data$x,
      # data$branch.length)

      x <- data$x
      # LR <- area.branch$LR
      # inv.R <- area.branch$inv.R
      edge.path <- data$edge.path
      branch.length <- data$branch.length
      n.descen <- data$n.descendants

    } else if(any(isFALSE(identical(names(x), rownames(edge.path))) #,
                  # isFALSE(identical(names(x), names(branch.length))),
                  # isFALSE(identical(names(x), names(n.descen)))
    )) {

      data <- phylo.pres(x, tree, full_tree_metr = full_tree_metr)
      # area.branch <- inv.range(data$x, data$branch.length)

      x <- data$x
      # LR <- area.branch$LR
      # inv.R <- area.branch$inv.R
      edge.path <- data$edge.path
      branch.length <- data$branch.length
      n.descen <- data$n.descendants
    }

  }


  ## function arguments
  #    .rast.ed.B(x, branch.length = bl.random, n.descen,
  #                              filename = temp[[i]], cores = cores)
  FUN_args <- list(edge.path = edge.path,
                   branch.length = branch.length,
                   n.descen = n.descen
                   # spp_seq = spp_seq,
                   # spp_seqLR = spp_seqLR,
                   # spp_seqINV = spp_seqINV,
                   # resu = resu,
                   # cores = cores
  )


  ## Null model (bootstrap structure)
  ed.ses <- SESraster::SESraster(x,
                                 FUN = ".rast.ed.B", FUN_args = FUN_args,
                                 # Fa_sample = "branch.length",
                                 # Fa_alg = "sample", Fa_alg_args =
                                 # list(replace=FALSE),
                                 # spat_alg = NULL, spat_alg_args = list(),
                                 spat_alg = spat_alg, spat_alg_args =
                                   spat_alg_args,
                                 aleats = aleats,
                                 filename = filename, ...)

  return(ed.ses)

}
