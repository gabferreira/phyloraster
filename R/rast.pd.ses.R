#' Calculate phylogenetic diversity for each raster cell
#'
#' @description Calculate the sum of the branch length for species present in
#' raster data.
#'
#' @inheritParams geo.phylo.ses
#'
#' @return SpatRaster
#'
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic
#' diversity. Biological conservation, 61(1), 1-10.
#'
#' @author Neander Marcel Heming and Gabriela Alves-Ferreira
#'
#' @keywords internal
.rast.pd.B <- function(x, edge.path, branch.length, filename = "", ...){

  # phylogenetic diversity
  rpd <- terra::app(x,
                    function(x, H1, branch.length){
                      if(all(is.na(x))) return(NA)

                      sum((crossprod(H1, x)>0) * branch.length)
                    }, H1 = edge.path, branch.length = branch.length,
                    filename = filename, ...)
  # rpd <- sum(x*branch.length,
  #            filename = filename, ...)

  terra::set.names(rpd, "PD")

  return(rpd)
}


#' Calculate phylogenetic diversity for raster data
#'
#' @description Calculate the sum of the branch length for species present in
#' each cell of the raster.
#'
#' @inheritParams geo.phylo.ses
#' @inheritParams phylo.pres
#'
#' @return SpatRaster
#'
#' @references Faith, D. P. (1992). Conservation evaluation and
#' phylogenetic diversity. Biological conservation, 61(1), 1-10.
#'
#' @author Neander Marcel Heming and Gabriela Alves-Ferreira
#'
#' @examples
#' \donttest{
#' library(terra)
#' library(phyloraster)
#' x <- rast(system.file("extdata", "rast.presab.tif",
#' package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex",
#' package="phyloraster"))
#' data <- phylo.pres(x[[1:3]], tree)
#' pd <- rast.pd(data$x, data$tree)
#' plot(pd)
#' }
#' @export
rast.pd <- function(x, tree,
                    edge.path, branch.length,
                    full_tree_metr = TRUE,
                    filename = "", ...){

  ## object checks
  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the
         calculations.")
  }
  ### initial argument check
  {
    miss4 <- arg.check(match.call(), c("inv.R", "edge.path",
                                       "branch.length",
                                       "n.descen"))
    miss.tree <- arg.check(match.call(), "tree")

    if(any(miss4) & miss.tree){

      stop("Either argument 'tree' or 'edge.path' and 'branch.length'
           need to be supplied")

    } else if(any(miss4)){

      data <- phylo.pres(x, tree, full_tree_metr = full_tree_metr)
      # area.branch <- inv.range(data$x, data$branch.length)

      x <- data$x
      # LR <- area.branch$LR
      # inv.R <- area.branch$inv.R
      edge.path <- data$edge.path
      branch.length <- data$branch.length
      # n.descen <- data$n.descendants

    } else if(any(#isFALSE(identical(names(x), names(branch.length))),
      # isFALSE(identical(names(x), names(inv.R))),
      isFALSE(identical(names(x), rownames(edge.path)))
      # isFALSE(identical(names(x), names(n.descen)))
    )) {

      data <- phylo.pres(x, tree, full_tree_metr = full_tree_metr)
      # area.branch <- inv.range(data$x, data$branch.length)

      x <- data$x
      # LR <- area.branch$LR
      # inv.R <- area.branch$inv.R
      edge.path <- data$edge.path
      branch.length <- data$branch.length
      # n.descen <- data$n.descendants
    }

  }


  ## phylogenetic diversity
  rpd <- .rast.pd.B(x,
                    edge.path, branch.length,
                    filename = filename, ...)

  return(rpd)

}


#' Standardized effect size for Phylogenetic diversity
#'
#' @description Calculates the standardized effect size for phylogenetic
#' diversity.
#'  See Details for more information.
#'
#' @inheritParams geo.phylo.ses
#' @inheritParams phylo.pres
#'
#' @seealso \code{\link{phylo.pres}}, \code{\link{inv.range}},
#' \code{\link{geo.phylo.ses}},
#' \code{\link{rast.ed.ses}}, \code{\link{rast.pd.ses}},
#' \code{\link{rast.we.ses}}, \code{\link{rast.pe.ses}},
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
#' @author Gabriela Alves-Ferreira and Neander Heming
#'
#' @examples
#' \donttest{
#' library(terra)
#' library(phyloraster)
#' library(SESraster)
#' x <- rast(system.file("extdata", "rast.presab.tif",
#' package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex",
#' package="phyloraster"))
#' data <- phylo.pres(x[[1:10]], tree)
#' t <- rast.pd.ses(data$x, edge.path = data$edge.path,
#'                 branch.length = data$branch.length, aleats = 3)
#' plot(t)
#' }
#' @export
rast.pd.ses <- function(x, tree,
                        edge.path, branch.length,
                        full_tree_metr = TRUE,
                        spat_alg = "bootspat_str",
                        spat_alg_args = list(rprob = NULL,
                                             rich = NULL,
                                             fr_prob = NULL),
                        aleats = 10,
                        # cores = 1,
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

      # stop("Either argument 'tree' or all 'inv.R', 'branch.length', and
      # 'n.descen' need to be supplied")
      stop("Either argument 'tree' or both 'edge.path' and 'branch.length'
      need to be supplied")

    } else if(any(miss4)){

      data <- phylo.pres(x, tree, full_tree_metr = full_tree_metr)
      # area.branch <- inv.range(data$x, data$branch.length)

      x <- data$x
      # LR <- area.branch$LR
      # inv.R <- area.branch$inv.R
      edge.path <- data$edge.path
      branch.length <- data$branch.length
      # n.descen <- data$n.descendants

    } else if(any(isFALSE(identical(names(x), rownames(edge.path)))
                  # isFALSE(identical(names(x), names(inv.R))),
                  # isFALSE(identical(names(x), names(branch.length)))
                  # isFALSE(identical(names(x), names(n.descen)))
    )) {

      data <- phylo.pres(x, tree, full_tree_metr = full_tree_metr)
      # area.branch <- inv.range(data$x, data$branch.length)

      x <- data$x
      # LR <- area.branch$LR
      # inv.R <- area.branch$inv.R
      edge.path <- data$edge.path
      branch.length <- data$branch.length
      # n.descen <- data$n.descendants
    }

  }


  ## function arguments
  #    .rast.ed.B(x, branch.length = bl.random, n.descen,
  #                              filename = temp[[i]], cores = cores)
  FUN_args <- list(edge.path = edge.path,
                   branch.length = branch.length
                   # n.descen = n.descen
  )


  ## Null model (bootstrap structure)
  pd.ses <- SESraster::SESraster(x,
                                 FUN = ".rast.pd.B", FUN_args = FUN_args,
                                 # Fa_sample = "branch.length",
                                 # Fa_alg = "sample", Fa_alg_args =
                                 #list(replace=FALSE),
                                 # spat_alg = NULL, spat_alg_args = list(),
                                 spat_alg = spat_alg, spat_alg_args =
                                   spat_alg_args,
                                 aleats = aleats,
                                 # cores = cores,
                                 filename = filename, ...)

  return(pd.ses)

}
