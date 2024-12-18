#' Calculate weighted endemism for each raster cell
#'
#' @description Calculate the sum of the inverse of the range size for species
#' present in each raster cell.
#'
#' @inheritParams geo.phylo.ses
#'
#' @return SpatRaster
#'
#' @author Neander Marcel Heming and Gabriela Alves-Ferreira
#'
#' @references Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J.,
#' VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and endemism
#' in conservation. In: Systematics and Conservation Evaluation (eds Forey PL,
#' Humphries CJ, Vane-Wright RI), p. 438. Oxford University Press, Oxford.
#' @references Crisp, M., Laffan, S., Linder, H., Monro, A. (2001). Endemism in
#' the Australian flora. Journal of Biogeography, 28, 183–198.
#'
#' @keywords internal
.rast.we.B <- function(x, inv.R, filename = "", ...){

  # weighted endemism
  rend <- sum(x*inv.R,
              filename = filename, ...)

  terra::set.names(rend, "WE") # layer name

  return(rend)
}

#' Calculate weighted endemism for raster data
#'
#' @description Calculate the weighted endemism for species present in raster
#' data.
#'
#' @inheritParams geo.phylo.ses
#'
#' @return SpatRaster
#'
#'
#' @references Laffan, S. W., Rosauer, D. F., Di Virgilio, G., Miller, J. T.,
#' González‐Orozco, C. E., Knerr, N., ... & Mishler, B. D. (2016).
#' Range‐weighted
#' metrics of species and phylogenetic turnover can better resolve
#' biogeographic
#' transition zones. Methods in Ecology and Evolution, 7(5), 580-588.
#' @references Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J.,
#' VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and endemism
#' in
#'  conservation. In: Systematics and Conservation Evaluation (eds Forey PL,
#'  Humphries CJ, Vane-Wright RI), p. 438. Oxford University Press, Oxford.
#' @references Crisp, M., Laffan, S., Linder, H., Monro, A. (2001). Endemism
#' in the Australian flora. Journal of Biogeography, 28, 183–198.
#'
#' @author Neander Marcel Heming and Gabriela Alves Ferreira
#' @examples
#' \donttest{
#' library(terra)
#' library(phyloraster)
#' x <- rast(system.file("extdata", "rast.presab.tif",
#' package="phyloraster"))
#' inv.R <- inv.range(x)
#' we <- rast.we(x, inv.R)
#' plot(we)
#'}
#' @export
rast.we <- function(x, inv.R,
                    filename = "", ...){

  ## object checks
  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  ### initial argument check
  {
    miss4 <- arg.check(match.call(), c("inv.R", "branch.length", "n.descen"))
    # miss.tree <- arg.check(match.call(), "tree")

    # if(any(miss4) & miss.tree){
    #
    #   stop("Argument 'inv.R'need to be supplied")
    #
    # } else
    if(any(miss4)){

      # data <- phylo.pres(x, tree)
      # area.branch <- inv.range(data$x, data$branch.length)
      inv.R <- inv.range(x)

      # x <- data$x
      # LR <- area.branch$LR
      # inv.R <- inv.range(x)$inv.R
      # branch.length <- data$branch.length
      # n.descen <- data$n.descendants

    } else if(any(#isFALSE(identical(names(x), names(LR))),
      isFALSE(identical(names(x), names(inv.R)))
      # isFALSE(identical(names(x), names(branch.length))),
      # isFALSE(identical(names(x), names(n.descen)))
    )) {

      # data <- phylo.pres(x, tree)
      # area.branch <- inv.range(data$x, data$branch.length)
      inv.R <- inv.range(x)

      # x <- data$x
      # LR <- area.branch$LR
      # inv.R <- inv.range(x)$inv.R
      # rs <- range_size(x)
      # branch.length <- data$branch.length
      # n.descen <- data$n.descendants
    }

  }

  ## weighted endemism
  .rast.we.B(x,
             inv.R,
             filename = filename, ...)
}



#' Calculate weighted endemism standardized for species richness
#'
#' @description Calculates the standardized effect size for weighted endemism.
#' See Details for more information.
#'
#' @inheritParams geo.phylo.ses
#'
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
#' @author Gabriela Alves-Ferreira and Neander Heming
#'
#' @examples
#' \donttest{
#' library(terra)
#' library(SESraster)
#' x <- terra::rast(system.file("extdata", "rast.presab.tif",
#' package="phyloraster"))
#' t <- rast.we.ses(x[[1:10]], aleats = 3)
#' plot(t)
#'}
#' @export
rast.we.ses <- function(x,
                        inv.R,
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
    miss4 <- arg.check(match.call(), c("inv.R", "branch.length", "n.descen"))
    # miss.tree <- arg.check(match.call(), "tree")

    # if(any(miss4) & miss.tree){
    #
    #   stop("Argument 'inv.R' need to be supplied")
    #
    # } else if(any(miss4)){
    if(any(miss4)){

      # data <- phylo.pres(x, tree)
      # area.branch <- inv.range(data$x, data$branch.length)

      # x <- data$x
      # LR <- area.branch$LR
      inv.R <- inv.range(x)
      # branch.length <- data$branch.length
      # n.descen <- data$n.descendants

    } else if(any(#isFALSE(identical(names(x), names(LR))),
      # isFALSE(identical(names(x), names(n.descen))),
      # isFALSE(identical(names(x), names(branch.length))),
      isFALSE(identical(names(x), names(inv.R))))) {

      # data <- phylo.pres(x, tree)
      # area.branch <- inv.range(data$x, data$branch.length)

      # x <- data$x
      # LR <- area.branch$LR
      inv.R <- inv.range(x)
      # branch.length <- data$branch.length
      # n.descen <- data$n.descendants
    }

  }


  ## vectorization setup
  # nspp <- terra::nlyr(x)
  # spp_seq <- seq_len(nspp)
  # spp_seqLR <- spp_seq + nspp
  # spp_seqINV <- spp_seq + 2*nspp
  # spp_seqINV <- spp_seq + nspp
  # resu <- setNames(rep(NA, 5), c("SR", "PD", "ED", "PE", "WE"))

  ## function arguments
  FUN_args <- list(
    inv.R = inv.R
    # branch.length = branch.length,
    # n.descen = n.descen,
    # spp_seq = spp_seq,
    # # spp_seqLR = spp_seqLR,
    # spp_seqINV = spp_seqINV,
    # resu = resu,
    # cores = cores
    )


  we.ses <- SESraster::SESraster(x,
                                 FUN = ".rast.we.B", FUN_args = FUN_args,
                                 # Fa_sample = "branch.length",
                                 # Fa_alg = "sample",
                                 # Fa_alg_args = list(replace=FALSE),
                                 # spat_alg = NULL,
                                 # spat_alg_args = list(),
                                 spat_alg = spat_alg,
                                 spat_alg_args = spat_alg_args,
                                 aleats = aleats,
                                 # cores = cores,
                                 filename = filename, ...)

  return(we.ses)

}
