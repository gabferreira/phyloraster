#' Calculate weighted endemism for each raster cell
#'
#' @description Calculate the sum of the inverse of the range size for species
#' present in each raster cell.
#'
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1)
#' for a set of species.
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster with
#'  that many cores is created and used.
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @inheritParams .vec.geo.phylo
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
# #' @export
# #' @examples
.rast.we.B <- function(x, spp_seq, spp_seqINV, cores = 1, filename = "", ...){

  # weighted endemism
  rend <- terra::app(x,
                     .vec.wpe,
                     spp_seq, spp_seqINV,
                     cores = cores, filename = filename, ...)

  terra::set.names(rend, "WE") # layer name

  return(rend)
}


#' Calculate weighted endemism for raster data
#'
#' @description Calculate the weighted endemism for species present in raster data.
#'
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1)
#' for a set of species.
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster with
#'  that many cores is created and used.
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed down from a calling function.
#'
#' @inheritParams geo.phylo
#'
#' @return SpatRaster
#'
#'
#' @references Laffan, S. W., Rosauer, D. F., Di Virgilio, G., Miller, J. T.,
#' González‐Orozco, C. E., Knerr, N., ... & Mishler, B. D. (2016). Range‐weighted
#' metrics of species and phylogenetic turnover can better resolve biogeographic
#' transition zones. Methods in Ecology and Evolution, 7(5), 580-588.
#' @references Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J.,
#' VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and endemism in
#'  conservation. In: Systematics and Conservation Evaluation (eds Forey PL,
#'  Humphries CJ, Vane-Wright RI), p. 438. Oxford University Press, Oxford.
#' @references Crisp, M., Laffan, S., Linder, H., Monro, A. (2001). Endemism
#' in the Australian flora. Journal of Biogeography, 28, 183–198.
#'
#' @author Neander Marcel Heming and Gabriela Alves Ferreira
#'
#' @examples
#' library(terra)
#' library(phyloraster)
#' x <- rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' inv.R <- inv.range(x)
#' we <- rast.we(x, inv.R$inv.R)
#' plot(we)
#'
#' @export
rast.we <- function(x, inv.R,
                    cores = 1, filename = "", ...){

  ## object checks
  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  ### initial argument check
  {
    miss4 <- arg.check(match.call(), c("LR", "inv.R", "branch.length", "n.descen"))
    # miss.tree <- arg.check(match.call(), "tree")

    # if(any(miss4) & miss.tree){
    #
    #   stop("Argument 'inv.R'need to be supplied")
    #
    # } else
    if(any(miss4)){

      # data <- phylo.pres(x, tree)
      # area.branch <- inv.range(data$x, data$branch.length)
      # inv.R <- inv.range(x)

      # x <- data$x
      # LR <- area.branch$LR
      inv.R <- inv.range(x)$inv.R
      # branch.length <- data$branch.length
      # n.descen <- data$n.descendants

    } else if(any(#isFALSE(identical(names(x), names(LR))),
      isFALSE(identical(names(x), names(inv.R)))
      # isFALSE(identical(names(x), names(branch.length))),
      # isFALSE(identical(names(x), names(n.descen)))
    )) {

      # data <- phylo.pres(x, tree)
      # area.branch <- inv.range(data$x, data$branch.length)
      # inv.R <- inv.range(x)

      # x <- data$x
      # LR <- area.branch$LR
      inv.R <- inv.range(x)$inv.R
      # branch.length <- data$branch.length
      # n.descen <- data$n.descendants
    }

  }

  ## vectorization setup
  nspp <- terra::nlyr(x)
  spp_seq <- seq_len(nspp)
  # spp_seqLR <- spp_seq + nspp
  spp_seqINV <- spp_seq + nspp
  # resu <- setNames(rep(NA, 5), c("SR", "PD", "ED", "PE", "WE"))

  # .rast.we.B(x, spp_seq, spp_seqINV, cores = 1, filename = "", ...){

  ## weighted endemism
  .rast.we.B(c(x, inv.R),
             spp_seq, spp_seqINV,
             cores = cores, filename = filename, ...)
}



#' Calculate weighted endemism standardized for species richness
#'
#' @description Calculates the standardized effect size for weighted endemism.
#' See Details for more information.
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
#' @references Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J.,
#' VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and endemism
#' in conservation. In: Systematics and Conservation Evaluation (eds Forey PL,
#' Humphries CJ, Vane-Wright RI), p. 438. Oxford University Press, Oxford.
#' @references Crisp, M., Laffan, S., Linder, H., Monro, A. (2001). Endemism in the
#' Australian flora. Journal of Biogeography, 28, 183–198.
#'
#' @author Neander Marcel Heming and Gabriela Alves-Ferreira
#'
#' @examples
#' library(terra)
#' library(phyloraster)
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' t <- rast.we.ses(x, aleats = 10, random = "spat")
#' plot(t)
#'
#' @export
rast.we.ses <- function(x,
                        inv.R,
                        spat_alg = "bootspat_str",
                        spat_alg_args = list(rprob = NULL,
                                             rich = NULL,
                                             fr_prob = NULL),
                        random = c("spat"),
                        aleats = 10,
                        cores = 1, filename = "", ...){

  ## object checks
  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  ### initial argument check
  {
    miss4 <- arg.check(match.call(), c("LR", "inv.R", "branch.length", "n.descen"))
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
      inv.R <- inv.range(x)$inv.R
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
      inv.R <- inv.range(x)$inv.R
      # branch.length <- data$branch.length
      # n.descen <- data$n.descendants
    }

  }

  require(SESraster)

  ## vectorization setup
  nspp <- terra::nlyr(x)
  spp_seq <- seq_len(nspp)
  # spp_seqLR <- spp_seq + nspp
  # spp_seqINV <- spp_seq + 2*nspp
  spp_seqINV <- spp_seq + nspp
  # resu <- setNames(rep(NA, 5), c("SR", "PD", "ED", "PE", "WE"))

  ## function arguments
  FUN_args = list(
    # branch.length = branch.length,
    # n.descen = n.descen,
    spp_seq = spp_seq,
    # spp_seqLR = spp_seqLR,
    spp_seqINV = spp_seqINV,
    # resu = resu,
    cores = cores)

  # if(random == "tip"){
  #
  #   we.ses <- SESraster::SESraster(c(x, inv.R = inv.R),
  #                                  FUN = ".rast.we.B", FUN_args = FUN_args,
  #                                  Fa_sample = "branch.length",
  #                                  Fa_alg = "sample", Fa_alg_args = list(replace = FALSE),
  #                                  spat_alg = NULL, spat_alg_args = list(),
  #                                  # spat_alg = spat_alg, spat_alg_args = spat_alg_args,
  #                                  aleats = aleats,
  #                                  cores = cores, filename = filename, ...)
  #
  # } else
  if(random == "spat"){

    we.ses <- SESraster::SESraster(c(x, inv.R = inv.R),
                                   FUN = ".rast.we.B", FUN_args = FUN_args,
                                   # Fa_sample = "branch.length",
                                   # Fa_alg = "sample", Fa_alg_args = list(replace=FALSE),
                                   # spat_alg = NULL, spat_alg_args = list(),
                                   spat_alg = spat_alg, spat_alg_args = spat_alg_args,
                                   aleats = aleats,
                                   cores = cores, filename = filename, ...)
  }  else {

    stop("Choose a valid randomization method! The methods currently available are: 'tip', 'spat'.")

  }
  return(we.ses)

}
