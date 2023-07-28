#' Calculate weighted endemism for a vector
#'
#' @description This function calculates phylogenetic endemism for a vector
#'
#' @inheritParams .vec.geo.phylo
#'
#' @return numeric
#'
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C.,
#' & Cook, L. G. (2009). Phylogenetic endemism: a new approach for identifying
#' geographical concentrations of evolutionary history. Molecular ecology,
#' 18(19), 4061-4072.
#'
#' @author Neander Marcel Heming and Gabriela Alves-Ferreira
#'
# #' @export
.vec.we <- function(x, resu = as.double(NA)){

  if(all(is.na(x))){
    return(resu)
  }

  resu[] <- sum(x, na.rm = TRUE)

}

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
# #' @export
# #' @examples
.rast.we.B <- function(x, inv.R,  cores = 1, filename = "", ...){

  # weighted endemism
  rend <- sum(x*inv.R,
              filename = filename, ...)

  terra::set.names(rend, "WE") # layer name

  return(rend)
}

#' Calculate weighted endemism for raster data
#'
#' @description Calculate the weighted endemism for species present in raster data.
#'
#' @inheritParams geo.phylo.ses
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
#' we <- rast.we(x, inv.R)
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
             cores = cores, filename = filename, ...)
}



#' Calculate weighted endemism standardized for species richness
#'
#' @description Calculates the standardized effect size for weighted endemism.
#' See Details for more information.
#'
#' @inheritParams geo.phylo.ses
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
                        # random = c("spat"),
                        aleats = 10,
                        cores = 1, filename = "", ...){

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

  requireNamespace("SESraster")

  ## vectorization setup
  # nspp <- terra::nlyr(x)
  # spp_seq <- seq_len(nspp)
  # spp_seqLR <- spp_seq + nspp
  # spp_seqINV <- spp_seq + 2*nspp
  # spp_seqINV <- spp_seq + nspp
  # resu <- setNames(rep(NA, 5), c("SR", "PD", "ED", "PE", "WE"))

  ## function arguments
  FUN_args = list(
    inv.R = inv.R,
    # branch.length = branch.length,
    # n.descen = n.descen,
    # spp_seq = spp_seq,
    # # spp_seqLR = spp_seqLR,
    # spp_seqINV = spp_seqINV,
    # resu = resu,
    cores = cores)


  we.ses <- SESraster::SESraster(x,
                                 FUN = ".rast.we.B", FUN_args = FUN_args,
                                 # Fa_sample = "branch.length",
                                 # Fa_alg = "sample", Fa_alg_args = list(replace=FALSE),
                                 # spat_alg = NULL, spat_alg_args = list(),
                                 spat_alg = spat_alg, spat_alg_args = spat_alg_args,
                                 aleats = aleats,
                                 cores = cores, filename = filename, ...)

  return(we.ses)

}
