#' Calculate phylogenetic endemism for a raster
#'
#' Calculate phylogenetic endemism using rasters as input and output.
#'
#' @inheritParams geo.phylo.ses
#'
#' @return SpatRaster
#'
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C.,
#'  & Cook, L. G. (2009). Phylogenetic endemism: a new approach for identifying
#'  geographical concentrations of evolutionary history. Molecular ecology,
#'  18(19), 4061-4072.
#'
#' @keywords internal
.rast.pe.B <- function(x, inv.R, branch.length, filename = "", ...){

  # phylogenetic endemism
  rpe <- sum(x*inv.R*branch.length,
             filename = filename, ...)

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
#' @description Calculate the sum of the inverse of the range size multiplied
#' by the branch length for the species present in raster data.
#'
#' @inheritParams geo.phylo.ses
#' @inheritParams phylo.pres
#'
#' @author Gabriela Alves-Ferreira and Neander Marcel Heming
#'
#' @references Laffan, S. W., Rosauer, D. F., Di Virgilio, G., Miller, J. T.,
#' González‐Orozco, C. E., Knerr, N., ... & Mishler, B. D. (2016).
#' Range‐weighted
#' metrics of species and phylogenetic turnover can better resolve
#' biogeographic
#' transition zones. Methods in Ecology and Evolution, 7(5), 580-588.
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D.,
#' Donnellan, S. C.
#'  and Cook, L. G. (2009). Phylogenetic endemism: a new approach for
#'  identifying
#'   geographical concentrations of evolutionary history. Molecular ecology,
#'   18(19), 4061-4072.
#'
#' @return SpatRaster
#'
#' @examples
#' \donttest{
#' library(terra)
#' library(phyloraster)
#' x <- rast(system.file("extdata", "rast.presab.tif",
#' package = "phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex",
#' package = "phyloraster"))
#' pe <- rast.pe(x = x[[1:3]], tree)
#' plot(pe)
#'}
#' @export
rast.pe <- function(x, tree,
                    inv.R,
                    branch.length,
                    full_tree_metr = FALSE,
                    filename = "", ...){

  ## object checks
  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed
         for the calculations.")
  }

  ### initial argument check
  {
    miss4 <- arg.check(match.call(), c("rs", "branch.length", "cellSz",
                                       "inv.R", "branch.length", "n.descen"))
    miss.tree <- arg.check(match.call(), "tree")

    if(any(miss4) & miss.tree){

      stop("Either argument 'tree' or 'inv.R' need to be supplied")

    } else if(any(miss4)){

      data <- phylo.pres(x, tree, full_tree_metr = full_tree_metr)
      # area.branch <- inv.range(data$x, data$branch.length)

      x <- data$x
      # range.BL <- area.branch$range.BL
      # inv.R <- area.branch$inv.R
      inv.R <- inv.range(x)
      branch.length <- data$branch.length
      # n.descen <- data$n.descendants
      # rs <- range_size(x)
      # cellSz <- terra::cellSize(x)

    } else if(any(isFALSE(identical(names(x), names(inv.R))),
                  isFALSE(identical(names(x), names(branch.length)))
                  # isFALSE(identical(names(x), names(n.descen)))
    )) {

      data <- phylo.pres(x, tree, full_tree_metr = full_tree_metr)
      # area.branch <- inv.range(data$x, data$branch.length)

      x <- data$x
      # range.BL <- area.branch$range.BL
      # inv.R <- area.branch$inv.R
      inv.R <- inv.range(x)
      branch.length <- data$branch.length
      # n.descen <- data$n.descendants
      # rs <- range_size(x)
      # cellSz <- terra::cellSize(x)
    }

  }

  # plot(.rast.pe.B(xeSZ*x, branch.length))

  ## run function
  .rast.pe.B(x,
             inv.R, branch.length,
             filename = filename, ...)
}


#' Phylogenetic endemism standardized for specie richness
#'
#' @description Calculates the standardized effect size for phylogenetic
#' endemism.
#' See Details for more information.
#'
#' @inheritParams geo.phylo.ses
#'
#' @return SpatRaster
#'
#' @details The spatial randomization (spat) keeps the richness exact and
#' samples
#'  species presences proportionally to their observed frequency (i.e. number
#'  of occupied pixels). The randomization will not assign values to cells with
#'  nodata. The phylogenetic randomization shuffles taxa names across all taxa
#'  included in phylogeny.
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
#'
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C.,
#'  & Cook, L. G. (2009). Phylogenetic endemism: a new approach for identifying
#'   geographical concentrations of evolutionary history.
#'   Molecular ecology, 18(19),
#'   4061-4072.
#'
#' @author Gabriela Alves-Ferreira and Neander Heming
#'
#' @examples
#' \donttest{
#' library(terra)
#' library(phyloraster)
#' library(SESraster)
#' x <- terra::rast(system.file("extdata", "rast.presab.tif",
#' package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex",
#' package="phyloraster"))
#' data <- phylo.pres(x[[1:3]], tree)
#' range.BL <- inv.range(data$x)
#' t <- rast.pe.ses(x = data$x,  tree, aleats = 3,
#' random = "spat")
#' plot(t)
#'}
#' @export
rast.pe.ses <- function(x, tree,
                        branch.length, inv.R,
                        # rs, cellSz,
                        spat_alg = "bootspat_str",
                        spat_alg_args = list(rprob = NULL,
                                             rich = NULL,
                                             fr_prob = NULL),
                        random = c("tip", "spat")[2],
                        aleats = 10,
                        cores = 1, filename = "", ...){

  requireNamespace("SESraster")

  ## object checks
  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  ### initial argument check
  {
    miss4 <- arg.check(match.call(), c("inv.R", "branch.length",
                                       "n.descen"))
    miss.tree <- arg.check(match.call(), "tree")

    if(any(miss4) & miss.tree){

      stop("Either argument 'tree' or all 'inv.R', and 'branch.length'
           need to be supplied")

    } else if(any(miss4)){

      data <- phylo.pres(x, tree)
      # area.branch <- inv.range(data$x, data$branch.length)

      x <- data$x
      # range.BL <- area.branch$range.Bl
      inv.R <- inv.range(data$x)
      branch.length <- data$branch.length
      # n.descen <- data$n.descendants

    } else if(any(isFALSE(identical(names(x), names(inv.R))),
                  isFALSE(identical(names(x), names(branch.length)))
                  # isFALSE(identical(names(x), names(n.descen)))
    )) {

      data <- phylo.pres(x, tree)
      # area.branch <- inv.range(data$x, data$branch.length)

      x <- data$x
      # range.BL <- area.branch$range.BL
      # inv.R <- area.branch$inv.R
      inv.R <- inv.range(data$x)
      branch.length <- data$branch.length
      # n.descen <- data$n.descendants
    }

  }


  ## function arguments

  # .rast.pe.B(xeSZ*x, branch.length,
  #            cores = cores, filename = filename)

  FUN_args <- list(
    branch.length = branch.length,
    inv.R = inv.R
    # n.descen = n.descen,
    # spp_seq = spp_seq,
    # spp_seqrange.BL = spp_seqrange.BL,
    # spp_seqINV = spp_seqINV,
    # resu = resu,
    # cores = cores
    )

  if(random == "tip"){

    ses <- SESraster::SESraster(x,
                                FUN = ".rast.pe.B", FUN_args = FUN_args,
                                Fa_sample = "branch.length",
                                Fa_alg = "sample", Fa_alg_args =
                                  list(replace = FALSE),
                                spat_alg = NULL, spat_alg_args = list(),
                                # spat_alg = spat_alg,
                                # spat_alg_args = spat_alg_args,
                                aleats = aleats,
                                cores = cores, filename = filename, ...)

  } else if(random == "spat"){

    ses <- SESraster::SESraster(x,
                                FUN = ".rast.pe.B", FUN_args = FUN_args,
                                # Fa_sample = "branch.length",
                                # Fa_alg = "sample",
                                # Fa_alg_args = list(replace=FALSE),
                                # spat_alg = NULL, spat_alg_args = list(),
                                spat_alg = spat_alg,
                                spat_alg_args = spat_alg_args,
                                aleats = aleats,
                                cores = cores, filename = filename, ...)
  }  else {

    stop("Choose a valid randomization method! The methods currently available
         are: 'tip', 'spat'.")

  }
  return(ses)

}

