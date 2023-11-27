#' Calculate phylogenetic diversity for each raster cell
#'
#' @description Calculate the sum of the branch length for species present in
#' raster data.
#'
#' @inheritParams geo.phylo.ses
#'
#' @return SpatRaster
#'
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
#' pd <- rast.pd(data$x, edge.path = data$edge.path,
#' branch.length = data$branch.length)
#' plot(pd)
#' }
#' @export
rast.pd <- function(x, tree,
                    edge.path, branch.length,
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

      data <- phylo.pres(x, tree)
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

      data <- phylo.pres(x, tree)
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


#' Phylogenetic diversity standardized for species richness
#'
#' @description Calculates the standardized effect size for phylogenetic
#' diversity.
#'  See Details for more information.
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
#' @seealso \code{\link{phylo.pres}}, \code{\link{inv.range}},
#' \code{\link{geo.phylo.ses}},
#' \code{\link{rast.ed.ses}}, \code{\link{rast.pd.ses}},
#' \code{\link{rast.we.ses}}, \code{\link{rast.pe.ses}},
#' \code{\link[SESraster]{bootspat_str}},
#' \code{\link[SESraster]{bootspat_naive}},
#' \code{\link[SESraster]{bootspat_ff}},
#' \code{\link[SESraster]{SESraster}}
#'
#' @return SpatRaster. The function returns the observed phylogenetic diversity,
#' the mean of the simulations calculated over n times, the standard deviation
#' of
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
#'                 branch.length = data$branch.length, aleats = 3,
#'                 random = "spat")
#' plot(t)
#' }
#' @export
rast.pd.ses <- function(x, tree,
                        edge.path, branch.length,
                        spat_alg = "bootspat_str",
                        spat_alg_args = list(rprob = NULL,
                                             rich = NULL,
                                             fr_prob = NULL),
                        random = c("tip", "spat")[2],
                        aleats = 10,
                        # cores = 1,
                        filename = "", ...){

  requireNamespace("SESraster")

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

      data <- phylo.pres(x, tree)
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

      data <- phylo.pres(x, tree)
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
  if(random == "tip"){

    pd.ses <- SESraster::SESraster(x,
                                   FUN = ".rast.pd.B", FUN_args = FUN_args,
                                   Fa_sample = "branch.length",
                                   Fa_alg = "sample", Fa_alg_args =
                                   list(replace=FALSE),
                                   spat_alg = NULL, spat_alg_args = list(),
                                   # spat_alg = spat_alg, spat_alg_args =
                                   # spat_alg_args,
                                   aleats = aleats,
                                   filename = filename, ...)

  } else if(random == "spat"){

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

  } else {

    stop("Choose a valid randomization method! The methods currently available
    are: 'tip' and 'spat'.")

  }

  return(pd.ses)

}
