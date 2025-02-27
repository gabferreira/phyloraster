#' Evaluate if the rasters generated in the function fits on available memory
#'
#' @description
#' Tests if the amount of RAM required is available to process a SpatRaster
#'
#' @inheritParams terra::mem_info
#'
#' @return logical
#' @keywords internal
.fit.memory <- function(x, n=1){
  # x rasters will be generated in this function, let's see if there is enough
  # memory in the user's pc
  sink(nullfile())    # suppress output
  mi <- terra::mem_info(x, n)[5] != 0 # proc in memory = T TRUE means that it
  # fits in the pc's memory, so you wouldn't have to use temporary files
  sink()
  return(mi)
}

#' Calculate species richness for raster data
#'
#' @description Calculate the species richness for raster data.
#'
#' @usage rast.sr(x, filename = "", ...)
#'
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1)
#' for a set of species.
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed down from a calling
#' function.
#'
#' @author Gabriela Alves Ferreira and Neander Marcel Heming
#' @return SpatRaster
#' @export
#' @examples
#' \donttest{
#' x <- terra::rast(system.file("extdata", "rast.presab.tif",
#' package="phyloraster"))
#' rse <- phyloraster::rast.sr(x)
#' terra::plot(rse)
#' }
rast.sr <- function(x, filename = "", ...){

  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  # richness
  rsr <- terra::app(x, sum, na.rm = T, filename = filename, ...)

  names(rsr) <- "SR"

  return(rsr)
}

#' Calculate phylogenetic community metrics for raster data
#'
#' Calculate species richness, phylogenetic diversity, evolutionary
#' distinctiveness,
#' phylogenetic endemism and weighted endemism using rasters as input
#'
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1)
#' for a set of species. The layers (species) must be sorted according to the
#' tree order. See the phylo.pres function.
#' @param inv.R SpatRaster. inverse of range size calculated from
#' \code{\link{inv.range}}
#' @param branch.length numeric. A named numerical vector containing the branch
#' length for each species.
#' @param n.descen numeric. A Named numeric vector of number of descendants for
#' each branch
#'
#' @inheritParams geo.phylo
#'
#' @return SpatRaster with one layer for each metric
#'
#' @details Community metrics calculated:
#' \itemize{
##'    \item{Phylogenetic diversity (Faith 1992)}
##'    \item{Species Richness}
##'    \item{Evolutionary distinctiveness by fair-proportion
##'    (Isaac et al. 2007)}
##'    \item{Phylogenetic endemism (Rosauer et al. 2009)}
##'    \item{Weighted endemism (Crisp et al. 2001, Williams et al. 1994)}
##'}
#'
#' @author Neander Marcel Heming
#'
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C.
#' and Cook, L. G. (2009). Phylogenetic endemism: a new approach for identifying
#' geographical concentrations of evolutionary history.
#' Molecular ecology, 18(19),
#' 4061-4072.
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic
#' diversity. Biological conservation, 61(1), 1-10.
#' @references Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J. and
#' VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and endemism in
#' conservation. In: Systematics and Conservation Evaluation (eds Forey PL,
#' Humphries C.J., Vane-Wright RI), p. 438. Oxford University Press, Oxford.
#' @references Crisp, M., Laffan, S., Linder, H. and Monro, A. (2001). Endemism
#' in the Australian flora. Journal of Biogeography, 28, 183–198.
#' @references Isaac, N. J., Turvey, S. T., Collen, B., Waterman, C. and
#' Baillie, J. E. (2007). Mammals on the EDGE: conservation priorities based on
#' threat and phylogeny. PLoS ONE 2, e296.
#' @references Laffan, S. W., Rosauer, D. F., Di Virgilio, G., Miller, J. T.,
#' González‐Orozco, C. E., Knerr, N., ... & Mishler, B. D. (2016).
#' Range‐weighted metrics of species and phylogenetic turnover can better
#' resolve biogeographic transition zones. Methods in Ecology and Evolution,
#' 7(5), 580-588.
#'
#' @keywords internal
.rast.geo.phylo <- function(x,
                            inv.R,
                            edge.path, branch.length, n.descen,
                            # spp_seq, spp_seqrange.BL, spp_seqINV,
                            # resu = stats::setNames(as.double(rep(NA, 5)),
                            #c("SR", "PD", "ED", "PE", "WE")),
                            filename = "", ...){

  mi <- .fit.memory(x, 1) ## proc in memory = TRUE means that it fits in the
  # pc's memory, so you wouldn't have to use temporary files

  # temporary files
  temp <- paste0(tempfile(), 1:5, "g.tif")  # to store the xe raster

  geop <- terra::rast(list(rast.sr(x,
                                   overwrite=TRUE, filename =
                                     ifelse(mi, "", temp[1])), # SR
                           .rast.pd.B(x, edge.path, branch.length,
                                      overwrite=TRUE, filename =
                                        ifelse(mi, "", temp[2])), # PD
                           .rast.ed.B(x, edge.path, branch.length, n.descen,
                                      overwrite=TRUE, filename =
                                        ifelse(mi, "", temp[3])), # ED
                           .rast.pe.B(x, inv.R, branch.length,
                                      metric = "pe",
                                      overwrite=TRUE, filename =
                                        ifelse(mi, "", temp[4])), # PE
                           .rast.we.B(x, inv.R,
                                      overwrite=TRUE, filename =
                                        ifelse(mi, "", temp[5])) # WE
  ))

  if(filename != ""){
    geop <- terra::writeRaster(geop, filename, ...)
  }

  unlink(temp)
  return(geop)
}

#' Calculate phylogenetic community metrics for raster data
#'
#' Calculate species richness, phylogenetic diversity, evolutionary
#' distinctiveness,
#' phylogenetic endemism and weighted endemism using rasters as input.
#'
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1)
#' for a set of species. The layers (species) will be sorted according to the
#' tree order. See the phylo.pres function.
#'
#' @inheritParams phylo.pres
#' @param inv.R SpatRaster. Inverse of range size. See \code{\link{inv.range}}
#' @param edge.path matrix. Matrix representing the paths through the tree from
#' root to each tip. See \code{\link{phylo.pres}}
#' @param branch.length numeric. A Named numeric vector of branch length for
#' each species. See \code{\link{phylo.pres}}
#' @param n.descen numeric. A Named numeric vector of number of descendants for
#' each branch. See \code{\link{phylo.pres}}
#' @inheritParams terra::app
#' @param ... additional arguments passed for terra::app
#'
#' @return SpatRaster with one layer for each metric
#'
#' @details Community metrics calculated:
#' \itemize{
##'    \item{Phylogenetic diversity (Faith 1992)}
##'    \item{Species Richness}
##'    \item{Evolutionary distinctiveness by fair-proportion
##'    (Isaac et al. 2007)}
##'    \item{Phylogenetic endemism (Rosauer et al. 2009)}
##'    \item{Weighted endemism (Crisp et al. 2001, Williams et al. 1994)}
##'}
#'
#' @seealso \code{\link{phylo.pres}}, \code{\link{inv.range}},
#' \code{\link{rast.ed}}, \code{\link{rast.pd}},
#' \code{\link{rast.we}}, \code{\link{rast.pe}}, \code{\link{rast.sr}},
#' \code{\link{geo.phylo.ses}},
#'
#' @author Neander Marcel Heming
#'
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C.
#' and Cook, L. G. (2009). Phylogenetic endemism: a new approach for identifying
#' geographical concentrations of evolutionary history. Molecular ecology,
#' 18(19), 4061-4072.
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic
#' diversity. Biological conservation, 61(1), 1-10.
#' @references Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J.
#' and VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and
#' endemism in conservation. In: Systematics and Conservation Evaluation
#' (eds Forey PL, Humphries C.J., Vane-Wright RI), p. 438. Oxford University
#' Press, Oxford.
#' @references Crisp, M., Laffan, S., Linder, H. and Monro, A. (2001).
#' Endemism in the Australian flora. Journal of Biogeography, 28, 183–198.
#' @references Isaac, N. J., Turvey, S. T., Collen, B., Waterman, C. and
#' Baillie, J. E. (2007). Mammals on the EDGE: conservation priorities based on
#' threat and phylogeny. PLoS ONE 2, e296.
#' @references Laffan, S. W., Rosauer, D. F., Di Virgilio, G., Miller, J. T.,
#' González‐Orozco, C. E., Knerr, N., ... & Mishler, B. D. (2016).
#' Range‐weighted metrics of species and phylogenetic turnover can better
#' resolve biogeographic transition zones. Methods in Ecology and Evolution,
#' 7(5), 580-588.
#'
#' @examples
#' \donttest{
#' library(terra)
#' library(phyloraster)
#' x <- terra::rast(system.file("extdata", "rast.presab.tif",
#' package="phyloraster"))[[1:10]]
#' tree <- ape::read.tree(system.file("extdata", "tree.nex",
#' package="phyloraster"))
#' data <- phylo.pres(x, tree)
#' inv.R <- inv.range(data$x)
#' t <- geo.phylo(data$x, inv.R = inv.R, edge.path = data$edge.path,
#' branch.length = data$branch.length, n.descen = data$n.descendants)
#' terra::plot(t)
#' }
#' @export
geo.phylo <- function(x, tree,
                      inv.R, edge.path, branch.length, n.descen,
                      full_tree_metr = TRUE,
                      filename = "", ...){

  ## object checks
  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  ### initial argument check
  {
    miss4 <- arg.check(match.call(), c("inv.R", "branch.length", "n.descen"))
    miss.tree <- arg.check(match.call(), "tree")

    if(any(miss4) & miss.tree){

      stop("Either argument 'tree' or all 'inv.R',
      'edge.path', 'branch.length',
           and 'n.descen' need to be supplied")

    } else if(any(miss4)){

      data <- phylo.pres(x, tree, full_tree_metr = full_tree_metr)

      x <- data$x
      inv.R <- inv.range(x)
      edge.path <- data$edge.path
      branch.length <- data$branch.length
      n.descen <- data$n.descendants

    } else if(any(isFALSE(identical(names(x), names(inv.R))) #,
                  # isFALSE(identical(names(x), names(branch.length))),
                  # isFALSE(identical(names(x), names(n.descen)))
                  )) {

      data <- phylo.pres(x, tree, full_tree_metr = full_tree_metr)

      x <- data$x
      inv.R <- inv.range(x)
      edge.path <- data$edge.path
      branch.length <- data$branch.length
      n.descen <- data$n.descendants

    }

  }

  ## vectorization setup
  resu <- stats::setNames(rep(NA, 5), c("SR", "PD", "ED", "PE", "WE"))

  ## run function
  .rast.geo.phylo(x,
                  inv.R = inv.R,
                  edge.path = edge.path,
                  branch.length = branch.length,
                  n.descen = n.descen,
                  resu = resu,
                  filename = filename, ...)
}


#' Calculate phylogenetic community metrics and their standardized effect sizes
#' for raster data
#'
#' @description Calculates the standardized effect size for phylogenetic
#' community metrics. See Details for more information.
#'
#' @inheritParams geo.phylo
#' @inheritParams SESraster::SESraster
#' @inheritParams phylo.pres
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
#' @seealso \code{\link{phylo.pres}},
#' \code{\link{inv.range}},
#' \code{\link{geo.phylo}},
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
#' Null model analysis of species co-occurrence patterns. –
#' Ecology 81: 2606–2621.
#' @references Heming, N. M., Mota, F. M. M. and Alves-Ferreira, G. 2023.
#' SESraster: raster randomization for null hypothesis testing.
#' https://CRAN.R-project.org/package=SESraster.
#'
#' @author Neander Marcel Heming
#'
#' @examples
#' \donttest{
#' library(terra)
#' library(phyloraster)
#' require("SESraster")
#' x <- terra::rast(system.file("extdata", "rast.presab.tif",
#' package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex",
#' package="phyloraster"))
#' tses <- geo.phylo.ses(x = x,
#'                      tree = tree,
#'                       spat_alg = "bootspat_str",
#'                       spat_alg_args = list(rprob = NULL,
#'                                            rich = NULL,
#'                                            fr_prob = NULL),
#'                       aleats = 2)
#' terra::plot(tses)
#' }
#' @export
geo.phylo.ses <- function(x, tree,
                          inv.R, edge.path, branch.length, n.descen,
                          full_tree_metr = TRUE,
                          spat_alg = "bootspat_str",
                          spat_alg_args = list(rprob = NULL,
                                               rich = NULL,
                                               fr_prob = NULL),
                          aleats = 10,
                          cores = 1, filename = "", ...){

  ## object checks
  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  ### initial argument check
  {
    miss4 <- arg.check(match.call(), c("inv.R", "branch.length", "n.descen"))
    miss.tree <- arg.check(match.call(), "tree")

    if(any(miss4) & miss.tree){

      stop("Either argument 'tree' or all 'inv.R', 'edge.path', 'branch.length',
           and 'n.descen' need to be supplied")

    } else if(any(miss4)){

      data <- phylo.pres(x, tree, full_tree_metr = full_tree_metr)

      x <- data$x
      inv.R <- inv.range(x)
      edge.path <- data$edge.path
      branch.length <- data$branch.length
      n.descen <- data$n.descendants

    } else if(any(isFALSE(identical(names(x), names(inv.R))) #,
                  # isFALSE(identical(names(x), names(branch.length))),
                  # isFALSE(identical(names(x), names(n.descen)))
                  )) {

      data <- phylo.pres(x, tree, full_tree_metr = full_tree_metr)

      x <- data$x
      inv.R <- inv.range(x)
      edge.path <- data$edge.path
      branch.length <- data$branch.length
      n.descen <- data$n.descendants
    }

  }

  requireNamespace("SESraster")

  ## vectorization setup
  resu <- stats::setNames(rep(NA, 5), c("SR", "PD", "ED", "PE", "WE"))

  ## function arguments
  FUN_args <- list(inv.R = inv.R,
                  edge.path = edge.path,
                  branch.length = branch.length,
                  n.descen = n.descen,
                  resu = resu,
                  cores = cores)

    ses <- SESraster::SESraster(x,
                                FUN = ".rast.geo.phylo",
                                FUN_args = FUN_args,
                                # Fa_sample = "branch.length",
                                # Fa_alg = "sample",
                                # Fa_alg_args = list(replace=FALSE),
                                # spat_alg = NULL, spat_alg_args = list(),
                                spat_alg = spat_alg,
                                spat_alg_args = spat_alg_args,
                                aleats = aleats,
                                cores = cores, filename = filename, ...)
  return(ses)

}
