#' Calculate phylogenetic endemism for a raster
#'
#' Calculate phylogenetic endemism using rasters as input and output.
#'
#' @inheritParams geo.phylo.ses
#' @param branch.length.alt numeric. Branch length calculated by using an
#' alternative phylogeny with non-zero branch lengths converted to a constant
#' value (1) and rescaled so the sum of all branch lengths is 1.
#' @param metric character. Names of the biodiversity metrics to calculate.
#' Available options are: "pe", "pe.alt", "rpe", or "all". See details.
#' @param filename	 character. Output filename
#' @param overwrite	 logical. If TRUE, filename is overwritten
#'
#' @return SpatRaster
#'
#' @export
#' @details Metrics available are:
#' - pe: Phylogenetic endemism (Rosauer et al., 2009)
#' - pe.alt: Alternate Phylogenetic endemism  (Mishler et al., 2014)
#' - rpe: Relative Phylogenetic endemism (Mishler et al., 2014)
#' - all: Calculate all available metrics
#' Alternate phylogenetic endemism (PE.alt, Mishler et al., 2014) is
#' calculated using an alternate phylogeny with non-zero branch lengths
#' converted to a constant value (here we use 1) and rescaled so the sum of all
#' branch lengths is 1.
#' Relative phylogenetic endemism (RPE, Mishler et al., 2014) is the ratio
#'  of phylogenetic endemism (PE, Rosauer et al., 2009) measured on the
#'  original tree versus PE measured on a alternate tree (PE.alt).
#'
#' @references Mishler, B. D., Knerr, N., González-Orozco, C. E., Thornhill,
#' A. H., Laffan, S. W. and Miller, J. T. 2014. Phylogenetic measures of
#' biodiversity and neo- and paleo-endemism in Australian Acacia. –
#' Nat. Commun. 5: 4473.
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C.,
#'  & Cook, L. G. (2009). Phylogenetic endemism: a new approach for identifying
#'  geographical concentrations of evolutionary history. Molecular ecology,
#'  18(19), 4061-4072.
#'
#' @author Gabriela Alves-Ferreira and Neander Heming
#'
#' @keywords internal
.rast.pe.B <- function(x, inv.R, branch.length, branch.length.alt,
                       metric = c("pe", "pe.alt", "rpe", "all")[1],
                       filename = "", overwrite = TRUE, ...){

  if(metric == "pe"){
    # phylogenetic endemism
    resu <- sum(x*inv.R*branch.length,
                filename = filename, overwrite = overwrite, ...)
    terra::set.names(resu, "PE")

  } else if(metric == "pe.alt"){
    # phylogenetic endemism altered
    resu <- sum(x*inv.R*branch.length.alt,
                filename = filename, overwrite = overwrite, ...)
    terra::set.names(resu, "PE.alt")

  } else if(metric == "rpe"){

    # phylogenetic endemism
    pe <- sum(x*inv.R*branch.length,
                filename = filename, overwrite = overwrite, ...)
    # phylogenetic endemism altered
    pe.alt <- sum(x*inv.R*branch.length.alt,
                filename = filename, overwrite = overwrite, ...)

    # relative phylogenetic endemism
    resu <- pe / pe.alt
    terra::set.names(resu, "RPE")

  } else if(metric == "all"){
    # phylogenetic endemism
    pe <- sum(x*inv.R*branch.length,
              filename = filename, overwrite = overwrite, ...)

    # phylogenetic endemism altered
    pe_alt <- sum(x*inv.R*branch.length.alt,
                  filename = filename, overwrite = overwrite, ...)

    # relative phylogenetic endemism
    rpe <- pe / pe_alt

    resu <- c(pe, pe_alt, rpe)
    terra::set.names(resu, c("PE", "PE.alt", "RPE"))
  }
  return(resu)
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
#' Range‐weighted metrics of species and phylogenetic turnover can better
#' resolve biogeographic transition zones. Methods in Ecology and Evolution,
#' 7(5), 580-588.
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C.
#'  and Cook, L. G. (2009). Phylogenetic endemism: a new approach for
#'  identifying geographical concentrations of evolutionary history.
#'  Molecular ecology, 18(19), 4061-4072.
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
                    full_tree_metr = TRUE,
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
             metric = "pe",
             filename = filename, ...)
}

#' Standardized effect size for Phylogenetic endemism
#'
#' @description Calculates the standardized effect size for phylogenetic
#' endemism.
#' See Details for more information.
#'
#' @inheritParams geo.phylo.ses
#' @inheritParams phylo.pres
#' @param branch.length.alt numeric. Branch length calculated by using an
#' alternative phylogeny with non-zero branch lengths converted to a constant
#' value (1) and rescaled so the sum of all branch lengths is 1.
#' @param metric 	character. Names of biodiversity metrics to calculate
#'  (pe, pe_alt, rpe, all). See details.
#' @param filename	 character. Output filename
#' @param overwrite	 logical. If TRUE, filename is overwritten
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
#'  Biodiversity metrics available are:
#' - pe: Phylogenetic endemism (Rosauer et al., 2009)
#' - pe.alt: Alternate Phylogenetic endemism  (Mishler et al., 2014)
#' - rpe: Relative Phylogenetic endemism (Mishler et al., 2014)
#' - all: Calculate all available metrics
#' Alternate phylogenetic endemism (PE.alt, Mishler et al., 2014) is
#' calculated using an alternate phylogeny with non-zero branch lengths
#' converted to a constant value (here we use 1) and rescaled so the sum of all
#' branch lengths is 1.
#' Relative phylogenetic endemism (RPE, Mishler et al., 2014) is the ratio
#'  of phylogenetic endemism (PE, Rosauer et al., 2009) measured on the
#'  original tree versus PE measured on a alternate tree (PE.alt).
#'
#'
#' @references Gotelli, N. J. 2000.
#' Null model analysis of species co-occurrence patterns.
#' Ecology 81: 2606–2621.
#' @references Heming, N. M., Mota, F. M. M. and Alves-Ferreira, G.
#' 2023. SESraster: raster randomization for null hypothesis testing.
#' https://CRAN.R-project.org/package=SESraster.
#' @references Mishler, B. D., Knerr, N., González-Orozco, C. E., Thornhill,
#' A. H., Laffan, S. W. and Miller, J. T. 2014. Phylogenetic measures of
#' biodiversity and neo- and paleo-endemism in Australian Acacia. –
#' Nat. Commun. 5: 4473.
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C.,
#'  & Cook, L. G. (2009). Phylogenetic endemism: a new approach for identifying
#'  geographical concentrations of evolutionary history. Molecular ecology,
#'  18(19), 4061-4072.
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
#' t <- rast.pe.ses(x = data$x, data$tree, aleats = 99, metric = "all")
#' plot(t)
#'}
#' @export
rast.pe.ses <- function(x, tree,
                        branch.length, branch.length.alt,
                        inv.R,
                        full_tree_metr = TRUE,
                        # rs, cellSz,
                        spat_alg = "bootspat_str",
                        spat_alg_args = list(rprob = NULL,
                                             rich = NULL,
                                             fr_prob = NULL),
                        metric = c("pe","pe.alt", "rpe", "all")[4],
                        aleats = 10,
                        cores = 1, filename = "",
                        overwrite = TRUE,
                        ...){

  requireNamespace("SESraster")
  message("Please cite SESraster when using spatial null models.
          See: citation(SESraster)")

  ## object checks
  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  ### initial argument check
  {
    miss4 <- arg.check(match.call(), c("inv.R", "branch.length",
                                       "branch.length.alt",
                                       "n.descen"))
    miss.tree <- arg.check(match.call(), "tree")

    if(any(miss4) & miss.tree){

      stop("Either argument 'tree' or all 'inv.R', 'branch.length', and
      'branch.length.alt' need to be supplied")

    } else if(any(miss4)){

      data <- phylo.pres(x, tree, full_tree_metr = full_tree_metr)
      # area.branch <- inv.range(data$x, data$branch.length)

      x <- data$x
      # range.BL <- area.branch$range.Bl
      inv.R <- inv.range(data$x)
      branch.length <- data$branch.length
      branch.length.alt <- data$branch.length.alt
      # n.descen <- data$n.descendants

    } else if(any(isFALSE(identical(names(x), names(inv.R))),
                  isFALSE(identical(names(x), names(branch.length)))
                  # isFALSE(identical(names(x), names(n.descen)))
    )) {

      data <- phylo.pres(x, tree, full_tree_metr = full_tree_metr)
      # area.branch <- inv.range(data$x, data$branch.length)

      x <- data$x
      # range.BL <- area.branch$range.BL
      # inv.R <- area.branch$inv.R
      inv.R <- inv.range(data$x)
      branch.length <- data$branch.length
      branch.length.alt <- data$branch.length.alt
      # n.descen <- data$n.descendants
    }

  }

  # metric
  if(metric == "pe"){

    ## function arguments

    # .rast.pe.B(xeSZ*x, branch.length, branch.length.alt , metric,
    #            cores = cores, filename = filename)

    FUN_args <- list(
      branch.length = branch.length,
      branch.length.alt = branch.length.alt,
      inv.R = inv.R,
      metric = "pe"
      # n.descen = n.descen,
      # spp_seq = spp_seq,
      # spp_seqrange.BL = spp_seqrange.BL,
      # spp_seqINV = spp_seqINV,
      # resu = resu,
      # cores = cores
    )

    ## Null model (bootstrap structure)
    ses <- SESraster::SESraster(x,
                                FUN = ".rast.pe.B", FUN_args = FUN_args,
                                # Fa_sample = "branch.length",
                                # Fa_alg = "sample",
                                # Fa_alg_args = list(replace=FALSE),
                                # spat_alg = NULL, spat_alg_args = list(),
                                spat_alg = spat_alg,
                                spat_alg_args = spat_alg_args,
                                aleats = aleats,
                                cores = cores, filename = filename,
                                overwrite = overwrite, ...)

    # masking to avoid values outside the extent of the original raster
    ses_masked <- terra::mask(ses, ses$Observed.PE)

    names(ses_masked) <- c("Observed.PE", "Null.Mean.PE",
                           "Null.SD.PE", "SES.PE",
                           "p.lower.PE", "p.upper.PE")

  } else if(metric == "pe.alt"){

    ## function arguments

    # .rast.pe.B(xeSZ*x, branch.length, branch.length.alt , metric,
    #            cores = cores, filename = filename)

    FUN_args <- list(
      branch.length = branch.length,
      branch.length.alt = branch.length.alt,
      inv.R = inv.R,
      metric = "pe.alt"
      # n.descen = n.descen,
      # spp_seq = spp_seq,
      # spp_seqrange.BL = spp_seqrange.BL,
      # spp_seqINV = spp_seqINV,
      # resu = resu,
      # cores = cores
    )

    ## Null model (bootstrap structure)
    ses <- SESraster::SESraster(x,
                                FUN = ".rast.pe.B", FUN_args = FUN_args,
                                # Fa_sample = "branch.length",
                                # Fa_alg = "sample",
                                # Fa_alg_args = list(replace=FALSE),
                                # spat_alg = NULL, spat_alg_args = list(),
                                spat_alg = spat_alg,
                                spat_alg_args = spat_alg_args,
                                aleats = aleats,
                                cores = cores, filename = filename,
                                overwrite = overwrite, ...)

    # masking to avoid values outside the extent of the original raster
    ses_masked <- terra::mask(ses, ses$Observed.PE.alt)

    names(ses_masked) <- c("Observed.PE.alt", "Null.Mean.PE.alt",
                           "Null.SD.PE.alt", "SES.PE.alt",
                           "p.lower.PE.alt", "p.upper.PE.alt")


  } else if(metric == "rpe"){

    ## function arguments

    # .rast.pe.B(xeSZ*x, branch.length, branch.length.alt , metric,
    #            cores = cores, filename = filename)

    FUN_args <- list(
      branch.length = branch.length,
      branch.length.alt = branch.length.alt,
      inv.R = inv.R,
      metric = "rpe"
      # n.descen = n.descen,
      # spp_seq = spp_seq,
      # spp_seqrange.BL = spp_seqrange.BL,
      # spp_seqINV = spp_seqINV,
      # resu = resu,
      # cores = cores
    )

    ## Null model (bootstrap structure)
    ses <- SESraster::SESraster(x,
                                FUN = ".rast.pe.B", FUN_args = FUN_args,
                                # Fa_sample = "branch.length",
                                # Fa_alg = "sample",
                                # Fa_alg_args = list(replace=FALSE),
                                # spat_alg = NULL, spat_alg_args = list(),
                                spat_alg = spat_alg,
                                spat_alg_args = spat_alg_args,
                                aleats = aleats,
                                cores = cores, filename = filename,
                                overwrite = overwrite, ...)

    # masking to avoid values outside the extent of the original raster
    ses_masked <- terra::mask(ses, ses$Observed.RPE)

    names(ses_masked) <- c("Observed.RPE", "Null.Mean.RPE", "Null.SD.RPE",
                           "SES.RPE",  "p.lower.RPE", "p.upper.RPE")

  } else if(metric == "all"){

    ## function arguments

    # .rast.pe.B(xeSZ*x, branch.length, branch.length.alt , metric,
    #            cores = cores, filename = filename)

    FUN_args <- list(
      branch.length = branch.length,
      branch.length.alt = branch.length.alt,
      inv.R = inv.R,
      metric = "all"
      # n.descen = n.descen,
      # spp_seq = spp_seq,
      # spp_seqrange.BL = spp_seqrange.BL,
      # spp_seqINV = spp_seqINV,
      # resu = resu,
      # cores = cores
    )

    ## Null model (bootstrap structure)
    ses <- SESraster::SESraster(x,
                                FUN = ".rast.pe.B", FUN_args = FUN_args,
                                # Fa_sample = "branch.length",
                                # Fa_alg = "sample",
                                # Fa_alg_args = list(replace=FALSE),
                                # spat_alg = NULL, spat_alg_args = list(),
                                spat_alg = spat_alg,
                                spat_alg_args = spat_alg_args,
                                aleats = aleats,
                                cores = cores, filename = filename,
                                overwrite = overwrite, ...)

    # masking to avoid values outside the extent of the original raster
    ses_masked <- terra::mask(ses, ses$Observed.PE)

    names(ses_masked) <- c("Observed.PE", "Null.Mean.PE", "Null.SD.PE",
                           "SES.PE", "p.lower.PE" , "p.upper.PE",
                           "Observed.PE.alt", "Null.Mean.PE.alt",
                           "Null.SD.PE.alt", "SES.PE.alt",
                           "p.lower.PE.alt", "p.upper.PE.alt",
                           "Observed.RPE", "Null.Mean.RPE", "Null.SD.RPE",
                           "SES.RPE",  "p.lower.RPE", "p.upper.RPE")

  }

  return(ses_masked)

}
