#' Evaluate if the rasters generated in the function fits on available memory
#'
#' @description
#' Tests if the amount of RAM that is required is available to process a SpatRaster
#'
#' @inheritParams terra::mem_info
#'
#' @return logical
# #' @export
# #' @examples
.fit.memory <- function(x, n=1){
  # x rasters will be generated in this function, let's see if there is enough memory in the user's pc
  sink(nullfile())    # suppress output
  mi <- terra::mem_info(x, n)[5] != 0 # proc in memory = T TRUE means that it fits in the pc's memory, so you wouldn't have to use temporary files
  sink()
  return(mi)
}

#' Calculate species richness for raster data
#'
#' @description Calculate the species richness for raster data.
#' @usage rast.se(x, filename = NULL, cores = 1, ...)
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for a set of species.
#' @param filename character. Output filename.
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster with that many cores is created and used.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @author Gabriela Alves Ferreira and Neander Marcel Heming
#' @return SpatRaster
#' @export
#' @examples
#' \dontrun{
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' rse <- rast.se(x)
#' terra::plot(rse)
#' }
rast.sr <- function(x, filename = "", cores = 1, ...){

  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  # richness
  rsr <- terra::set.names(terra::app(x, sum, na.rm = TRUE, cores = cores,
                                     filename = filename, ...),
                          "SR")

  return(rsr)
}


#' Calculate phylogenetic community metrics for raster data
#'
#' Calculate species richness, phylogenetic diversity, evolutionary distinctiveness,
#' phylogenetic endemism and weighted endemism using rasters as input and output.
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1)
#' for a set of species. The layers (species) must be sorted according to the
#' tree order. See the phylo.pres function.
#' @param area.branch description
#' @param inv.R description
#' @param branch.length description
#' @param n.descen description
#' @inheritParams terra::app
#' @param ... additional arguments passed for terra::app
#'
#' @return SpatRaster with one layer for each metric
#'
#' @details Community metrics calculated:
#' \itemize{
##'    \item{Phylogenetic diversity (Faith 1992)}
##'    \item{Richness}
##'    \item{Evolutionary distinctiveness by fair-proportion (Isaac et al. 2007)}
##'    \item{Phylogenetic endemism (Rosauer et al. 2009)}
##'    \item{Weighted endemism (Crisp et al. 2001, Williams et al. 1994)}
##'}
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C. and Cook, L. G. (2009). Phylogenetic endemism: a new approach for identifying geographical concentrations of evolutionary history. Molecular ecology, 18(19), 4061-4072.
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological conservation, 61(1), 1-10.
#' @references Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J. and VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and endemism in conservation. In: Systematics and Conservation Evaluation (eds Forey PL, Humphries CJ, Vane-Wright RI), p. 438. Oxford University Press, Oxford.
#' @references Crisp, M., Laffan, S., Linder, H. and Monro, A. (2001). Endemism in theAustralian flora. Journal of Biogeography, 28, 183–198.
#' @references Isaac, N. J., Turvey, S. T., Collen, B., Waterman, C. and Baillie, J. E. (2007). Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE 2, e296.
#' @references Laffan, S. W., Rosauer, D. F., Di Virgilio, G., Miller, J. T., González‐Orozco, C. E., Knerr, N., ... & Mishler, B. D. (2016). Range‐weighted metrics of species and phylogenetic turnover can better resolve biogeographic transition zones. Methods in Ecology and Evolution, 7(5), 580-588.
#'
.rast.geo.phylo <- function(x, branch.length, n.descen,
                            spp_seq, spp_seqLR, spp_seqINV,
                            resu = setNames(rep(NA, 5), c("SR", "PD", "ED", "PE", "WE")),
                            cores = 1, filename = "", ...){

  terra::app(x3,
             function(x, branch.length, n.descen,
                      spp_seq, spp_seqLR, spp_seqINV,
                      resu){
               if(all(is.na(x))){
                 return(resu)
               }
               ### separating raster layers of each cell
               # xa <- x[spp_seq]
               # xLR <- x[spp_seqLR]
               # xINV <- x[spp_seqINV]

               ### computing metrics
               se <- sum(x[spp_seq], na.rm = T)
               pd <- .vec.pd(x[spp_seq], branch.length)[[1]] # [[1]] return only the first raster
               ed <- .evol.distin(x[spp_seq], branch.length, n.descen)[[1]] # [[1]] return only the first raster
               pe <- .vec.wpe(x, spp_seq, spp_seqLR)
               we <- .vec.wpe(x, spp_seq, spp_seqINV)

               resu[] <- c(se, pd, ed, pe, we)
               return(resu)
             },
             branch.length = branch.length,
             n.descen = n.descen,
             spp_seq = spp_seq,
             spp_seqLR = spp_seqLR,
             spp_seqINV = spp_seqINV,
             resu = resu,
             cores = cores, filename = filename, ...)

}

#' Calculate phylogenetic community metrics for raster data
#'
#' Calculate species richness, phylogenetic diversity, evolutionary distinctiveness,
#' phylogenetic endemism and weighted endemism using rasters as input and output.
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1)
#' for a set of species. The layers (species) must be sorted according to the
#' tree order. See the phylo.pres function.
#' @param area.branch description
#' @param inv.R description
#' @param branch.length description
#' @param n.descen description
#' @inheritParams terra::app
#' @param ... additional arguments passed for terra::app
#'
#' @return SpatRaster with one layer for each metric
#'
#' @details Community metrics calculated:
#' \itemize{
##'    \item{Phylogenetic diversity (Faith 1992)}
##'    \item{Richness}
##'    \item{Evolutionary distinctiveness by fair-proportion (Isaac et al. 2007)}
##'    \item{Phylogenetic endemism (Rosauer et al. 2009)}
##'    \item{Weighted endemism (Crisp et al. 2001, Williams et al. 1994)}
##'}
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C. and Cook, L. G. (2009). Phylogenetic endemism: a new approach for identifying geographical concentrations of evolutionary history. Molecular ecology, 18(19), 4061-4072.
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological conservation, 61(1), 1-10.
#' @references Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J. and VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and endemism in conservation. In: Systematics and Conservation Evaluation (eds Forey PL, Humphries CJ, Vane-Wright RI), p. 438. Oxford University Press, Oxford.
#' @references Crisp, M., Laffan, S., Linder, H. and Monro, A. (2001). Endemism in theAustralian flora. Journal of Biogeography, 28, 183–198.
#' @references Isaac, N. J., Turvey, S. T., Collen, B., Waterman, C. and Baillie, J. E. (2007). Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE 2, e296.
#' @references Laffan, S. W., Rosauer, D. F., Di Virgilio, G., Miller, J. T., González‐Orozco, C. E., Knerr, N., ... & Mishler, B. D. (2016). Range‐weighted metrics of species and phylogenetic turnover can better resolve biogeographic transition zones. Methods in Ecology and Evolution, 7(5), 580-588.
#'
#' @examples
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
#' data <- phylo.pres(x, tree)
#' area.branch <- phyloraster::inv.range(data$x, data$branch.length, LR=T)
#' area.branch <- inv.range(data$x, data$branch.length)
#' t <- geo.phylo(x=data$x, LR=area.branch$LR, inv.R=area.branch$inv.R,
#'                branch.length=data$branch.length, n.descen=data$n.descendants)
#' terra::plot(t)
#'
#' @export
geo.phylo <- function(x, LR, inv.R,
                      branch.length, n.descen,
                      cores = 1, filename = "", ...){

  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  nspp <- terra::nlyr(x)
  spp_seq <- seq_len(nspp)
  spp_seqLR <- spp_seq + nspp
  spp_seqINV <- spp_seq + 2*nspp
  resu <- setNames(rep(NA, 5), c("SR", "PD", "ED", "PE", "WE"))
  x3 <- c(x, LR, inv.R)

  .rast.geo.phylo(x3, branch.length, n.descen,
                  spp_seq, spp_seqLR, spp_seqINV,
                  resu = setNames(rep(NA, 5), c("SR", "PD", "ED", "PE", "WE")),
                  cores = 1, filename = "", ...)
}
