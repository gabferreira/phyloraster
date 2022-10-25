#' Calculate phylogenetic diversity (Faith 1992) for a vector
#'
#' @param x numeric. A Named numeric vector of presence-absence
#' @param branch.length numeric. A Named numeric vector of branch length for each specie
#'
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological conservation, 61(1), 1-10.
#' @return numeric
#' @export
#' @examples

.vec.pd <- function(x, branch.length){
  x[is.na(x)] <- 0

  if(sum(x)== 0) {
    return(c(NA,NA))
  }

  if(sum(x) != 0) {
    pres <- x == 1 # only species present
    species <- names(branch.length)
    n.species <- length(species)
    pd <- sum(branch.length[pres]) # pd Faith 1992
    pd1 <- pd # tive que adicionar pq a funcao nao funcionava se eu retornasse apenas um raster
    # pdr <- pd/sum(pres) # pd relative to richness
  }
  return(c(pd = pd, pd1 = pd1))
}

#' Calculate Evolutionary distinctiveness for a vector
#'
#' This function calculates evolutionary distinctiveness according to the fair-proportion index (Isaac et al., 2007).
#'
#' @param x numeric. A Named numeric vector of presence-absence
#' @param branch.length branch.length numeric. A Named numeric vector of branch length for each specie
#' @param n.descen numeric. A Named numeric vector of number of descendants for each branch
#' @references Isaac, N. J., Turvey, S. T., Collen, B., Waterman, C. & Baillie, J. E. (2007). Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE 2, e296.
#' @return numeric
#'
#' @return
#' @export
#'
#' @examples

.evol.distin <- function(x, branch.length, n.descen){

  x[is.na(x)] <- 0

  if(sum(x) == 0) {
    return(c(NA,NA))
  }

  if(sum(x) != 0){
    pres <- x == 1
    species <- names(x)
    ed <- sum(branch.length[pres]/n.descen[pres])
    ed1 <- ed # tive que adicionar pq a funcao nao funcionava se eu retornasse apenas um raster
  }
  return(c(ed = ed, ed1 = ed1))
}

#' Calculate phylogenetic diversity (PD. Faith, 1992) for a raster
#'
#' Calculate phylogenetic diversity using rasters as input and output. This function follows Faith (1992).
#'
#' @inheritParams rast.we
#' @param x SpatRaster. A presence-absence SpatRaster with the layers ordered according to the tree order.
#' @param branch.length numeric. A Named numeric vector containing the branch length of each specie.
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @return SpatRaster
#' @export
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological conservation, 61(1), 1-10.
#' @examples
#' \dontrun{
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
#' data <- phylo.pres(ras, tree)
#' rast.pd(data$x, data$branch.length)
#' }

rast.pd <- function(x, branch.length, filename = NULL, cores = 1, ...){

  if(!all.equal(names(x), names(branch.length))){
    stop("Species names are not in the same order on 'x' and 'branch.length' arguments! See 'phylogrid::phylo.pres' function.")
  } else {
    # phylogenetic diversity and richness-relative phylogenetic diversity
    if(cores>1){
      rpd <- terra::app(x, fun = .vec.pd,
                        branch.length = branch.length, cores = cores)
      rpd <- rpd[[1]] # select only the first raster
      names(rpd) <- c("PD")
    } else{
      rpd <- terra::app(x, fun = .vec.pd,
                        branch.length = branch.length)
      rpd <- rpd[[1]] # select only the first raster
      names(rpd) <- c("PD")
    }
  }

  if (!is.null(filename)){ # to save the rasters when the output filename is provide
    rpd <- terra::writeRaster(rpd, filename, ...)
  }

  return(rpd)
}


#' Calculate Evolutionary distinctiveness for a raster
#'
#' This function calculates evolutionary distinctiveness according to the fair-proportion index (Isaac et al., 2007).
#'
#' @inheritParams rast.pd
#' @param x SpatRaster. A presence-absence SpatRaster with the layers ordered according to the tree order.
#' @param branch.length numeric. A Named numeric vector containing the branch length of each specie.
#' @param n.descen numeric. A Named numeric vector of number of descendants for each branch
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed for fun.
#' @references Isaac, N. J., Turvey, S. T., Collen, B., Waterman, C. & Baillie, J. E. (2007). Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE 2, e296.
#'
#' @return SpatRaster
#' @export
#'
#' @examples
#' \dontrun{
#' # raster
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' # phylogenetic tree
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
#' data <- phylogrid::phylo.pres(x, tree)
#' ed <- phylogrid::rast.ed(data$x, data$branch.length, data$n.descen, cores = 2)
#' }
#'

rast.ed <- function(x, branch.length, n.descen, filename = NULL, cores = 1, ...){

  if(!all.equal(names(x), names(branch.length))){
    stop("Species names are not in the same order on 'x' and 'branch.length' arguments! See 'phylogrid::phylo.pres' function.")
  }

  # Evolutionary distinctiveness

  if(cores>1){
    red <- terra::app(x, fun = .evol.distin,
                      branch.length, n.descen, cores = 1)
    red <- red[[1]] # separando so pro primeiro raster devido o problema da funcao app de nao calcular quando eh pra retornar apenas um raster
    names(red) <- c("ED")

  } else {
    red <- terra::app(x, fun = .evol.distin,
                      branch.length, n.descen)
    red <- red[[1]] # separando so pro primeiro raster devido o problema da funcao app de nao calcular quando eh pra retornar apenas um raster
    names(red) <- "ED"
  }

  if (!is.null(filename)){ # to save the rasters when the output filename is provide
    red <- terra::writeRaster(red, filename, ...)
  }

  return(red)
}

#' Calculate phylogenetic endemism (PE. Rosauer et al. 2009) for a raster
#'
#' Calculate phylogenetic endemism following Rosauer et al. (2009) using rasters as input and output.
#'
#' @inheritParams rast.we
#' @param x SpatRaster. A presence-absence SpatRaster with the layers ordered according to the tree order
#' @param branch.length numeric. A Named numeric vector containing the branch length of each specie
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @return SpatRaster
#' @export
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C., & Cook, L. G. (2009). Phylogenetic endemism: a new approach for identifying geographical concentrations of evolutionary history. Molecular ecology, 18(19), 4061-4072.
#' @examples
#' \dontrun{
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
#' data <- phylo.pres(ras, tree)
#' rast.pe(data$x, data$branch.length)
#' }

rast.pe <- function(x, branch.length, filename = NULL, cores = 1, ...){

  if(!all.equal(names(x), names(branch.length))){
    stop("Species names are not in the same order on 'x' and 'branch.length' arguments! See 'phylogrid::phylo.pres' function.")
  } else{
    area.branch <- phylogrid::inv.range(x, branch.length)
    # phylogenetic endemism
    # message("Calculating the phylogenetic endemism")
    if(cores>1){
      rpe <- terra::app(area.branch$LR,
                        function(x){
                          if(all(is.na(x))){
                            return(NA)
                          }
                          sum(x, na.rm = T)
                        }, cores = cores)
      rpe <- terra::app(rpe, function(x, m){ # to reescale the values from 0 to 1
        (x/m)
      }, m = terra::minmax(rpe)[2,])
    } else {
      rpe <- terra::app(area.branch$LR,
                        function(x){
                          if(all(is.na(x))){
                            return(NA)
                          }
                          sum(x, na.rm = T)
                        })
      rpe <- terra::app(rpe, function(x, m){ # to reescale the values from 0 to 1
        (x/m)
      }, m = terra::minmax(rpe)[2,])
    }

    names(rpe) <- c("PE")
  }

  if (!is.null(filename)){ # to save the rasters when the path is provide
    rpe <- terra::writeRaster(rpe, filename, ...)
  }

  return(rpe)
}

#' Calculate Weighted Endemism (WE. Williams et al. 1994, Crisp et al. 2001) for a raster
#'
#' Calculate weighted endemism following (Williams et al. 1994, Crisp et al. 2001) using rasters as input and output
#'
#' @param x SpatRaster. A presence-absence SpatRaster with 0 representing absence and 1 representing presence.
#' @param filename character. Output filename.
#' @param cores numeric. A positive integer indicating the number of clusters used for parallelization
#' @param ... additional arguments to be passed passed down from a calling function.
#' @return SpatRaster
#' @export
#' @references Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J., VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and endemism in conservation. In: Systematics and Conservation Evaluation (eds Forey PL, Humphries CJ, Vane-Wright RI), p. 438. Oxford University Press, Oxford.
#' @references Crisp, M., Laffan, S., Linder, H., Monro, A. (2001). Endemism in theAustralian flora. Journal of Biogeography, 28, 183–198.
#' @examples
#' \dontrun{
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' rast.we(ras)
#' }

rast.we <- function(x, filename = NULL, cores = 1, ...){

  temp <- vector("list", length = 2) # to create a temporary vector with the raster number
  temp[[1]] <- paste0(tempfile(), ".tif")  # to store the first raster

  rs <- range.size(x)
  # inverse of range size
  inv.R <- terra::ifel(x == 0, 0, 1/(x * rs),
                       filename = temp[[1]], overwrite = TRUE) # calculating the inverse of range size
  # weighted endemism
  if(cores>1) {
    # calculating we
    rend <- terra::app(inv.R,
                       function(x){
                         if(all(is.na(x))){
                           return(NA)}
                         sum(x, na.rm = T)
                       }, cores = cores)

    rend <- terra::app(rend, function(x, m){ # to reescale the values
      (x/m)
    }, m = terra::minmax(rend)[2,]) # 2 is the max
  } else {
    rend <- terra::app(inv.R,
                       function(x){
                         if(all(is.na(x))){
                           return(NA)}
                         sum(x, na.rm = T)
                       })

    rend <- terra::app(rend, function(x, m){ # to reescale the values
      (x/m)
    }, m = terra::minmax(rend)[2,])
  }

  names(rend) <- c("WE")

  if (!is.null(filename)){ # to save the rasters when the path is provide
    rend <- terra::writeRaster(rend, filename, ...)
  }
  unlink(temp[[1]]) # delete the archive

  return(rend)
}

#' Calculate richness, phylogenetic diversity, evolutionary distinctiveness, phylogenetic endemism and weighted endemism for a raster
#'
#' Calculate richness, phylogenetic diversity, evolutionary distinctiveness, phylogenetic endemism and weighted endemism using rasters as input and output.
#'
#' @inheritParams rast.ed
#' @return SpatRaster
#' @export
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C., & Cook, L. G. (2009). Phylogenetic endemism: a new approach for identifying geographical concentrations of evolutionary history. Molecular ecology, 18(19), 4061-4072.
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological conservation, 61(1), 1-10.
#' @references Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J., VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and endemism in conservation. In: Systematics and Conservation Evaluation (eds Forey PL, Humphries CJ, Vane-Wright RI), p. 438. Oxford University Press, Oxford.
#' @references Crisp, M., Laffan, S., Linder, H., Monro, A. (2001). Endemism in theAustralian flora. Journal of Biogeography, 28, 183–198.
#' @references Isaac, N. J., Turvey, S. T., Collen, B., Waterman, C. & Baillie, J. E. (2007). Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE 2, e296.
#' @examples
#' \dontrun{
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
#' geo.phylo(ras, tree, metric = 'phylo.diversity')
#' }

geo.phylo <- function(x, tree, metric = c('richness', 'phylo.diversity', 'evol.distinct', 'phylo.endemism', 'weigh.endemism'), filename = NULL, cores = 1, ...){

  {
    # reordering the raster according to tree
    # getting branch length
    # getting ancestor number
    pp <- phylo.pres(x, tree)
    branch.length <- pp$branch.length # branch length
    x <- pp$x # raster reordered
    n.descen <- pp$n.descen # number of descendants by each branch
  }

  {
    ir <- inv.range(x, branch.length) # calculating inverse of range size
    area.inv <- ir$inv.R # subletting only the inverse of range size
    area.tips <- ir$LR # inverse of range size multiplied by branch length
  }

  # if(!all.equal(names(x), names(branch.length))){
  #   stop('Species names are not in the same order on 'x' and 'branch.length' arguments! See 'phylogrid::phylo.pres' function.')
  # }

  # calculating PD, SE, WE, PE, and ED
  if (cores>1){

    if(metric == 'phylo.diversity'){
      # phylogenetic diversity
      resu <- phylogrid::rast.pd(x, branch.length, cores = cores)
      if (!is.null(filename)){ # to save the rasters when the path is provide
        resu <- terra::writeRaster(resu, filename, ...)
      }

    } else if (metric == 'evol.distinct'){
      # evolutionary distinctiveness
      resu <- phylogrid::rast.ed(x, branch.length, n.descen, cores = cores)
      if (!is.null(filename)){ # to save the rasters when the path is provide
        resu <- terra::writeRaster(resu, filename, ...)
      }

    } else if (metric == 'richness'){
      # richness
      resu <- terra::app(x, sum, na.rm = TRUE, cores = cores)
      if (!is.null(filename)){ # to save the rasters when the path is provide
        resu <- terra::writeRaster(resu, filename, ...)
      }

    } else if (metric == 'weigh.endemism'){
      #  weighted endemism
      resu <- phylogrid::rast.we(x, cores = cores)
      if (!is.null(filename)){ # to save the rasters when the path is provide
        resu <- terra::writeRaster(resu, filename, ...)
      }

    } else if (metric == 'phylo.endemism'){
      # phylogenetic endemism
      resu <- phylogrid::rast.pe(x, branch.length, cores = cores)
      if (!is.null(filename)){ # to save the rasters when the path is provide
        resu <- terra::writeRaster(resu, filename, ...)
      }

    } else {
      stop("Choose a valid metric! The community metrics currently available are: 'richness', 'phylo.diversity', 'evol.distinct', 'phylo.endemism', 'weigh.endemism'.")
    }
  } else {
    if(metric == 'phylo.diversity'){
      # phylogenetic diversity
      resu <- phylogrid::rast.pd(x, branch.length, cores = cores)
      if (!is.null(filename)){ # to save the rasters when the path is provide
        resu <- terra::writeRaster(resu, filename, ...)
      }

    } else if (metric == 'evol.distinct'){
      # evolutionary distinctiveness
      resu <- phylogrid::rast.ed(x, branch.length, n.descen, cores = cores)
      if (!is.null(filename)){ # to save the rasters when the path is provide
        resu <- terra::writeRaster(resu, filename, ...)
      }

    } else if (metric == 'richness'){
      # richness
      resu <- terra::app(x, sum, na.rm = TRUE, cores = cores)
      if (!is.null(filename)){ # to save the rasters when the path is provide
        resu <- terra::writeRaster(resu, filename, ...)
      }

    } else if (metric == 'weigh.endemism'){
      # weighted endemism
      resu <- phylogrid::rast.we(x, cores = cores)
      if (!is.null(filename)){ # to save the rasters when the path is provide
        resu <- terra::writeRaster(resu, filename, ...)
      }

    } else if (metric == 'phylo.endemism'){
      # phylogenetic endemism
      resu <- phylogrid::rast.pe(x, branch.length, cores = cores)
      if (!is.null(filename)){ # to save the rasters when the path is provide
        resu <- terra::writeRaster(resu, filename, ...)
      }

    } else {
      stop("Choose a valid metric! The community metrics currently available are: 'richness', 'phylo.diversity', 'evol.distinct', 'phylo.endemism', 'weigh.endemism'.")
    }
  }
  return(resu)
}
