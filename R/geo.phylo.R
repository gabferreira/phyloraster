#' Calculate phylogenetic diversity (Faith 1992) for a vector
#'
#' @description This function calculates the sum of the branch length for a set of species for one sample.
#' @param x numeric. A named numerical vector of presence-absence for one sample.
#' @param branch.length numeric. A named numerical vector containing the branch length for each species.
#' @author Neander Marcel Heming and Gabriela Alves-Ferreira
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological conservation, 61(1), 1-10.
#' @return numeric
# #' @export

.vec.pd <- function(x, branch.length){

  x[is.na(x)] <- 0 # 0 for all value = NA

  if(sum(x)== 0) { # return NA if x = 0
    return(c(NA,NA))
  }

  if(sum(x) != 0) { # if the sum of x is non-zero then do this:
    pres <- x == 1 # only species present in the vector
    # species <- names(branch.length)
    # n.species <- length(species)
    pd <- sum(branch.length[pres]) # pd Faith 1992
    pd1 <- pd # terra::app function does not work when it returns only one raster
  }

  return(c(pd = pd, pd1 = pd1))

}

#' Calculate Evolutionary Distinctiveness for a vector
#'
#' @description This function calculates evolutionary distinctiveness for a set of species using the fair-proportion index (Isaac et al., 2007).
#' @usage .evol.distin(x, branch.length, n.descen)
#' @param x numeric. A Named numeric vector of presence-absence
#' @param branch.length numeric. A Named numeric vector of branch length for each specie
#' @param n.descen numeric. A Named numeric vector of number of descendants for each branch
#' @author Gabriela Alves-Ferreira and Neander Marcel Heming
#' @references Isaac, N. J., Turvey, S. T., Collen, B., Waterman, C. and Baillie, J. E. (2007). Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE 2, e296.
#' @return numeric
# #' @export

.evol.distin <- function(x, branch.length, n.descen){

  x[is.na(x)] <- 0 # 0 for all value = NA

  if(sum(x) == 0) { # return NA if x = 0
    return(c(NA,NA))
  }

  if(sum(x) != 0){ # if the sum of x is non-zero then do this:
    pres <- x == 1 # only species present in the vector
    # species <- names(x)
    ed <- sum(branch.length[pres]/n.descen[pres]) # evolutionary distinctiveness
    ed1 <- ed # terra::app function does not work when this intern function returns only one raster
  }

  return(c(ed = ed, ed1 = ed1))

}

#' Calculate species richness for each raster cell
#'
#' @description Calculate the species richness for each raster cell.
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
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' rse <- rast.se(ras)
#' terra::plot(rse)
#' }
rast.se <- function(x, filename = NULL, cores = 1, ...){

  # 1 rasters will be generated in this function, let's see if there is enough memory in the user's pc
  sink(nullfile())    # suppress output
  mi <- terra::mem_info(x, 1)[5] != 0 # proc in memory = T TRUE means that it fits in the pc's memory, so you wouldn't have to use temporary files
  sink()

  temp <- vector("list", length = 1) # to create a temporary vector with the raster number
  temp[[1]] <- paste0(tempfile(), ".tif")  # to store the first raster

  # richness
  rse <- terra::app(x, sum, na.rm = TRUE, cores = cores,
                    filename = ifelse(mi, "", temp[[1]]))
  names(rse) <- c("SE")

  if (!is.null(filename)){ # to save the rasters when the output filename is provide
    rse <- terra::writeRaster(rse, filename)
  }

  return(rse)

}

#' Calculate phylogenetic diversity for each raster cell
#'
#' @description Calculate the sum of the branch length for species present in each cell of the raster.
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for a set of species. The layers (species) must be sorted according to the tree order. See the phylo.pres function.
#' @param branch.length numeric. A Named numerical vector containing the branch length for a set of species.
#' @param filename character. Output filename.
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster with that many cores is created and used.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @author Neander Marcel Heming and Gabriela Alves-Ferreira
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological conservation, 61(1), 1-10.
#' @return SpatRaster
#' @export
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

    # 1 rasters will be generated in this function, let's see if there is enough memory in the user's pc
    sink(nullfile())    # suppress output
    mi <- terra::mem_info(x, 1)[5] != 0 # proc in memory = T TRUE means that it fits in the pc's memory, so you wouldn't have to use temporary files
    sink()

    temp <- vector("list", length = 1) # to create a temporary vector with the raster number
    temp[[1]] <- paste0(tempfile(), ".tif")  # to store the first raster

    # phylogenetic diversity
    rpd <- terra::app(x, fun = .vec.pd,
                      branch.length = branch.length, cores = cores,
                      filename = ifelse(mi, "", temp[[1]]))
    rpd <- rpd[[1]] # select only the first raster
    names(rpd) <- c("PD")
  }

  if (!is.null(filename)){ # to save the rasters when the output filename is provide
    rpd <- terra::writeRaster(rpd, filename)
  }

  return(rpd)
}


#' Calculate Evolutionary distinctiveness for each raster cell
#'
#' @description This function calculates evolutionary distinctiveness according to the fair-proportion index.
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for a set of species. The layers (species) must be sorted according to the tree order. See the phylo.pres function.
#' @param branch.length numeric. A Named numeric vector containing the branch length of each specie.
#' @param n.descen numeric. A Named numeric vector of number of descendants for each branch
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster with that many cores is created and used.
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed for fun.
#' @author Gabriela Alves-Ferreira and Neander Marcel Heming
#' @references Isaac, N. J., Turvey, S. T., Collen, B., Waterman, C. and Baillie, J. E. (2007). Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE 2, e296.
#' @return SpatRaster
#' @export
#' @examples
#' \dontrun{
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' # phylogenetic tree
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
#' data <- phylogrid::phylo.pres(x, tree)
#' ed <- phylogrid::rast.ed(data$x, data$branch.length, data$n.descen, cores = 2)
#' }
#'

rast.ed <- function(x, branch.length, n.descen, cores = 1, filename = NULL, ...){

  if(!all.equal(names(x), names(branch.length))){

    stop("Species names are not in the same order on 'x' and 'branch.length' arguments! See 'phylogrid::phylo.pres' function.")

  } else {

    # 1 rasters will be generated in this function, let's see if there is enough memory in the user's pc
    sink(nullfile())    # suppress output
    mi <- terra::mem_info(x, 1)[5] != 0 # proc in memory = T TRUE means that it fits in the pc's memory, so you wouldn't have to use temporary files
    sink()

    temp <- vector("list", length = 1) # to create a temporary vector with the raster number
    temp[[1]] <- paste0(tempfile(), ".tif")  # to store the first raster

    # evolutionary distinctiveness
    red <- terra::app(x, fun = .evol.distin,
                      branch.length, n.descen, cores = cores,
                      filename = ifelse(mi, "", temp[[1]]))
    red <- red[[1]] # only the first raster
    names(red) <- "ED" # layer name
  }

  if(!is.null(filename)){ # to save the rasters when the output filename is provide
    red <- terra::writeRaster(red, filename)
  }

  return(red)

}

#' Calculate phylogenetic endemism for each raster cell
#'
#' @description Calculate the sum of the inverse of the range size multiplied by the branch length for the species present in each raster cell.
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for a set of species. The layers (species) must be sorted according to the tree order. See the phylo.pres function.
#' @param branch.length numeric. A Named numeric vector containing the branch length of each specie
#' @param rescale logical. If TRUE, the values are scaled from 0 to 1. The default is FALSE.
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster with that many cores is created and used.
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @author Gabriela Alves-Ferreira and Neander Marcel Heming
#' @references Laffan, S. W., Rosauer, D. F., Di Virgilio, G., Miller, J. T., González‐Orozco, C. E., Knerr, N., ... & Mishler, B. D. (2016). Range‐weighted metrics of species and phylogenetic turnover can better resolve biogeographic transition zones. Methods in Ecology and Evolution, 7(5), 580-588.
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C. and Cook, L. G. (2009). Phylogenetic endemism: a new approach for identifying geographical concentrations of evolutionary history. Molecular ecology, 18(19), 4061-4072.
#' @return SpatRaster
#' @export
#' @examples
#' \dontrun{
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
#' data <- phylo.pres(ras, tree)
#' rast.pe(data$x, data$branch.length, rescale = FALSE, cores = 1)
#' }

rast.pe <- function(x, branch.length, rescale = FALSE, cores = 1, filename = NULL, ...){

  if(!all.equal(names(x), names(branch.length))){

    stop("Species names are not in the same order on 'x' and 'branch.length' arguments! See 'phylogrid::phylo.pres' function.")

  } else {

    # 3 rasters will be generated in this function, let's see if there is enough memory in the user's pc
    sink(nullfile())    # suppress output
    mi <- terra::mem_info(x, 3)[5] != 0 # proc in memory = T TRUE means that it fits in the pc's memory, so you wouldn't have to use temporary files
    sink()

    temp <- vector("list", length = 3) # to create a temporary vector with the raster number
    temp[[1]] <- paste0(tempfile(), ".tif")  # to store the first raster
    temp[[2]] <- paste0(tempfile(), ".tif")  # to store the second raster
    temp[[3]] <- paste0(tempfile(), ".tif")  # to store the third raster

    area.branch <- phylogrid::inv.range(x, branch.length,
                                        filename = ifelse(mi, "", temp[[1]])) # calculate the inverse of range size multiplied by branch length of each species

    # phylogenetic endemism
    rpe <- terra::app(area.branch$LR,
                      function(x){
                        if(all(is.na(x))){
                          return(NA)
                        }
                        sum(x, na.rm = T)
                      }, cores = cores, filename = ifelse(mi, "", temp[[2]]))

  }

  if(rescale == TRUE){
    rpe <- terra::app(rpe, function(x, m){ # rescale the values from 0 to 1
      (x/m)
    }, m = terra::minmax(rpe)[2,], filename = ifelse(mi, "", temp[[3]]))

  }

  names(rpe) <- "PE" # layer name


  if (!is.null(filename)){ # to save the rasters when the path is provide
    rpe <- terra::writeRaster(rpe, filename)
  }

  return(rpe)

}

#' Calculate weighted endemism for each raster cell
#'
#' @description Calculate the sum of the inverse of the range size for species present in each raster cell.
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for a set of species.
#' @param rescale logical. If TRUE, the values are scaled from 0 to 1. The default is FALSE.
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster with that many cores is created and used.
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @author Neander Marcel Heming and Gabriela Alves Ferreira
#' @return SpatRaster
#' @references Laffan, S. W., Rosauer, D. F., Di Virgilio, G., Miller, J. T., González‐Orozco, C. E., Knerr, N., ... & Mishler, B. D. (2016). Range‐weighted metrics of species and phylogenetic turnover can better resolve biogeographic transition zones. Methods in Ecology and Evolution, 7(5), 580-588.
#' @references Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J., VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and endemism in conservation. In: Systematics and Conservation Evaluation (eds Forey PL, Humphries CJ, Vane-Wright RI), p. 438. Oxford University Press, Oxford.
#' @references Crisp, M., Laffan, S., Linder, H., Monro, A. (2001). Endemism in theAustralian flora. Journal of Biogeography, 28, 183–198.
#' @export
#' @examples
#' \dontrun{
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' rast.we(ras)
#' }
#'
rast.we <- function(x, rescale = FALSE, cores = 1, filename = NULL, ...){

  # 1 rasters will be generated in this function, let's see if there is enough memory in the user's pc
  sink(nullfile())    # suppress output
  mi <- terra::mem_info(x, 3)[5] != 0 # proc in memory = T TRUE means that it fits in the pc's memory, so you wouldn't have to use temporary files
  sink()

  temp <- vector("list", length = 3) # to create a temporary vector with the raster number
  temp[[1]] <- paste0(tempfile(), ".tif")  # to store the first raster
  temp[[2]] <- paste0(tempfile(), ".tif")  # to store the second raster
  temp[[3]] <- paste0(tempfile(), ".tif")  # to store the third raster

  rs <- range_size(x) # return an vector with species's range
  cell.area <- terra::cellSize(terra::rast(x), filename = temp[[1]]) # to calculate cell size

  # inverse of range size
  inv.R <- terra::ifel(x == 0, 0, cell.area/(x * rs),
                       filename = ifelse(mi, "", temp[[1]]), overwrite = TRUE) # calculating the inverse of range size

  # weighted endemism
  # calculating we
  rend <- terra::app(inv.R,
                     function(x){
                       if(all(is.na(x))){
                         return(NA)}
                       sum(x, na.rm = T)
                     }, cores = cores, filename = ifelse(mi, "", temp[[2]]))

  if(rescale == TRUE){
    rend <- terra::app(rend, function(x, m){ # rescale the values
      (x/m)
    }, m = terra::minmax(rend)[2,], cores = cores, filename = ifelse(mi, "", temp[[3]])) # 2 is the max
  }

  names(rend) <- "WE" # layer name

  if (!is.null(filename)){ # to save the rasters when the path is provide
    rend <- terra::writeRaster(rend, filename = filename)
  }

  unlink(temp[[1]]) # delete the files

  return(rend)
}

#' Calculate community metrics for each raster cell
#'
#' Calculate richness, phylogenetic diversity, evolutionary distinctiveness, phylogenetic endemism and weighted endemism using rasters as input and output.
#' @param x SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for a set of species. The layers (species) must be sorted according to the tree order. See the phylo.pres function.
#' @param tree phylo. A dated tree.
#' @param metric character. Name of metric to use, available metrics are: 'richness', 'phylo.diversity', 'evol.distinct', 'phylo.endemism' and 'weigh.endemism'. See Details for more information.
#' @param rescale logical. If TRUE, the values are scaled from 0 to 1. The default is FALSE.
#' @param cores positive integer. If cores > 1, a 'parallel' package cluster with that many cores is created and used.
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed for fun.
#' @return SpatRaster
#' @export
#' @details Community metrics available:
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
#' @examples
#' \dontrun{
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
#' geo.phylo(ras, tree, metric = 'phylo.diversity')
#' geo.phylo(ras, tree, metric = 'phylo.endemism')
#' geo.phylo(ras, metric = 'weigh.endemism')
#' }
#'
geo.phylo <- function(x, tree, metric = c('richness', 'phylo.diversity',
                                          'evol.distinct', 'phylo.endemism',
                                          'weigh.endemism'), rescale = FALSE, cores = 1,
                      filename = NULL, ...){

  ### calculating PD, SE, WE, PE, and ED

  if(metric == 'phylo.diversity'){

    # 1 rasters will be generated in this function, let's see if there is enough memory in the user's pc
    sink(nullfile())    # suppress output
    mi <- terra::mem_info(x, 1)[5] != 0 # proc in memory = T TRUE means that it fits in the pc's memory, so you wouldn't have to use temporary files
    sink()

    ### preparing data
    # reordering the raster according to tree
    # getting branch length
    pp <- phylo.pres(x, tree)

    # phylogenetic diversity
    resu <- phylogrid::rast.pd(pp$x, pp$branch.length, cores = cores)
    if (!is.null(filename)){ # to save the rasters when the path is provide
      resu <- terra::writeRaster(resu, filename)
    }

  } else if (metric == 'evol.distinct'){

    # 1 rasters will be generated in this function, let's see if there is enough memory in the user's pc
    sink(nullfile())    # suppress output
    mi <- terra::mem_info(x, 1)[5] != 0 # proc in memory = T TRUE means that it fits in the pc's memory, so you wouldn't have to use temporary files
    sink()

    ### preparing data
    # reordering the raster according to tree
    # getting branch length
    # getting ancestor number
    pp <- phylo.pres(x, tree)

    # evolutionary distinctiveness
    resu <- phylogrid::rast.ed(pp$x, pp$branch.length, pp$n.descen, cores = cores)
    if (!is.null(filename)){ # to save the rasters when the path is provide
      resu <- terra::writeRaster(resu, filename)
    }

  } else if (metric == 'richness'){

    # 1 rasters will be generated in this function, let's see if there is enough memory in the user's pc
    sink(nullfile())    # suppress output
    mi <- terra::mem_info(x, 1)[5] != 0 # proc in memory = T TRUE means that it fits in the pc's memory, so you wouldn't have to use temporary files
    sink()

    # richness
    resu <- terra::app(x, sum, na.rm = TRUE, cores = cores)
    if (!is.null(filename)){ # to save the rasters when the path is provide
      resu <- terra::writeRaster(resu, filename)
    }

  } else if (metric == 'weigh.endemism'){


    # 2 rasters will be generated in this function, let's see if there is enough memory in the user's pc
    sink(nullfile())    # suppress output
    mi <- terra::mem_info(x, 2)[5] != 0 # proc in memory = T TRUE means that it fits in the pc's memory, so you wouldn't have to use temporary files
    sink()

    temp <- vector("list", length = 2) # to create a temporary vector with the raster number
    temp[[1]] <- paste0(tempfile(), ".tif")  # to store the first raster
    temp[[2]] <- paste0(tempfile(), ".tif")  # to store the second raster

    rs <- range_size(x) # return an vector with species's range
    cell.area <- terra::cellSize(terra::rast(x), filename = temp[[1]]) # to calculate cell size

    # inverse of range size
    inv.R <- terra::ifel(x == 0, 0, cell.area/(x * rs),
                         filename = ifelse(mi, "", temp[[1]]), overwrite = TRUE) # calculating the inverse of range size

    # weighted endemism
    # calculating we
    resu <- terra::app(inv.R,
                       function(x){
                         if(all(is.na(x))){
                           return(NA)}
                         sum(x, na.rm = T)
                       }, cores = cores, filename = ifelse(mi, "", temp[[2]]))

    names(resu) <- "WE" # layer name

    if (!is.null(filename)){ # to save the rasters when the path is provide
      rend <- terra::writeRaster(rend, filename = filename)
    }

    unlink(temp[[1]]) # delete the files

  }
    else if (metric == 'phylo.endemism'){

    # 1 rasters will be generated in this function, let's see if there is enough memory in the user's pc
    sink(nullfile())    # suppress output
    mi <- terra::mem_info(x, 1)[5] != 0 # proc in memory = T TRUE means that it fits in the pc's memory, so you wouldn't have to use temporary files
    sink()

    ### preparing data
    # reordering the raster according to tree
    # getting branch length
    pp <- phylogrid::phylo.pres(x, tree)

    # phylogenetic endemism
    resu <- phylogrid::rast.pe(pp$x, pp$branch.length, cores = cores, rescale = rescale)
    if (!is.null(filename)){ # to save the rasters when the path is provide
      resu <- terra::writeRaster(resu, filename)
    }

  } else {
    stop("Choose a valid metric! The community metrics currently available are: 'richness', 'phylo.diversity', 'evol.distinct', 'phylo.endemism', 'weigh.endemism'.")
  }
  return(resu)
}

