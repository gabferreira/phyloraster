#' Function to evaluate if the rasters generated in the function fits on memory
#'
#' @param x number of rasters generated in the function
#'
#' @return logical
# #' @export
#' @examples
.fit.memory <- function(x){
  # x rasters will be generated in this function, let's see if there is enough memory in the user's pc
  sink(nullfile())    # suppress output
  mi <- terra::mem_info(x, 1)[5] != 0 # proc in memory = T TRUE means that it fits in the pc's memory, so you wouldn't have to use temporary files
  sink()
  return(mi)
}

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
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' rse <- rast.se(x)
#' terra::plot(rse)
#' }
rast.se <- function(x, filename = NULL, cores = 1, ...){

  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

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
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
#' data <- phylo.pres(x, tree)
#' rast.pd(data$x, data$branch.length)
#' }
rast.pd <- function(x, branch.length, filename = NULL, cores = 1, ...){

  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  if(!all.equal(names(x), names(branch.length))){

    stop("Species names are not in the same order on 'x' and 'branch.length' arguments! See 'phyloraster::phylo.pres' function.")

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
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' # phylogenetic tree
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
#' data <- phyloraster::phylo.pres(x, tree)
#' ed <- phyloraster::rast.ed(data$x, data$branch.length, data$n.descen, cores = 2)
#' }
rast.ed <- function(x, branch.length, n.descen, cores = 1, filename = NULL, ...){

  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  if(!all.equal(names(x), names(branch.length))){

    stop("Species names are not in the same order on 'x' and 'branch.length' arguments! See 'phyloraster::phylo.pres' function.")

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
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
#' data <- phylo.pres(x, tree)
#' rast.pe(data$x, data$branch.length, cores = 1)
#' }

rast.pe <- function(x, branch.length, cores = 1, filename = NULL, ...){

  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  if(!all.equal(names(x), names(branch.length))){

    stop("Species names are not in the same order on 'x' and 'branch.length' arguments! See 'phyloraster::phylo.pres' function.")

  } else {

    # 3 rasters will be generated in this function, let's see if there is enough memory in the user's pc
    sink(nullfile())    # suppress output
    mi <- terra::mem_info(x, 3)[5] != 0 # proc in memory = T TRUE means that it fits in the pc's memory, so you wouldn't have to use temporary files
    sink()

    temp <- vector("list", length = 3) # to create a temporary vector with the raster number
    temp[[1]] <- paste0(tempfile(), ".tif")  # to store the first raster
    temp[[2]] <- paste0(tempfile(), ".tif")  # to store the second raster
    temp[[3]] <- paste0(tempfile(), ".tif")  # to store the third raster

    area.branch <- phyloraster::inv.range(x, branch.length, LR = T,
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

  # if(rescale == TRUE){
  #   rpe <- terra::app(rpe, function(x, m){ # rescale the values from 0 to 1
  #     (x/m)
  #   }, m = terra::minmax(rpe)[2,], filename = ifelse(mi, "", temp[[3]]))
  #
  # }

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
#' x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' rast.we(x)
#' }
#'
rast.we <- function(x, cores = 1, filename = NULL, ...){

  if(!terra::is.lonlat(x)){
    stop("Geographic coordinates are needed for the calculations.")
  }

  # 2 rasters will be generated in this function, let's see if there is enough memory in the user's pc
  sink(nullfile())    # suppress output
  mi <- terra::mem_info(x, 2)[5] != 0 # proc in memory = T TRUE means that it fits in the pc's memory, so you wouldn't have to use temporary files
  sink()

  # if(!mi){
    temp <- vector("list", length = 2) # to create a temporary vector with the raster number
    temp[[1]] <- paste0(tempfile(), ".tif")  # to store the first raster
    temp[[2]] <- paste0(tempfile(), ".tif")  # to store the second raster
  # }

  # inverse of range size
  inv.R <- inv.range(x, filename = ifelse(mi, "", temp[[1]])) # calculate the inverse of range size multiplied by branch length of each species

  # weighted endemism
  # calculating we
  rend <- terra::app(inv.R$inv.R,
                     function(x){
                       if(all(is.na(x))){
                         return(NA)}
                       sum(x, na.rm = T)
                     }, cores = cores, filename = ifelse(mi, "", temp[[2]]))

  # if(rescale == TRUE){
  #   rend <- terra::app(rend, function(x, m){ # rescale the values
  #     (x/m)
  #   }, m = terra::minmax(rend)[2,], cores = cores, filename = ifelse(mi, "", temp[[3]])) # 2 is the max
  # }

  names(rend) <- "WE" # layer name

  if (!is.null(filename)){ # to save the rasters when the path is provide
    rend <- terra::writeRaster(rend, filename = filename)
  }

  unlink(temp[[1]]) # delete the files

  return(rend)
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
#' branch.length = data$branch.length
#' n.descen = data$n.descendants
#' area.branch <- phyloraster::inv.range(data$x, data$branch.length, LR = T)
#'
#' t <- geo.phylo(x, area.branch$LR, area.branch$inv.R, data$branch.length, data$n.descendants)
#' terra::plot(t)
#'
#' @export
geo.phylo <- function(x, area.branch, inv.R,
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
  x3 <- c(x, area.branch, inv.R)

  terra::app(x3,
             function(x, branch.length, n.descen, spp_seq, spp_seqLR,
                      spp_seqINV, resu){
               if(all(is.na(x))){
                 return(resu)
               }
               ### separating raster layers of each cell
               xa <- x[spp_seq]
               xLR <- x[spp_seqLR]
               xINV <- x[spp_seqINV]

               ### computing metrics
               se <- sum(xa, na.rm = T)
               pd <- .vec.pd(xa, branch.length)[[1]] # [[1]] return only the first raster
               ed <- .evol.distin(xa, branch.length, n.descen)[[1]] # [[1]] return only the first raster
               pe <- sum(xLR, na.rm = T)
               we <- sum(xINV, na.rm = T)

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

