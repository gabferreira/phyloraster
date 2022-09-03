#' Calculate phylogenetic diversity (Faith 1992) for a vector
#'
#' @param pres.ab numeric. A Named numeric vector of presence-absence
#' @param branch.length numeric. A Named numeric vector of branch length for each specie
#'
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological conservation, 61(1), 1-10.
#' @return numeric
# #' @export
#' @examples
.vec.pd <- function(pres.ab, branch.length){
  pres.ab[is.na(pres.ab)] <- 0

  if(sum(pres.ab)== 0) {
    return(c(NA,NA))
  }

  if(sum(pres.ab) != 0) {
    pres <- pres.ab == 1 # only species present
    species <- names(branch.length)
    n.species <- length(species)
    pd <- sum(branch.length[pres]) # pd Faith 1992
    pdr <- pd/sum(pres) # pd relative to richness
  }
  return(c(pd = pd, pdr = pdr))
}

#' Calculate phylogenetic diversity, phylogenetic endemism and WE for a raster
#'
#' Calculate phylogenetic diversity, phylogenetic endemism and weighted endemism using rasters as input and output.
#'
#' @param x SpatRaster. A presence-absence SpatRaster with the layers ordered according to the tree order.
#' @param area.inv SpatRaster. A presence-absence SpatRaster with the inverse of the area for each specie.
#' @param area.tips SpatRaster. A presence-absence SpatRaster with the inverse of the area vs branch length.
#' @param branch.length numeric. A Named numeric vector containing the branch length of each specie.
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @return SpatRaster
#' @export
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C., & Cook, L. G. (2009). Phylogenetic endemism: a new approach for identifying geographical concentrations of evolutionary history. Molecular ecology, 18(19), 4061-4072.
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological conservation, 61(1), 1-10.
#' @references Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J., VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and endemism in conservation. In: Systematics and Conservation Evaluation (eds Forey PL, Humphries CJ, Vane-Wright RI), p. 438. Oxford University Press, Oxford.
#' @references Crisp, M., Laffan, S., Linder, H., Monro, A. (2001). Endemism in theAustralian flora. Journal of Biogeography, 28, 183–198.
#' @examples
#' \dontrun{
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
#' data <- phylo.pres(ras, tree)
#' range <- inv.range(data$x, data$branch.length)
#' geo.phylo(data$x, range$inv.R, range$LR, data$branch.length)
#' }

geo.phylo <- function(x, area.inv, area.tips, branch.length, filename = NULL, ...){

  if(!all.equal(names(x), names(branch.length))){
    stop("Species names are not in the same order on 'x' and 'branch.length' arguments! See 'phylogrid::phylo.pres' function.")
  }

  if(all.equal(names(x), names(branch.length))){ # to check wether names match

    # phylogenetic diversity and richness-relative phylogenetic diversity
    rpd <- terra::app(x, fun = .vec.pd,
                      branch.length = branch.length)

    # weighted endemism
    {
      rend <- terra::app(area.inv,
                         function(x){
                           if(all(is.na(x))){
                             return(NA)}
                           sum(x, na.rm = T)
                         })


      rend <- terra::app(rend, function(x, m){ # to reescale the values
        (x/m)
      }, m = terra::minmax(rend)[2,])
    }


    # phylogenetic endemism
    {
      rpe <- terra::app(area.tips,
                        function(x){
                          if(all(is.na(x))){
                            return(NA)
                          }
                          sum(x, na.rm = T)
                        })
      rpe <- terra::app(rpe, function(x, m){ # to reescale the values
        (x/m)
      }, m = terra::minmax(rpe)[2,])
    }

    # raster stack
    gp <- c(rpd, rend, rpe)
    names(gp) <- c("PD", "PDR", "WE", "PE")
  }

  if (!is.null(filename)){ # to save the rasters when the path is provide
    gp <- terra::writeRaster(gp, filename, ...)
  }

  return(gp)
}

#' Calculate phylogenetic diversity, phylogenetic endemism and weighted endemism for a raster- version 2
#'
#' Calculate phylogenetic diversity, phylogenetic endemism and weighted endemism using rasters as input and output.
#'
#' @param rast.branch SpatRaster. a presence-absence SpatRaster with the layers ordered according to the tree order and a vector with the branch lengths. These objects can be obtained in the phylo.pres function
#' @param range.branch SpatRaster. a presence-absence SpatRaster with the inverse of the area and the inverse of the area vs branch length. These rasters can be obtained in the inv.range function
#' @param filename character. Output filename.
#' @return SpatRaster
# #' @export
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C., & Cook, L. G. (2009). Phylogenetic endemism: a new approach for identifying geographical concentrations of evolutionary history. Molecular ecology, 18(19), 4061-4072.
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological conservation, 61(1), 1-10.
#' @references Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J., VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and endemism in conservation. In: Systematics and Conservation Evaluation (eds Forey PL, Humphries CJ, Vane-Wright RI), p. 438. Oxford University Press, Oxford.
#' @references Crisp, M., Laffan, S., Linder, H., Monro, A. (2001). Endemism in theAustralian flora. Journal of Biogeography, 28, 183–198.
#' @examples
#' \dontrun{
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
#' data <- phylo.pres(ras, tree)
#' range <- inv.range(data$x, data$branch.length)
#' geo.phylo2(data, range)
#' }
geo.phylo2 <- function(rast.branch, range.branch, filename = NULL){
  gp2 <- geo.phylo(rast.branch$x,
                   range.branch$inv.R, range.branch$LR,
                  rast.branch$branch.length)
  return(gp2)
}


#' Calculate phylogenetic diversity (PD. Faith, 1992) for a raster
#'
#' Calculate phylogenetic diversity using rasters as input and output. This function follows Faith (1992).
#'
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
#' range <- inv.range(data$x, data$branch.length)
#' rast.pd(data$x, data$branch.length)
#' }
rast.pd <- function(x, branch.length, filename = NULL, ...){

  if(!all.equal(names(x), names(branch.length))){
    stop("Species names are not in the same order on 'x' and 'branch.length' arguments! See 'phylogrid::phylo.pres' function.")
  }

  if(all.equal(names(x), names(branch.length))){ # to check wether names match

    # phylogenetic diversity and richness-relative phylogenetic diversity
    rpd <- terra::app(x, fun = .vec.pd,
                      branch.length = branch.length)
    names(rpd) <- c("PD", "PDR")
  }

  if (!is.null(filename)){ # to save the rasters when the output filename is provide
    rpd <- terra::writeRaster(rpd, filename, ...)
  }

  return(rpd)
}

#' Calculate phylogenetic endemism (PE) for a raster
#'
#' Calculate phylogenetic endemism using rasters as input and output.
#'
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
rast.pe <- function(x, branch.length, filename = NULL, ...){

  if(!all.equal(names(x), names(branch.length))){
    stop("Species names are not in the same order on 'x' and 'branch.length' arguments! See 'phylogrid::phylo.pres' function.")
  }

  if(all.equal(names(x), names(branch.length))){ # to check wether names match

    area.branch <- phylogrid::inv.range(x, branch.length)

    # phylogenetic endemism
    message("Calculating the phylogenetic endemism")
    {
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

#' Calculate Weighted Endemism (WE) for a raster
#'
#' Calculate weighted endemism using rasters as input and output.
#'
#' @param x SpatRaster. A presence-absence SpatRaster with 0 representing absence and 1 representing presence.
#' @param range.size Named numeric vector. A vector containing range size. See the function phylogrid::range.size.
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @return SpatRaster
#' @export
#' @references Williams, P.H., Humphries, C.J., Forey, P.L., Humphries, C.J., VaneWright, R.I. (1994). Biodiversity, taxonomic relatedness, and endemism in conservation. In: Systematics and Conservation Evaluation (eds Forey PL, Humphries CJ, Vane-Wright RI), p. 438. Oxford University Press, Oxford.
#' @references Crisp, M., Laffan, S., Linder, H., Monro, A. (2001). Endemism in theAustralian flora. Journal of Biogeography, 28, 183–198.
#' @examples
#' \dontrun{
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' rs <- phylogrid::range.size(ras)
#' rast.we(ras, rs)
#' }
rast.we <- function(x, range.size, filename = NULL, ...){

  temp <- vector("list", length = 1) # to create a temporary vector with the raster number
  temp[[1]] <- paste0(tempfile(), ".tif")  # to store the first raster

  # inverse of range size
  inv.R <- terra::ifel(x == 0, 0, 1/(x * range.size),
                       filename = temp[[1]], overwrite = TRUE) # calculating the inverse of range size
  # weighted endemism
  { # calculating we
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
