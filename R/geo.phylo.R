#' Calculate phylogenetic diversity (Faith 1992) for a vector
#'
#' @param pres.ab numeric. A Named num vector of presence-absence
#' @param branch.length numeric. A Named num vector of branch length for each specie
#'
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

#' Calculate PD, PE and WE for a raster
#'
#' Calculate phylogenetic diversity, phylogenetic endemism and weighted endemism using rasters as input and output
#'
#' @param pres.reord SpatRaster. A presence-absence SpatRaster with the layers ordered according to the tree order
#' @param area.inv SpatRaster. A presence-absence SpatRaster with the inverse of the area for each specie
#' @param area.tips SpatRaster. A presence-absence SpatRaster with the inverse of the area vs branch length
#' @param branch.length numeric. A Named numeric vector containing the branch length of each specie
#' @param filename character. A character specifying the path where rasters can be saved
#' @param ... arguments to be passed to the function
#' @return SpatRaster
#' @export
#' @examples
#' \dontrun{
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
#' data <- phylo.pres(ras, tree)
#' range <- inv.range(data$pres.reord, data$branch.length)
#' geo.phylo(data$pres.reord, range$inv.R, range$LR, data$branch.length)
#' }

geo.phylo <- function(pres.reord, area.inv, area.tips, branch.length, filename = NULL, ...){

  if(!all.equal(names(pres.reord), names(branch.length))){
    stop("Species names are not in the same order on 'pres.reord' and 'branch.length' arguments! See 'phylogrid::phylo.pres' function.")
  }

  if(all.equal(names(pres.reord), names(branch.length))){ # to check wether names match

    # phylogenetic diversity and richness-relative phylogenetic diversity
    rpd <- terra::app(pres.reord, fun = .vec.pd,
                      branch.length = branch.length)

    # weighted endemism
    {
      rend <- terra::app(area.inv,
                         function(x){
                           if(all(is.na(x))){
                             return(NA)}
                           sum(x, na.rm = T)
                         })


      rend <- terra::app(rend, function(x, m){
        (x/m)
      }, m = terra::minmax(rend)[2,]) # to reescale values
    }


    ## phylogenetic endemism
    {
      rpe <- terra::app(area.tips,
                        function(x){
                          if(all(is.na(x))){
                            return(NA)
                          }
                          sum(x, na.rm = T)
                        })
      rpe <- terra::app(rpe, function(x, m){
        (x/m)
      }, m = terra::minmax(rpe)[2,])
    }

    # Juntando todos os rasters
    gp <- c(rpd, rend, rpe)
    names(gp) <- c("PD", "PDR", "WE", "PE")
  }

  if (!is.null(filename)){ # to save the rasters when the path is provide
    gp <- terra::writeRaster(gp, filename, ...)
  }

  return(gp)
}


#' Calculate PD, PE and WE for a raster- version 2
#'
#' Calculate phylogenetic diversity, phylogenetic endemism and weighted endemism using rasters as input and output
#'
#' @param rast.branch SpatRaster. a presence-absence SpatRaster with the layers ordered according to the tree order and a vector with the branch lengths. These objects can be obtained in the phylo.pres function
#' @param range.branch SpatRaster. a presence-absence SpatRaster with the inverse of the area and the inverse of the area vs branch length. These rasters can be obtained in the inv.range function
#' @param filename character. Output filename
#' @return SpatRaster
#' @export
#' @examples
#' \dontrun{
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
#' data <- phylo.pres(ras, tree)
#' range <- inv.range(data$pres.reord, data$branch.length)
#' geo.phylo2(data, range)
#' }
geo.phylo2 <- function(rast.branch, range.branch, filename = NULL){
  gp2 <- geo.phylo(rast.branch$pres.reord,
                   range.branch$inv.R, range.branch$LR,
                  rast.branch$branch.length)
  return(gp2)
}
