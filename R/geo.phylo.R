#' Calculate phylogenetic diversity (Faith 1992) for a vector
#'
#' @param pres.ab a Named num vector of presence-absence
#' @param branch.length a Named num vector of branch length for each specie
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
#' @inheritParams .vec.pd
#' @param pres.reord a presence-absence SpatRaster with the layers ordered according to the tree order
#' @param area.inv a presence-absence SpatRaster with the inverse of the area for each specie
#' @param area.tips a presence-absence SpatRaster with the inverse of the area vs branch length
#' @param branch.length a Named numeric vector containing the branch length of each specie
#' @param filename a character specifying the path where rasters can be saved
#' @param ... arguments to be passed to the function
#' @return SpatRaster
#' @export
#' @examples
geo.phylo <- function(pres.reord, area.inv, area.tips, branch.length, filename = NULL, ...){
  ### mudar para all equal
  if((names(pres.reord) != names(branch.length))){
    stop("Species names are not in the same order on 'pres.reord' and 'branch.length' arguments! See 'phylogrid::phylo.pres' function.")
  }

  if((names(pres.reord) == names(branch.length))){ # to check wether names match
    rpd <- terra::app(pres.reord, fun = .vec.pd, branch.length = branch.length) # phylogenetic diversity and richness-relative phylogenetic diversity
    rend <- terra::app(area.inv,
                       function(x){
                         if(all(is.na(x))){
                           return(NA)}
                         sum(x, na.rm = T)
                       }) # weighted endemism
    rend <- terra::app(rend, function(x, m){
      (x/m)
    }, m = minmax(rend)[2,])
    rpe <- terra::app(area.tips,
               function(x){
                 if(all(is.na(x))){
                   return(NA)
                   }
                 sum(x, na.rm = T)
               }) # phylogenetic endemism
    rpe <- terra::app(rpe, function(x, m){
      (x/m)
    }, m = minmax(rpe)[2,])
    gp <- c(rpd, rend, rpe)
    names(gp) <- c("pd", "pdr", "we", "pe")
  }

  if (!is.null(filename)){
    gp <- terra::writeRaster(gp, filename, ...)
  }

  return(gp)
}
