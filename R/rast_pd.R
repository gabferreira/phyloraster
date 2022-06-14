#' Calculate phylogenetic diversity (Faith 1992) for a vector
#'
#' @param pres.ab a Named num vector of presence-absence
#' @param branch.length a Named num vector of branch length for each specie
#'
#' @return numeric
# #' @export
#'
#' @examples
.vec.pd <- function(pres.ab, branch.length){
  pres.ab[is.na(pres.ab)] <- 0
  if(sum(pres.ab)== 0) {
    return(c(NA,NA))
  }
  pres <- pres.ab == 1 # only species present
  species <- names(branch.length)
  n.species <- length(species)
  pd <- sum(branch.length[pres]) # pd Faith 1992
  pdr <- pd/sum(pres) # pd relative to richness
  return(c(pd = pd, pdr = pdr))
}

#' Data preparation
#'
#' Reorder stack according to tree order, get branch length, inverse area, and inverse area x branch lengths.
#'
#' @param pres.stack a raster of presence-absence. It can be an object of class 'raster' or 'SpatRaster'
#' @param tree an object of class 'phylo'
#'
#' @return SpatRaster and numeric vector
#' @export
#'
#' @examples
phylo.pres <- function(pres.stack, tree){
  spat.names <- as.character(names(pres.stack))
  tree.phy4 <- phylobase::phylo4(tree) # to extract tipLabels
  labels <- as.character(phylobase::tipLabels(tree.phy4))
  on.tree <- intersect(spat.names,labels)   # raster species that are on the tree
  subtree <- ape::keep.tip(tree, on.tree)

  stack.reord <- pres.stack[[subtree[["tip.label"]]]] # reorder the species stack according to the species order in the tree
  stack.reord.t <- terra::rast(stack.reord) # SpatRaster
  #names(stack.reord.t) == subtree[["tip.label"]] # to check wether names match
  subtree <- phylobase::phylo4(subtree)
  bl <- as.numeric(phylobase::edgeLength(subtree, 1:phylobase::nTips(subtree)))

  srt <- terra::app(stack.reord.t, function(x){
    ifelse(x == 0, NA, x) # to transform 0 in NA
  })
  srt <- terra::cellSize(srt) # to obtain cell size
  names(srt) <- names(stack.reord.t)

  # the function below added all the presences to know how many pixels the species is present
  ## That is, this function determines the size of the species distribution in number of pixels
  area <- sapply(1:terra::nlyr(srt),
                 function(x, i){
                   terra::expanse(x[[i]])
                 }, x = srt)
  area.inv <- terra::app(srt, function(x, a){
    x/a # to calculate the inverse of area size
  }, a = area)

  area.bl <- terra::app(area.inv, function(x, branch.length){
    x * branch.length
  }, branch.length = bl)

  pp <- list(stack.reord.t, area.inv, area.bl, bl)
  names(pp) <- c("pres.reord", "area.inv", "area.bl", "branch.length")
  return(pp)
}

#' Calculate PD, PE and WE for a raster

#' Calculate phylogenetic diversity, phylogenetic endemism and weighted endemism using rasters as input and output
#' @param pres.reord a presence-absence SpatRaster with the layers ordered according to the tree order
#' @param area.inv a presence-absence SpatRaster with the inverse of the area for each specie
#' @param area.tips a presence-absence SpatRaster with the inverse of the area vs branch length
#' @param branch.length a numeric vector containing the branch length of each specie
#' @param filename a character specifying the path where rasters can be saved
#'
#' @return SpatRaster
#' @export
#'
#' @examples
geo.phylo <- function(pres.reord, area.inv, area.tips, branch.length, filename = NULL){
  rpd <- terra::app(pres.reord, fun = .vec.pd, branch.length = branch.length) # phylogenetic diversity and richness-relative phylogenetic diversity
  rend <- sum(area.inv, na.rm = T) # weighted endemism
  rpe <- sum(area.tips, na.rm = T) # phylogenetic endemism
  gp <- c(rpd, rend, rpe)
  names(gp) <- c("pd", "pdr", "te", "pe")

  if (!is.null(filename)){
    gp <- terra::writeRaster(gp, filename)
  }
  return(gp)
}
