#' Load an example dataset with presence-absence data of 33 tree frogs and a phylogenetic tree for this species
#'
#' This function load a phylogenetic tree, a raster and a data.frame with presence-absence of
#' 33 Australian tree frogs from Rosauer (2017). We also provide distribution shapefiles for each species according to the IUCN..
#'
#' @return data.frame, SpatRaster, SpatVector and phylo
#' @export
#' @source Rosauer, 2017. Available on: <https://github.com/DanRosauer/phylospatial/tree/master/PhyloEndemism_in_R/Tree%20Frog%20Data>
#' @source IUCN. 2022. The IUCN Red List of Threatened Species (spatial data). Version 2022-1. https://www.iucnredlist.org. Accessed on [29 September 2022].
#' @export
load.data.rosauer <- function(){
  x <- phylogrid::dataR
  r <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
  s <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp", package="phylogrid"))
  return(list(presab = x, raster = r, IUCN_shapefile = s, tree = tree))
}
