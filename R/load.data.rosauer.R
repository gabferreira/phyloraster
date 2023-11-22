#' Load an example dataset with presence-absence data of 33 tree frogs and
#' a phylogenetic tree for this species
#'
#' @description This function load a phylogenetic tree, a raster and a
#' data.frame with presence-absence of 33 Australian tree frogs from Rosauer
#' (2017). We also provide distribution shapefiles for ten species according
#' to the IUCN.
#' @return data.frame, SpatRaster, SpatVector and phylo
#' @usage load.data.rosauer()
#' @export
#' @source Rosauer, 2017. Available on:
#' \href{https://github.com/DanRosauer/phylospatial/tree/master/PhyloEndemism_in_R/Tree%20Frog%20Data/}{Github}
#' @source IUCN. 2022. The IUCN Red List of Threatened Species (spatial data).
#' Version 2022-1. \href{https://www.iucnredlist.org}{IUCN}
#' @export
load.data.rosauer <- function(){
  x <- phyloraster::dataR
  r <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))
  s <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp",
                               package="phyloraster"))
  return(list(presab = x, raster = r, IUCN_shapefile = s, tree = tree))
}
