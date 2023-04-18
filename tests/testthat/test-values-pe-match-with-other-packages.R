test_that("results of the analyses replicate those of other packages", {

  library(epm)
  library(phylogrid)
  library(phyloregion)

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
  data <- phylogrid::phylo.pres(x, tree)

  # phylogrid
  pg <- phylogrid::rast.pe(data$x, data$branch.length)

  # phyloregion
  fdir <- system.file("extdata/rasters", package="phylogrid")
  files <- file.path(fdir, dir(fdir))
  com <- phyloregion::raster2comm(files)
  pr <- phyloregion::phylo_endemism(com$comm_dat, tree)
  m <- merge(com$poly_shp, data.frame(grids=names(pr), PE=pr), by="grids")
  m <- m[!is.na(m$PE),]

  # transform to raster to make this comparable
  r <- terra::rasterize(terra::vect(m), x, field = "PE")

  # epm
  datepm <- epm::createEPMgrid(x, resolution = 0.01)
  data$subtree <- as(data$subtree, "phylo")
  datepm <- epm::addPhylo(datepm, data$subtree)
  ep <- epm::gridMetrics(datepm, metric = "phyloWeightedEndemism")

  testthat::expect_equivalent(round(values(pg), 10), round(values(ep$grid$phyloWeightedEndemism)), 10)
  testthat::expect_equivalent(round(values(pg), 10), round(values(r$PE)), 10)
})
