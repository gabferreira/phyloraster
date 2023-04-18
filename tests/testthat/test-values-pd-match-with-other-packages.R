test_that("results of the analyses replicate those of other packages", {

  library(epm)
  library(phylogrid)
  library(phyloregion)

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
  data <- phylogrid::phylo.pres(x, tree)

  # phylogrid
  pg <- phylogrid::rast.pd(data$x, data$branch.length)

  # phyloregion
  fdir <- system.file("extdata/rasters", package="phylogrid")
  files <- file.path(fdir, dir(fdir))
  com <- phyloregion::raster2comm(files)
  pr <- phyloregion::PD(com$comm_dat, tree)
  m <- merge(com$poly_shp, data.frame(grids=names(pr), PD=pr), by="grids")
  m <- m[!is.na(m$PD),]

  # transform to raster to make this comparable
  r <- terra::rasterize(terra::vect(m), x, field = "PD")

  # epm
  datepm <- epm::createEPMgrid(x, resolution = 0.01)
  data$subtree <- as(data$subtree, "phylo")
  datepm <- epm::addPhylo(datepm, data$subtree)
  ep <- epm::gridMetrics(datepm, metric = "pd")

  testthat::expect_equivalent(round(values(pg), 10), round(values(ep$grid$pd)), 10)
  testthat::expect_equivalent(round(values(pg), 10), round(values(r$PD)), 10)
})
