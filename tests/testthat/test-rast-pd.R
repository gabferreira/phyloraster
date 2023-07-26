test_that("returned object classes are correct", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
  data <- phylo.pres(x, tree)

  # tests
  expect_s4_class(rast.pd(data$x, branch.length = data$branch.length), "SpatRaster")
})

test_that("Test that error is returned with wrong order of the species names", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
  data <- phylo.pres(x, tree)

  # metric PE
  expect_error(rast.pd(x, branch.length = data$branch.length[5:8]))
})

test_that("results of the analyses replicate those of other packages", {

  library(epm)
  library(phyloraster)
  library(phyloregion)

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
  data <- phylo.pres(x, tree)

  # phyloraster
  pg <- rast.pd(data$x, branch.length = data$branch.length)

  # # phyloregion
  # fdir <- system.file("extdata/rasters", package="phyloraster")
  # files <- file.path(fdir, dir(fdir))
  # com <- phyloregion::raster2comm(files)
  # pr <- phyloregion::PD(com$comm_dat, tree)
  # m <- merge(com$poly_shp, data.frame(grids=names(pr), PD=pr), by="grids")
  # m <- m[!is.na(m$PD),]
  #
  # # transform to raster to make this comparable
  # r <- terra::rasterize(terra::vect(m), x, field = "PD")

  # epm
  datepm <- epm::createEPMgrid(x, resolution = 0.01)
  data$tree <- as(data$tree, "phylo")
  datepm <- epm::addPhylo(datepm, data$tree)
  ep <- epm::gridMetrics(datepm, metric = "pd")

  testthat::expect_equivalent(round(terra::values(pg), 10), round(terra::values(ep$grid$pd)), 10)
})

test_that("error is returned when the raster does not have a longitude/latitude coordinate reference system (CRS)", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))

  w <- terra::project(x, "EPSG:2169")

  data <- phylo.pres(w, tree)

  # metric PE
  expect_error(rast.pd(w, branch.length = data$branch.length))

})

# test_that("error is returned with wrong order of the species names", {
#
#   x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#   tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
#   data <- phylo.pres(x, tree)
#
#   # metric PE
#   expect_error(rast.pd(x, branch.length = data$branch.length))
# })

