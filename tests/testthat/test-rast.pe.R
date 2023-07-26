test_that("returned values are correct", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))

  # getting fewer cells to test all values
  r <- terra::rast(extent = c(150.0157, 150.8157, -23.044, -22.8563))
  xcrop <- terra::crop(x, terra::ext(r))

  # metric PE
  pe.obs <- round(terra::values(rast.pe(xcrop, tree)), 7)
  pe.expect <- matrix(data = c(0.2277756, 0.2277756, 0.2277756, 0.2345368,
                               0.3594504, 0.3594504, 0.3594504,
                              0.3594504 ,0.2276112, 0.2276112, 0.1883587,
                              0.2276112, 0.3524346, 0.3591910,
                             0.3591910, 0.3591910), ncol = 1)
  expect_equivalent(pe.obs, pe.expect)
})

test_that("error is returned with wrong order of the species names", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
  data <- phylo.pres(x, tree)

  # metric PE
  expect_error(rast.pe(data$x[[1:4]], branch.length = data$branch.length[5:8]))
  expect_error(rast.pe(data$x[[1:4]]))
})

test_that("results of the analyses replicate those of other packages", {

  requireNamespace("epm")
  requireNamespace("phyloraster")
  requireNamespace("phyloregion")

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
  data <- phylo.pres(x, tree)

  # phyloraster
  pg <- rast.pe(x, tree)

  # # phyloregion
  # fdir <- system.file("extdata/rasters", package="phyloraster")
  # files <- file.path(fdir, dir(fdir))
  # com <- phyloregion::raster2comm(files)
  # pr <- phyloregion::phylo_endemism(com$comm_dat, tree)
  # m <- merge(com$poly_shp, data.frame(grids=names(pr), PE=pr), by="grids")
  # m <- m[!is.na(m$PE),]
  #
  # # transform to raster to make this comparable
  # r <- terra::rasterize(terra::vect(m), x, field = "PE")

  # epm
  require(epm)
  datepm <- epm::createEPMgrid(x, resolution = 0.01)
  data$tree <- as(data$tree, "phylo")
  datepm <- epm::addPhylo(datepm, data$tree)
  ep <- epm::gridMetrics(datepm, metric = "phyloWeightedEndemism")

  testthat::expect_equivalent(round(terra::values(pg), 10), round(terra::values(ep$grid$phyloWeightedEndemism)), 10)
  # testthat::expect_equivalent(round(terra::values(pg), 10), round(terra::values(m$PE), 10))
})

test_that("returned object classes are correct", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))

  # tests
  expect_s4_class(rast.pe(x, tree), "SpatRaster")
})

test_that("error is returned when the raster does not have a longitude/latitude coordinate reference system (CRS)", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))

  w <- terra::project(x, "EPSG:2169")

  data <- phylo.pres(w, tree)

  # metric PE
  expect_error(rast.pe(w, branch.length = data$branch.length))

})
