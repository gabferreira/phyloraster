test_that("returned object classes are correct", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))

  # tests
  expect_s4_class(rast.we(x), "SpatRaster")
})

test_that("results of the analyses replicate those of other packages", {

  library(epm)
  library(phyloraster)
  # library(phyloregion)

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))

  # phyloraster
  pg <- rast.we(x)

  # # phyloregion
  # fdir <- system.file("extdata/rasters", package="phyloraster")
  # files <- file.path(fdir, dir(fdir))
  # com <- phyloregion::raster2comm(files)
  # pr <- phyloregion::weighted_endemism(com$comm_dat)
  # m <- merge(com$poly_shp, data.frame(grids=names(pr), WE=pr), by="grids")
  # m <- m[!is.na(m$WE),]
  # x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  # r <- terra::rasterize(terra::vect(m), x, field = "WE")

  # epm
  library(epm)
  datepm <- epm::createEPMgrid(x, resolution = 0.01)
  ep <- epm::gridMetrics(datepm, metric = "weightedEndemism")

  testthat::expect_equal(matrix(terra::values(pg), ncol=1),
                         matrix(terra::values(ep$grid$weightedEndemism),  ncol=1),
                         tolerance = 0.002)
  # testthat::expect_equivalent(round(values(pg), 10), round(values(r$WE)), 10)
})

test_that("error is returned when the raster does not have a longitude/latitude coordinate reference system (CRS)", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))

  w <- terra::project(x, "EPSG:2169")

  # metric PE
  expect_error(rast.we(w))

})
