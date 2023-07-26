test_that("are the returned values correct?", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))

  # getting fewer cells to test all values
  r <- terra::rast()
  terra::ext(r) <- c(150.0157, 150.8157, -23.044, -22.8563)
  xcrop <- terra::crop(x, terra::ext(r))

  # metric SE richness
  se.obs <- terra::values(rast.sr(xcrop))
  expect_equivalent(se.obs, c(12, 12, 12, 13, 14, 14, 14, 14, 12, 12, 11, 12,
                         13, 14, 14, 14))

})

test_that("returned object classes are correct", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))

  # tests
  expect_s4_class(rast.sr(x), "SpatRaster")

})

test_that("error is returned when the raster does not have a longitude/latitude coordinate reference system (CRS)", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))

  w <- terra::project(x, "EPSG:2169")

  # metric PE
  expect_error(rast.sr(w))

})
