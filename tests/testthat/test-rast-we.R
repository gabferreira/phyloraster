test_that("returned object classes are correct", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044, -22.8563)))

  # tests
  expect_s4_class(rast.we(x), "SpatRaster")
})

test_that("results of the analyses replicate those of other packages", {

  library(phyloraster)

  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))

  # phyloraster
  pg <- rast.we(x)

  # epm
  ep <- terra::rast(system.file("extdata", "epm_WE.tif",
                                package="phyloraster"))

  testthat::expect_equal(matrix(terra::values(pg), ncol=1),
                         matrix(terra::values(ep),  ncol=1),
                         tolerance = 0.002)
})

test_that("error is returned when the raster does not have a longitude/latitude
          coordinate reference system (CRS)", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044, -22.8563)))

  w <- terra::project(x, "EPSG:2169")

  # metric PE
  expect_error(rast.we(w))

})

