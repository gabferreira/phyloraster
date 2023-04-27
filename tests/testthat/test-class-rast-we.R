test_that("returned object classes are correct", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))

  # tests
  expect_s4_class(phyloraster::rast.we(x), "SpatRaster")
})

