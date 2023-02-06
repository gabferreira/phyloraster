test_that("returned object classes are correct", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))

  # tests
  expect_s4_class(phylogrid::rast.se(x), "SpatRaster")
})
