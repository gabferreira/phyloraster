test_that("check if the object class is correct", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  sr <- phylogrid::spat.rand(x, random = "site")

  # tests
  expect_s4_class(sr, "SpatRaster")
})
