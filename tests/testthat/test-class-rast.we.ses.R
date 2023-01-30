test_that("check if the object class is correct", {

  # load data
  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  t <- phylogrid::rast.we.ses(ras, aleats = 3, random = "site")

  # tests
  expect_s4_class(t, "SpatRaster")
})
