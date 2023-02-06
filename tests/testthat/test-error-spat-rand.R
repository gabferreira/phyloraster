test_that("Test that error is returned with wrong argument names", {

  # load data and functions
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))

  # tests
  expect_error(phylogrid::spat.rand(data$IUCN_shapefile, random = "sites"))
  expect_error(phylogrid::spat.rand(data$presab, random = "sites"))
  expect_error(phylogrid::spat.rand(x, random = "sites"))
  expect_error(phylogrid::spat.rand(x, random = "specie"))
  expect_error(phylogrid::spat.rand(x, random = "Both"))
})
