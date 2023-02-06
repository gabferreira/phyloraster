test_that("Are the returned values correct?", {

  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))

  sr <- spat.rand(ras, random = "site")

  expect_equal(terra::minmax(sr)[1], 0)
  expect_equal(terra::minmax(sr)[2], 1)

})
