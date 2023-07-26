test_that("check if the returned object class is correct", {

  # load data
  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))

  # tests
  # expect_type(inv.range(ras), "list")
  expect_s4_class(inv.range(ras), "SpatRaster")
})
