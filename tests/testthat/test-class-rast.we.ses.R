test_that("check if the object class is correct", {

  # load data
  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))

  # tests
  expect_s4_class(phyloraster::rast.we.ses(ras[[1:5]], aleats = 3, random = "spat"),
                  "SpatRaster")
})
