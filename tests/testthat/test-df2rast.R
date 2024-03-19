test_that("check if the object class is correct", {

  # load data
  dat <- phyloraster::load.data.rosauer()

  # tests
  expect_s4_class(phyloraster::df2rast(dat$presab),
                  "SpatRaster")
  expect_s4_class(phyloraster::df2rast(dat$presab, CRS = "EPSG:2169"),
                  "SpatRaster")
})
