test_that("check if the object class is correct", {
  # load data
  dat <- phylogrid::load.data.rosauer()

  # tests
  expect_s4_class(phylogrid::df2rast(dat$presab), "SpatRaster")
})
