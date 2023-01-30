test_that("check if the returned object class is correct", {

  # load data
  dat <- phylogrid::load.data.rosauer()
  r <- phylogrid::df2rast(dat$presab)

  # tests
  expect_s4_class(r, "SpatRaster")
})
