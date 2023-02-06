test_that("Are the returned values correct", {

  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  phylogrid::range.size(ras, scale = TRUE)

  expect_equivalent(phylogrid::range.size(ras, scale = TRUE)[1:5], c(0.01500945, 0.58801338,
                                                                0.01352102, 0.08045583,
                                                                0.64801865))
})
