test_that("Are the returned values correct?", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))

  # metric SE richness
  r <- terra::values(phylogrid::rast.se(x))[4423:4460,]
  expect_equivalent(r, c(3, 3, 5, 4, 6, 6, 7, 7, 7, 7, 7, 6, 7, 7, 6, 7, 7, 7, 7, 7,
                         8, 8, 8, 8, 8, 8, 8, 9, 8, 9, 10, 11, 10, 9, 9, 10, 10, 12))
})
