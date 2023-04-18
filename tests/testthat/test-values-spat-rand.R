test_that("Are the returned values correct?", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  # getting fewer cells to test all values
  r <- terra::rast()
  terra::ext(r) <- c(150.0157, 150.8157, -23.044, -22.8563)
  xcrop <- terra::crop(x[[1:3]], terra::ext(r))

  set.seed(7)
  sr.obs <- terra::values(phylogrid::spat.rand(xcrop, random = "site"))
  # c(sr.obs[,3])
  sr.expect <- matrix(data = c(1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                              0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0,
                               0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1), ncol = 3, byrow = F)

  expect_equivalent(sr.obs, sr.expect)
})
