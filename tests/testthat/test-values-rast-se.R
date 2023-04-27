test_that("Are the returned values correct?", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))

  # getting fewer cells to test all values
  r <- terra::rast()
  terra::ext(r) <- c(150.0157, 150.8157, -23.044, -22.8563)
  xcrop <- terra::crop(x, terra::ext(r))

  # metric SE richness
  se.obs <- terra::values(phyloraster::rast.se(xcrop))
  expect_equivalent(se.obs, c(12, 12, 12, 13, 14, 14, 14, 14, 12, 12, 11, 12,
                         13, 14, 14, 14))
})
