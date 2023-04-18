test_that("Are the returned values correct?", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))

  # getting fewer cells to test all values
  r <- terra::rast()
  terra::ext(r) <- c(150.0157, 150.8157, -23.044, -22.8563)
  xcrop <- terra::crop(x, terra::ext(r))

  # metric SE richness
  riq.pres <- phylogrid::rast.se(xcrop)
  riq.fut <- phylogrid::rast.se(xcrop[[c(1:9)]]) # imagine we lost some species in the future
  dg.obs <- terra::values(phylogrid::delta.grid(riq.pres, riq.fut))
  c(dg.obs)
  dg.expect <- matrix(data = c(-8,  -8,  -8,  -9, -10, -10, -10, -10,  -8,  -8,
                               -8,  -8,  -9, -10, -10, -10), ncol = 1, byrow = F)

  expect_equivalent(dg.obs, dg.expect)
})
