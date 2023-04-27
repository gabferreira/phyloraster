test_that("Are the returned values correct?", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))

  # getting fewer cells to test all values
  r <- terra::rast()
  terra::ext(r) <- c(150.0157, 150.8157, -23.044, -22.8563)
  xcrop <- terra::crop(x, terra::ext(r))
  data <- phyloraster::phylo.pres(xcrop, tree)

  # metric PE
  pe.obs <- round(terra::values(phyloraster::rast.pe(data$x, data$branch.length)), 7)
  pe.expect <- matrix(data = c(0.2277756, 0.2277756, 0.2277756, 0.2345368,
                               0.3594504, 0.3594504, 0.3594504,
                              0.3594504 ,0.2276112, 0.2276112, 0.1883587,
                              0.2276112, 0.3524346, 0.3591910,
                             0.3591910, 0.3591910), ncol = 1)
  expect_equivalent(pe.obs, pe.expect)
})
