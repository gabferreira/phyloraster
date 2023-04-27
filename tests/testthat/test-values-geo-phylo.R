test_that("Are the returned values correct?", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))

  # getting fewer cells to test all values
  r <- terra::rast()
  terra::ext(r) <- c(150.0157, 150.8157, -23.044, -22.8563)
  xcrop <- terra::crop(x, terra::ext(r))

  # metric SE
  se.obs <- terra::values(phyloraster::geo.phylo(xcrop, metric = "richness"))
  se.expect <- matrix(data = c(12, 12, 12, 13, 14, 14, 14, 14, 12, 12,
                               11, 12, 13, 14, 14, 14),
                      ncol = 1, byrow = F)
  expect_equivalent(se.obs, se.expect)

  # metric PE
  pe.obs <- round(terra::values(phyloraster::geo.phylo(xcrop, tree, metric = "phylo.endemism")), 7)
  pe.expect <- matrix(data = c(0.2277756, 0.2277756, 0.2277756, 0.2345368,
                               0.3594504, 0.3594504, 0.3594504,
                               0.3594504 ,0.2276112, 0.2276112, 0.1883587,
                               0.2276112, 0.3524346, 0.3591910,
                               0.3591910, 0.3591910), ncol = 1)
  expect_equivalent(pe.obs, pe.expect)


  # metric PD
  pd.obs <- terra::values(phyloraster::geo.phylo(xcrop, tree, metric = "phylo.diversity"))
  pd.expect <- matrix(data = c(3.603842, 3.603842, 3.603842, 3.657917, 4.656865,
                               4.656865, 4.656865, 4.656865, 3.603842,
                               3.603842, 3.014827, 3.603842, 4.602790,
                               4.656865, 4.656865, 4.656865),
                      ncol = 1, byrow = F)
  expect_equivalent(pd.obs, pd.expect)

  # metric ED
  ed.obs <- terra::values(phyloraster::geo.phylo(xcrop, tree, metric = "evol.distinct"))
  ed.expect <- matrix(data = c(0.08343477, 0.08343477, 0.08343477, 0.08445505,
                               0.10715842, 0.10715842, 0.10715842,
                               0.10715842, 0.08343477, 0.08343477, 0.06751545,
                               0.08343477, 0.10613813, 0.10715842,
                               0.10715842, 0.10715842),
                      ncol = 1, byrow = F)
  expect_equivalent(ed.obs, ed.expect)
})

