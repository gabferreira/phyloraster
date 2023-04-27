test_that("Test that error is returned with wrong order of the species names", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
  data <- phyloraster::phylo.pres(x, tree)

  # metric PE
  expect_error(phyloraster::rast.pd(x, data$branch.length[5:8]))
})
