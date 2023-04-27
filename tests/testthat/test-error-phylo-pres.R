test_that("Test that error is returned with wrong input class", {

  # load data and functions
  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
  pp <- phyloraster::phylo.pres(ras, tree)

  # tests
  expect_error(phyloraster::phylo.pres(ras, pp$branch.length))
})
