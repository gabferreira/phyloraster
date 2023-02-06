test_that("Test that error is returned with wrong input class", {

  # load data and functions
  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
  pp <- phylogrid::phylo.pres(ras, tree)

  # tests
  expect_error(phylogrid::phylo.pres(ras, pp$branch.length))
})
