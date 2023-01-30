test_that("check if the returned object class is correct", {

  # load data and functions
  phylo.pres <- phylogrid::phylo.pres

  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
  pp <- phylo.pres(ras, tree)

  # tests
  expect_type(pp, "list")
})
