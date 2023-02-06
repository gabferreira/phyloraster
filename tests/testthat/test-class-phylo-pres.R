test_that("check if the returned object class is correct", {

  # load data
  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
  pp <- phylogrid::phylo.pres(ras, tree)

  # tests
  expect_s4_class(pp$x, "SpatRaster")
  expect_s4_class(pp$subtree, "phylo4")
  expect_type(pp$branch.length, "double")
})
