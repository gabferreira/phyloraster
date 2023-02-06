test_that("Test that names in the raster and in the branch lenght are in the same order", {

  # load data
  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
  datapres <- phylogrid::phylo.pres(ras, tree)

  # test
  expect_equal(all.equal(names(datapres$x), tree$tip.label), TRUE)
})
