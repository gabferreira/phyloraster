test_that("names are equal", {

  # load data
  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
  datapres <- phylogrid::phylo.pres(ras, tree)

  # test
  expect_equal(all.equal(names(datapres$x), tree$tip.label), TRUE)
})
