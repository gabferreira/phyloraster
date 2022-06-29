test_that("names are equal", {

  # load data
  phylo.pres <- phylogrid::phylo.pres
  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
  datapres <- phylo.pres(ras, tree)

  # test
  expect_equal(all.equal(names(datapres$pres.reord), tree$tip.label), TRUE)
})
