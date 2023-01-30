test_that("check if the object class is correct", {

  # load data
  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
  data <- phylogrid::phylo.pres(x = ras, tree = tree)
  sa <- phylogrid::tip.rand(data$branch.length, aleats = 3)

  # tests
  expect_type(sa, "list")
})
