test_that("check if the returned object class is correct", {

  # load data
  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
  data <- phylogrid::phylo.pres(ras, tree)
  ir <- phylogrid::inv.range(data$x, data$branch.length)

  # tests
  expect_type(ir, "list")
})
