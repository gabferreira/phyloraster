test_that("Test that error is returned when the argument n.descendents is missing", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
  data <- phylogrid::phylo.pres(x, tree)

  # metric PE
  expect_error(phylogrid::rast.ed(x, data$branch.length, data$n.decendents))
  expect_error(phylogrid::rast.ed(data$x, data$n.descendents))
  expect_error(phylogrid::rast.ed(data$x, tree))

})
