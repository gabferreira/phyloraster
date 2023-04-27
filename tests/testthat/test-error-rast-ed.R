test_that("Test that error is returned when the argument n.descendents is missing", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
  data <- phyloraster::phylo.pres(x, tree)

  # metric PE
  expect_error(phyloraster::rast.ed(x, data$branch.length, data$n.decendents))
  expect_error(phyloraster::rast.ed(data$x, data$n.descendents))
  expect_error(phyloraster::rast.ed(data$x, tree))

})
