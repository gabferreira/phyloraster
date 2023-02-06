test_that("Are the returned values correct?", {

  # load data
  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
  pp <- phylogrid::phylo.pres(ras, tree)

  # test
  expect_equivalent(pp$n.descendents[c(1:6, 31:33)], c(37, 37, 38, 38,
                                                       37, 39, 59, 59, 44))
  expect_equivalent(terra::values(pp$x)[6114:6130], c(0, 0, 0, 0, 0, 0, 0,
                                                      NaN, NaN, 1, 1, NaN, NaN,
                                                      NaN, NaN, NaN, NaN))
  expect_equivalent(pp$branch.length[c(1:6, 31:33)], c(0.589016, 0.589015, 0.395144,
                                                       0.395144, 0.589015, 0.490836,
                                                       0.267000, 0.267000, 0.998948))
})
