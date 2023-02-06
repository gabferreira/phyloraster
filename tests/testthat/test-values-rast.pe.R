test_that("Are the returned values correct?", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
  data <- phylogrid::phylo.pres(x, tree)

  # metric PE
  pe <- terra::values(phylogrid::rast.pe(data$x, data$branch.length))[2995:2999,]
  expect_equivalent(pe, c(7.026690e-12, 2.603572e-12, 1.183266e-12,
                          1.273825e-11, 1.142087e-11))
})
