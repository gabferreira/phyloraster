test_that("returned object classes are correct", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
  data <- phylogrid::phylo.pres(x, tree)

  # tests
  expect_s4_class(phylogrid::rast.pe(data$x, data$branch.length), "SpatRaster")
})
