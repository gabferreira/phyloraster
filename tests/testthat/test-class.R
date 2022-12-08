test_that("returned object classes are correct", {

  # load data and functions
  geo.phylo <- phylogrid::geo.phylo

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
  gp <- geo.phylo(x, tree, metric = 'phylo.diversity')

  # tests
  expect_s4_class(gp, "SpatRaster")
})

