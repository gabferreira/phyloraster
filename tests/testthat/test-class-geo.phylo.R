test_that("returned object classes are correct", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
  gp <- phyloraster::geo.phylo(x, tree, metric = 'phylo.diversity')

  # tests
  expect_s4_class(gp, "SpatRaster")
})

