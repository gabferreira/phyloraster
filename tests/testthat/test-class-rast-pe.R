test_that("returned object classes are correct", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
  data <- phyloraster::phylo.pres(x, tree)

  # tests
  expect_s4_class(phyloraster::rast.pe(data$x, data$branch.length), "SpatRaster")
})
