test_that("check if the object class is correct", {
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))

  # metric SE richness
  riq.pres <- phylogrid::rast.se(x)
  riq.fut <- phylogrid::rast.se(x[[c(1:15)]]) # imagine we lost some species in the future
  dg <- phylogrid::delta.grid(riq.pres, riq.fut)

  # tests
  expect_s4_class(dg, "SpatRaster")

})
