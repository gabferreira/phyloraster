test_that("check if the object class is correct", {
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))

  # metric SE richness
  riq.pres <- phyloraster::rast.se(x)
  riq.fut <- phyloraster::rast.se(x[[c(1:15)]]) # imagine we lost some species in the future
  dg <- phyloraster::delta.grid(riq.pres, riq.fut)

  # tests
  expect_s4_class(dg, "SpatRaster")

})
