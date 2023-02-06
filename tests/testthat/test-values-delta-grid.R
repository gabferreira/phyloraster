test_that("Are the returned values correct?", {
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))

  # metric SE richness
  riq.pres <- phylogrid::rast.se(x)
  riq.fut <- phylogrid::rast.se(x[[c(1:15)]]) # imagine we lost some species in the future
  dgvalues <- terra::values(phylogrid::delta.grid(riq.pres, riq.fut))[3197:3214,]
  expect_equivalent(dgvalues, c(-3, -1, -1, -1, -1, -1,
                               -1, -1, -2, -3, -2, -3,
                               -3, -4, -3, -4, -4, -5))
})
