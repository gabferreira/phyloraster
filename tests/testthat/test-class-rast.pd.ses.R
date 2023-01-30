test_that("check if the object class is correct", {

  # load data
  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
  data <- phylogrid::phylo.pres(ras, tree)
  t <- phylogrid::rast.pd.ses(data$x, data$branch.length, aleats = 3, random = "site")

  # tests
  expect_s4_class(t, "SpatRaster")
})
