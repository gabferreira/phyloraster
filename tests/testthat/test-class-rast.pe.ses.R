test_that("check if the object class is correct", {

  # load data
  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
  data <- phylo.pres(ras, tree)

  # tests
  expect_s4_class(rast.pe.ses(data$x[[1:3]], data$branch.length[1:3], aleats = 2, random = "tip"), "SpatRaster")
  expect_s4_class(rast.pe.ses(data$x[[1:3]], data$branch.length[1:3], aleats = 2, random = "spat"), "SpatRaster")

})

