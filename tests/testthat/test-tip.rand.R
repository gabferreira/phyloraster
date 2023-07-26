test_that("check if the object class is correct", {

  # load data
  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
  data <- phyloraster::phylo.pres(x = ras, tree = tree)
  sa <- phyloraster::tip.rand(data$branch.length, aleats = 3)

  # tests
  expect_type(sa, "list")
})
