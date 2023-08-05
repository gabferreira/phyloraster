test_that("check if the object class is correct", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
  data <- phyloraster::phylo.pres(x = x, tree = tree)
  sa <- phyloraster::tip.rand(data$branch.length, aleats = 3)

  # tests
  expect_type(sa, "list")
})
