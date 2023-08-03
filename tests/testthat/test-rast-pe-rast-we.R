test_that("PE function returns the same result as WE when the tree is a rake", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))

  tree$edge.length <- c(rep(1, 58)) # setting all branches to length equal to one

  # calculating WE and PE with all branch lengths = 1
  we <- rast.we(x)
  pe <- rast.pe(x, tree)

  # testing if the values are the same for WE and PE
  expect_equivalent(terra::values(we), terra::values(pe))

})
