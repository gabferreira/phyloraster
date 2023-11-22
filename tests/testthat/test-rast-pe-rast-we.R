test_that("PE function returns the same result as WE when the
          tree is a rake", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044,
                                   -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))
  tree <- ape::compute.brlen(ape::stree(terra::nlyr(x), "s"), 1)
  tree$tip.label <- names(x)

  # calculating WE and PE with all branch lengths = 1
  we <- rast.we(x)
  pe <- rast.pe(x, tree)

  # testing if the values are the same for WE and PE
  expect_equivalent(terra::values(we, na.rm=TRUE),
                    terra::values(pe, na.rm=TRUE))

})
