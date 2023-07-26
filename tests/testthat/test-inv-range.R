test_that("check if the returned object class is correct", {

  # load data
  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))

  # tests
  expect_type(inv.range(ras), "list")
  expect_s4_class(inv.range(ras)[[2]], "SpatRaster")
})

test_that("error is returned when names does not match in 'x' and 'branch.length'", {

  # load data
  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
  data <- phylo.pres(ras, tree)

  # tests
  expect_error(inv.range(ras, data$branch.length))
})
