test_that("returned object classes are correct", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157,
                                   -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))

  require(SESraster)
  ses <- rast.pe.ses(x, tree,
                     aleats = 2, metric = "all")

  # tests
  expect_s4_class(canape.rast(ses$p.upper.PE, ses$p.upper.PE.alt,
                              ses$p.upper.RPE, ses$p.lower.RPE), "SpatRaster")

})

test_that("error is returned when an argument is missing", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157,
                                   -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))

  require(SESraster)
  ses <- rast.pe.ses(x, tree, aleats = 10, metric = "all")

  # metric PE
  expect_error(canape.rast(ses$p.upper.PE,
                           # ses$p.upper.PE.alt,
                           ses$p.upper.RPE, ses$p.lower.RPE))

})
