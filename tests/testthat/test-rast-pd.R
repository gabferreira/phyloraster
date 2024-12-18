test_that("returned object classes are correct", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))
  data <- phylo.pres(x, tree)

  # tests
  expect_s4_class(rast.pd(data$x, edge.path = data$edge.path,
                          branch.length = data$branch.length),
                  "SpatRaster")
})

test_that("Test that error is returned with wrong order of the
          species names", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157,
                                   -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))
  data <- phylo.pres(x, tree)

  # metric PE
  expect_error(rast.pd(x, edge.path = data$edge.path,
                       branch.length = data$branch.length[5:8]))
})

test_that("results of the analyses replicate those of
          other packages", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))
  data <- phylo.pres(x, tree, full_tree_metr = T)

  # phyloraster
  pg <- rast.pd(data$x, data$tree, full_tree_metr = T)

  # epm
  ep <- terra::rast(system.file("extdata", "epm_PD.tif",
                                package="phyloraster"))

  testthat::expect_equal(matrix(terra::values(pg), ncol=1),
                         matrix(terra::values(ep),  ncol=1),
                         tolerance = 1.43e-07)
})

test_that("results of the analyses replicate those of
          other packages", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))
  data <- phylo.pres(x, tree, full_tree_metr = T)

  # phyloraster
  pg <- rast.pd(data$x, edge.path = data$edge.path,
                branch.length = data$branch.length)

  # epm
  ep <- terra::rast(system.file("extdata", "epm_PD.tif",
                                package="phyloraster"))

  testthat::expect_equal(matrix(terra::values(pg), ncol=1),
                         matrix(terra::values(ep),  ncol=1),
                         tolerance = 1.43e-07)

})

test_that("error is returned when the raster does not have a
          longitude/latitude coordinate reference system (CRS)", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157,
                                   -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))

  w <- terra::project(x, "EPSG:2169")

  data <- phylo.pres(w, tree)

  # metric PE
  expect_error(rast.pd(w, branch.length = data$branch.length))

})
