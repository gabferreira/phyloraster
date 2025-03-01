test_that("returned values are correct", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  # getting fewer cells to test all values
  # x <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))

  # getting fewer cells to test all values
  r <- terra::rast(extent = c(150.0157, 150.8157, -23.044, -22.8563))
  xcrop <- terra::crop(x, terra::ext(r))

  # metric PE
  pe.obs <- round(terra::values(rast.pe(x, tree, full_tree_metr = T)), 7)

  # PE from epm package
  epm.pe <- round(terra::values(terra::rast(system.file("extdata", "epm_PE.tif",
                                package="phyloraster"))), 7)

  expect_equivalent(pe.obs, epm.pe, tolerance = 0.00571)
})

test_that("error is returned with wrong order of the species names", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))
  data <- phylo.pres(x, tree, full_tree_metr = TRUE)

  # metric PE
  expect_error(rast.pe(data$x[[1:4]], branch.length = data$branch.length[5:8]))
})

test_that("error is returned when tree or edge.path and branch lenght are
          not supplied", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044, -22.8563)))

  # metric PE
  expect_error(rast.pe(data$x[[1:4]]))
})

test_that("results of the analyses replicate those of other packages", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))
  data <- phylo.pres(x, tree)

  # phyloraster
  pg <- rast.pe(x, tree)

  # # phyloregion
  # fdir <- system.file("extdata/rasters", package="phyloraster")
  # files <- file.path(fdir, dir(fdir))
  # com <- phyloregion::raster2comm(files)
  # pr <- phyloregion::phylo_endemism(com$comm_dat, tree)
  # m <- merge(com$poly_shp, data.frame(grids=names(pr), PE=pr), by="grids")
  # m <- m[!is.na(m$PE),]
  #
  # # transform to raster to make this comparable
  # r <- terra::rasterize(terra::vect(m), x, field = "PE")

  # epm
  # epm
  ep <- terra::rast(system.file("extdata", "epm_PE.tif",
                                package="phyloraster"))


  testthat::expect_equal(matrix(terra::values(pg), ncol=1),
                         matrix(terra::values(ep),  ncol=1),
                         tolerance = 0.01)
})

test_that("returned object classes are correct", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))

  # tests
  expect_s4_class(rast.pe(x, tree), "SpatRaster")


})

test_that("error is returned when the raster does not have a longitude/latitude
          coordinate reference system (CRS)", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))

  w <- terra::project(x, "EPSG:2169")

  data <- phylo.pres(w, tree)

  # metric PE
  expect_error(rast.pe(w, branch.length = data$branch.length))

})
