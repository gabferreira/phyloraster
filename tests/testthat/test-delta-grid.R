test_that("check if the object class is correct", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))

  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157,
                                   -23.044, -22.8563)))

  # metric SE richness
  riq.pres <- phyloraster::rast.sr(x)
  riq.fut <- phyloraster::rast.sr(x[[c(1:15)]]) # imagine we
  #lost some species in the future
  dg <- phyloraster::delta.grid(riq.pres, riq.fut)

  # tests
  expect_s4_class(dg, "SpatRaster")

})

test_that("Are the returned values correct?", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))

  # getting fewer cells to test all values
  xcrop <- terra::crop(x, terra::ext(c(150.0157, 150.8157,
                                       -23.044, -22.8563)))

  # metric SE richness
  riq.pres <- phyloraster::rast.sr(xcrop)
  riq.fut <- phyloraster::rast.sr(xcrop[[c(1:9)]]) # imagine
  #we lost some species in the future
  dg.obs <- terra::values(phyloraster::delta.grid(riq.pres, riq.fut))
  c(dg.obs)
  dg.expect <- matrix(data = c(-8,  -8,  -8,  -9, -10, -10,
                               -10, -10,  -8,  -8,
                               -8,  -8,  -9, -10, -10, -10),
                      ncol = 1, byrow = F)

  expect_equivalent(dg.obs, dg.expect)
})

test_that("error is returned when the extent do not match", {

  data <- phyloraster::load.data.rosauer()
  shp <- data$IUCN_shapefile

  rt <- shp2rast(shp, sps.col = "BINOMIAL", ymask = FALSE,
                              background = NA, resolution = 0.1)
  richrt <- phyloraster::rast.sr(rt)

  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  richx <- phyloraster::rast.sr(x)

  expect_error(delta.grid(richrt, richx))

})

test_that("Raster is saved when filename is provided", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))

  # getting fewer cells to test all values
  xcrop <- terra::crop(x, terra::ext(c(150.0157, 150.8157,
                                       -23.044, -22.8563)))

  # metric SE richness
  riq.pres <- phyloraster::rast.sr(xcrop)
  riq.fut <- phyloraster::rast.sr(xcrop[[c(1:9)]]) # imagine
  #we lost some species in the future

  temp <- tempfile(fileext = ".tif")

  # rasterizing with a mask of a country for example
  expect(phyloraster::delta.grid(riq.pres, riq.fut,
                                 filename = temp), ok = T)
  unlink(temp)
})


