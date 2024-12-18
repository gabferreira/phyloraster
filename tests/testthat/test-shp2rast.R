test_that("error is returned with arguments missing and wrong imput class", {

  # load data and functions
  shp <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp",
                                 package="phyloraster"))
  shp$BINOMIAL <- NULL
  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))

  # tests
  expect_error(phyloraster::shp2rast(x = shp, sps.col = "BINOMIAL",
                                     ymask = FALSE,
                                   resolution = 0.1, background = 0))
  expect_error(phyloraster::shp2rast(x = x, sps.col = "BINOMIAL",
                                     ymask = FALSE,
                                   resolution = 0.1, background = 0))
})

test_that("check if the object class is correct", {

  # load data
  shp <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp",
                                 package="phyloraster"))
  sr <- shp2rast(shp, sps.col = "BINOMIAL", ymask = FALSE, background = 0,
                 resolution = 0.5)

  # tests
  expect_s4_class(sr, "SpatRaster")
})

test_that("Are the returned values correct?", {

  shp <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp",
                                 package="phyloraster"))

  # getting fewer cells to test all values
  r <- terra::rast(vals=1)
  terra::ext(r) <- c(112.380470276, 115.55355835, -30.06026001, -28.55233326)
  # plot(shp)
  # plot(r, add=T, col="gray")
  shpc <- terra::crop(shp, terra::ext(r))
  shpr <- phyloraster::shp2rast(shpc, sps.col = "BINOMIAL",
                                ymask = FALSE, background = 0,
                                resolution = 0.25)
  # plot(shpr, add=T)
  ob.values <- as.vector(terra::values(shpr))

  # c(ob.values[,3])
  expec.values <- matrix(data = c(1, 1, 1, 1, 1, 1, 1,
                                  1, 1, 1, 1, 1, 1, 1, 1),
                         ncol = 3, byrow = F)
  expec.values <- c(1, 1, 1, 0, 1, 1,
                    1, 0, 0, 1, 1, 1,
                    0, 1, 1, 1, 0, 1,
                    1, 1, 0, 1, 1, 1)
  expect_equivalent(ob.values, expec.values)

})

test_that("function runs ok when a mask is applied", {

  # load data
  shp <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp",
                                 package="phyloraster"))

  # create a polygon to use as mask
  ## with an extent
  e <- terra::ext(113, 123, -43.64, -33.90)
  p <- terra::as.polygons(e, crs="")

  coun.crop <- terra::crop(p, terra::ext(shp)) # cut by the total extension of
  # the polygons
  coun.rast <- terra::rasterize(coun.crop,
                                terra::rast(terra::ext(shp), resolution = 0.5))

  # rasterizing with a mask of a country for example
  expect(phyloraster::shp2rast(shp, y = coun.rast, sps.col = "BINOMIAL",
                               ymask = TRUE,
                               background = 0), ok = T)
})

test_that("Raster is saved when filename is provided", {

  # load data
  shp <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp",
                                 package="phyloraster"))

  # create a polygon to use as mask
  ## with an extent
  e <- terra::ext(113, 123, -43.64, -33.90)
  p <- terra::as.polygons(e, crs="")

  coun.crop <- terra::crop(p, terra::ext(shp)) # cut by the total extension of
  # the polygons
  coun.rast <- terra::rasterize(coun.crop,
                                terra::rast(terra::ext(shp), resolution = 0.5))
  temp <- tempfile(fileext = ".tif")

  # rasterizing with a mask of a country for example
  expect(shp2rast(shp, y = coun.rast, sps.col = "BINOMIAL",
                               background = 0, filename = temp), ok = T)
  unlink(temp)
})

