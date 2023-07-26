test_that("error is returned with arguments missing and wrong imput class", {

  # load data and functions
  shp <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp", package="phyloraster"))
  shp$BINOMIAL <- NULL
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))

  # tests
  expect_error(phyloraster::shp2rast(x = shp, sps.col = "BINOMIAL", ymask = FALSE,
                                   resolution = 0.1, background = 0))
  expect_error(phyloraster::shp2rast(x = x, sps.col = "BINOMIAL", ymask = FALSE,
                                   resolution = 0.1, background = 0))
})

test_that("check if the object class is correct", {

  # load data
  shp <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp", package="phyloraster"))
  sr <- shp2rast(shp, sps.col = "BINOMIAL", ymask = FALSE, background = 0, resolution = 0.5)

  # tests
  expect_s4_class(sr, "SpatRaster")
})

test_that("Are the returned values correct?", {

  shp <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp",
                                 package="phyloraster"))

  # getting fewer cells to test all values
  r <- terra::rast()
  terra::ext(r) <- c(113.380470276, 114.55355835, -28.06026001, -27.65233326)
  shpc <- terra::crop(shp, terra::ext(r))

  ob.values <- terra::values(phyloraster::shp2rast(shpc, sps.col = "BINOMIAL",
                                                   ymask = FALSE, background = 0,
                                                   resolution = 0.1))

  # c(ob.values[,3])
  expec.values <- matrix(data = c(0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0),
                         ncol = 3, byrow = F)
  expect_equivalent(ob.values, expec.values)

})

test_that("function runs ok when a mask is applied", {

  # load data
  require(rnaturalearth)
  shp <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp",
                                package="phyloraster"))
  countries <- terra::vect(rnaturalearth::ne_countries()) # world map
  coun.crop <- terra::crop(countries, terra::ext(shp)) # cut by the total extension of the polygons
  coun.rast <- terra::rasterize(coun.crop,
                        terra::rast(terra::ext(shp), resolution = 0.5))

  # rasterizing with a mask of a country for example
  expect(shp2rast(shp, y = coun.rast, sps.col = "BINOMIAL", ymask = TRUE,
                  background = 0), ok = T)
})
