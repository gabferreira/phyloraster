test_that("Test that error is returned with arguments missing and wrong imput class", {

  # load data and functions
  shp <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp", package="phylogrid"))
  shp$BINOMIAL <- NULL
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))

  # tests
  expect_error(phylogrid::shp2rast(x = shp, sps.col = "BINOMIAL", ymask = FALSE,
                                   resolution = 0.1, background = 0))
  expect_error(phylogrid::shp2rast(x = x, sps.col = "BINOMIAL", ymask = FALSE,
                                   resolution = 0.1, background = 0))
})

