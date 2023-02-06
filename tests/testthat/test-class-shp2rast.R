test_that("check if the object class is correct", {

  # load data
  shp <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp", package="phylogrid"))
  sr <- shp2rast(shp, sps.col = "BINOMIAL", ymask = FALSE, background = 0, resolution = 0.5)

  # tests
  expect_s4_class(sr, "SpatRaster")
})
