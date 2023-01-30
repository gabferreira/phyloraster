test_that("check if the object class is correct", {

  # load data
  shp <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp", package="phylogrid"))
  spps <- shp$BINOMIAL
  sr <- phylogrid::shp2rast(shp, spps, resolution = 0.1)

  # tests
  expect_s4_class(sr, "SpatRaster")
})
