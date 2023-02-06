test_that("Are the returned values correct?", {

  shp <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp",
                                 package="phylogrid"))

  countries <- terra::vect(rnaturalearth::ne_countries()) # mapa mundi
  coun.crop <- terra::crop(countries, terra::ext(shp)) # cut by the total extension of the polygons
  coun.rast <- terra::rasterize(coun.crop,
                         terra::rast(terra::ext(shp), resolution = 0.5))
  # rasterizing with a mask of a country for example
  mvalues.mask <- terra::values(phylogrid::shp2rast(shp, y = coun.rast,
                                                    sps.col = "BINOMIAL", ymask = TRUE,
                                                    background = 0))[6819:6821,1:5]
  expect_equivalent(mvalues.mask, matrix(data = c(0, 0, NA, 1, 1, NA, 0, 0, NA, 0, 0, NA, 0, 0, NA),
                                    nrow = 3, ncol = 5))

  mvalues <- terra::values(phylogrid::shp2rast(shp, sps.col = "BINOMIAL",
                                               ymask = FALSE, background = 0,
                                               resolution = 0.1))[111722:111725,3:6]
  expect_equivalent(mvalues, matrix(data = c(rep(1, 4), rep(0, 4), rep(1, 8)),
                                    nrow = 4, ncol = 4))

})
