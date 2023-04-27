test_that("Are the returned values correct?", {

  shp <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp",
                                 package="phyloraster"))

  # getting fewer cells to test all values
  r <- terra::rast()
  terra::ext(r) <- c(113.380470276, 114.55355835, -28.06026001, -27.65233326)
  shpc <- terra::crop(shp, terra::ext(r))
  terra::plot(shpc)
  unique(shpc$BINOMIAL)

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
