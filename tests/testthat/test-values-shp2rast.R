test_that("Are the returned values correct?", {

  shp <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp",
                                 package="phylogrid"))

  mvalues <- terra::values(phylogrid::shp2rast(shp, sps.col = "BINOMIAL",
                                               ymask = FALSE, background = 0,
                                               resolution = 0.1))[111722:111725,3:6]
  expect_equivalent(mvalues, matrix(data = c(rep(1, 4), rep(0, 4), rep(1, 8)),
                                    nrow = 4, ncol = 4))

})
