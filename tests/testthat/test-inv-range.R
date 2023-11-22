test_that("check if the returned object class is correct", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))

  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157,
                                   -23.044, -22.8563)))

  # tests
  # expect_type(inv.range(ras), "list")
  expect_s4_class(inv.range(x), "SpatRaster")
})
