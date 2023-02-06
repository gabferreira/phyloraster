test_that("check if the object class is correct", {

  # load data
  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  t <- phylogrid::rast.we.ses(ras, aleats = 3, random = "species")

  # tests
  expect_s4_class(phylogrid::rast.we.ses(ras[[1:5]], aleats = 3, random = "species"),
                  "SpatRaster")
  expect_s4_class(phylogrid::rast.we.ses(ras[[1:5]], aleats = 3, random = "site"),
                  "SpatRaster")
  expect_s4_class(phylogrid::rast.we.ses(ras[[1:5]], aleats = 3, random = "both"),
                  "SpatRaster")
})
