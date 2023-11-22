test_that("Are the returned values correct", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))

  x <- x[[1:5]]

  rs.obs <- range_size(x, terra::cellSize(x))
  rs.expec <- c(5021016056, 196704341368, 4523101328,
                26914371304, 216777518096)

  expect_equivalent(rs.obs, rs.expec)
})
