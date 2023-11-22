test_that("check if the object class is correct", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044, -22.8563)))

  require("SESraster")
  # tests
  expect_s4_class(rast.we.ses(x, aleats = 3), "SpatRaster")

})

test_that("check if function corrects arguments with wrong names", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))
  data <- phylo.pres(x, tree)
  inv.R <- inv.range(data$x)
  # branch.length <- data$branch.length
  # names(branch.length) <- sample(names(branch.length))
  names(inv.R) <- sample(names(inv.R))

  require("SESraster")
  # tests
  expect_s4_class(rast.we.ses(x, inv.R=inv.R, aleats = 2), "SpatRaster")
  expect_s4_class(rast.we(x, inv.R=inv.R), "SpatRaster")
})

test_that("error is returned when the raster does not have a longitude/latitude
          coordinate reference system (CRS)", {

            x <- terra::rast(system.file("extdata", "rast.presab.tif",
                                         package="phyloraster"))
            # getting fewer cells to test all values
            x <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044,
                                             -22.8563)))

            # tree <- ape::read.tree(system.file("extdata", "tree.nex",
            # package="phyloraster"))
            # data <- phylo.pres(x, tree)
            # area.branch <- inv.range(data$x, data$branch.length)
            w <- terra::project(x, "EPSG:2169")

            require("SESraster")
            # tests
            expect_error(rast.we.ses(x = w,
                                     # tree = data$tree,
                                     # FUN_args = list(range.BL=
                                     # area.branch$range.BL,
                                     # inv.R=area.branch$inv.R,
                                     # branch.length=data$branch.length,
                                     # n.descen = data$n.descendants),
                                     spat_alg = "bootspat_str",
                                     spat_alg_args = list(rprob = NULL,
                                                          rich = NULL,
                                                          fr_prob = NULL),
                                     aleats = 5))
          })


test_that("error is returned when the user choose a randomization method not
          available", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044, -22.8563)))

  # tree <- ape::read.tree(system.file("extdata", "tree.nex",
  # package="phyloraster"))
  # data <- phylo.pres(x, tree)
  # area.branch <- inv.range(data$x, data$branch.length)

  require("SESraster")
  # tests
  expect_error(rast.we.ses(x = data$x,
                           # tree = data$tree,
                           # FUN_args = list(range.BL=area.branch$range.BL,
                           # inv.R=area.branch$inv.R,
                           # branch.length=data$branch.length,
                           # n.descen = data$n.descendants),
                           spat_alg = "bootspat_str",
                           spat_alg_args = list(rprob = NULL,
                                                rich = NULL,
                                                fr_prob = NULL),
                           random = "spatial",
                           aleats = 5))
})
