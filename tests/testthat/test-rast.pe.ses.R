test_that("check if function corrects names in branch", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))
  data <- phylo.pres(x, tree)

  require("SESraster")
  # tests
  expect_s4_class(rast.pe.ses(x, tree, aleats = 2, random = "spat"),
                  "SpatRaster")

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
  branch.length <- data$branch.length
  names(branch.length) <- sample(names(branch.length))
  names(inv.R) <- sample(names(inv.R))

  require("SESraster")
  # tests
  expect_s4_class(rast.pe.ses(x, tree, branch.length=branch.length,
                              inv.R=inv.R, aleats = 2, random = "spat"),
                  "SpatRaster")
  expect_s4_class(rast.pe(x, tree, branch.length=branch.length,
                          inv.R=inv.R), "SpatRaster")

})

test_that("error is returned when the raster does not have a longitude/latitude
          coordinate reference system (CRS)", {

            x <- terra::rast(system.file("extdata", "rast.presab.tif",
                                         package="phyloraster"))
            # getting fewer cells to test all values
            x <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044,
                                             -22.8563)))

            tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                               package="phyloraster"))
            # data <- phylo.pres(x, tree)
            # area.branch <- inv.range(data$x, data$branch.length)
            w <- terra::project(x, "EPSG:2169")

            data <- phyloraster::phylo.pres(w, tree)
            # branch.length <- data$branch.length
            # n.descen <- data$n.descendants
            # area.branch <- phyloraster::inv.range(data$x,
            # data$branch.length, LR = T)

            require("SESraster")
            # tests
            expect_error(rast.pe.ses(x = data$x,
                                     tree = data$tree,
                                     # FUN_args = list(range.BL=
                                     # area.branch$range.BL,
                                     # inv.R=area.branch$inv.R,
                                     # branch.length=data$branch.length,
                                     # n.descen = data$n.descendants),
                                     spat_alg = "bootspat_str",
                                     spat_alg_args = list(rprob = NULL,
                                                          rich = NULL,
                                                          fr_prob = NULL),
                                     aleats = 3))
          })

test_that("error is returned when only argument x is provided", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044, -22.8563)))

  require("SESraster")
  # tests
  expect_error(rast.pe.ses(x = data$x,
                           # tree = data$tree,
                           # FUN_args = list(range.BL=area.branch$range.BL,
                           # inv.R=area.branch$inv.R,
                           # branch.length=data$branch.length,
                           # n.descen = data$n.descendants),
                           spat_alg = "bootspat_str",
                           spat_alg_args = list(rprob = NULL,
                                                rich = NULL,
                                                fr_prob = NULL),
                           aleats = 5))
})

test_that("error is returned when the user choose a randomization
          method not available", {

            x <- terra::rast(system.file("extdata", "rast.presab.tif",
                                         package="phyloraster"))
            # getting fewer cells to test all values
            x <- terra::crop(x, terra::ext(c(150.0157, 150.8157,
                                             -23.044, -22.8563)))

            tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                               package="phyloraster"))
            # data <- phylo.pres(x, tree)
            # area.branch <- inv.range(data$x, data$branch.length)

            require("SESraster")
            # tests
            expect_error(rast.pe.ses(x = x,
                                     tree = tree,
                                     # FUN_args = list(range.BL=
                                     # area.branch$range.BL,
                                     # inv.R=area.branch$inv.R,
                                     # branch.length=data$branch.length,
                                     # n.descen = data$n.descendants),
                                     spat_alg = "bootspat_str",
                                     spat_alg_args = list(rprob = NULL,
                                                          rich = NULL,
                                                          fr_prob = NULL),
                                     random = "spatial",
                                     aleats = 3))
          })

test_that("function runs ok with the method 'tip'", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))
  # data <- phylo.pres(x, tree)
  # area.branch <- inv.range(data$x, data$branch.length)

  expect(rast.pe.ses(x = x,
                     tree = tree,
                     # FUN_args = list(range.BL=area.branch$range.BL,
                     # inv.R=area.branch$inv.R,
                     # branch.length=data$branch.length,
                     # n.descen = data$n.descendants),
                     spat_alg = "bootspat_str",
                     spat_alg_args = list(rprob = NULL,
                                          rich = NULL,
                                          fr_prob = NULL),
                     random = "tip",
                     aleats = 5), ok = T)
})

