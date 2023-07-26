test_that("returned object class is correct", {

            x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
            tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
            # data <- phylo.pres(x, tree)
            # area.branch <- inv.range(data$x, data$branch.length)

            expect_s4_class(geo.phylo.ses(x = x,
                                          tree = tree,
                                          # FUN_args = list(range.BL=area.branch$range.BL,
                                          # inv.R=area.branch$inv.R,
                                          # branch.length=data$branch.length,
                                          # n.descen = data$n.descendants),
                                          spat_alg = "bootspat_str",
                                          spat_alg_args = list(rprob = NULL,
                                                               rich = NULL,
                                                               fr_prob = NULL),
                                          aleats = 5), "SpatRaster")
})

test_that("error is returned when the raster does not have a longitude/latitude
          coordinate reference system (CRS)", {

            x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
            tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
            # data <- phylo.pres(x, tree)
            # area.branch <- inv.range(data$x, data$branch.length)
            w <- terra::project(x, "EPSG:2169")

            data <- phyloraster::phylo.pres(w, tree)
            branch.length <- data$branch.length
            n.descen <- data$n.descendants
            area.branch <- phyloraster::inv.range(data$x, data$branch.length, LR = T)

            expect_error(geo.phylo.ses(x = data$x,
                                          tree = data$tree,
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

test_that("error is returned when only argument x is provided", {

            x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))

            expect_error(geo.phylo.ses(x = data$x,
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

test_that("error is returned when the user choose a randomization method not available", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
  # data <- phylo.pres(x, tree)
  # area.branch <- inv.range(data$x, data$branch.length)

  expect_error(geo.phylo.ses(x = data$x,
                             tree = data$tree,
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

test_that("function runs ok with the method 'tip'", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
  # data <- phylo.pres(x, tree)
  # area.branch <- inv.range(data$x, data$branch.length)

  expect(geo.phylo.ses(x = x,
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
