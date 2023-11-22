test_that("returned object classes are correct", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157,
                                   -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))
  data <- phylo.pres(x, tree)
  # branch.length <- data$branch.length
  # n.descen <- data$n.descendants
  inv.R <- phyloraster::inv.range(data$x)

  t <- geo.phylo(data$x, inv.R = inv.R, edge.path = data$edge.path,
                 branch.length = data$branch.length,
                 n.descen = data$n.descendants)

  # tests
  expect_s4_class(t, "SpatRaster")
})

test_that("Are the returned values correct?", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))

  # getting fewer cells to test all values
  xcrop <- terra::crop(x, terra::ext(c(150.0157, 150.8157,
                                       -23.044, -22.8563)))

  # matching phylogenetic tree and the raster
  data <- phylo.pres(xcrop, tree)
  branch.length <- data$branch.length
  n.descen <- data$n.descendants
  inv.R <- phyloraster::inv.range(data$x)

  # metric SE
  t <- geo.phylo(data$x, inv.R = inv.R, edge.path = data$edge.path,
                 branch.length = data$branch.length,
                 n.descen = data$n.descendants)

  se.observed <- terra::values(t$SR)
  se.expect <- matrix(data = c(12, 12, 12, 13, 14,
                               14, 14, 14, 12, 12,
                               11, 12, 13, 14, 14, 14),
                      ncol = 1, byrow = F)

  expect_equivalent(terra::values(t$SR), se.expect)

  # metric PE
  pe.obs <- round(terra::values(t$PE),7)
  pe.expect <- c(0.3807860, 0.3807860, 0.3807860,
                 0.4651102, 0.4776858, 0.4776858,
                 0.4776858, 0.4776858, 0.3805112,
                 0.3805112, 0.3496615, 0.3805112,
                 0.3930776, 0.4773410, 0.4773410,
                 0.4773410)
  expect_equivalent(pe.obs, pe.expect)

  # metric WE

  # metric PE
  we.obs <- round(terra::values(t$WE),7)
  we.expect <- c(0.7544373, 0.7544373, 0.7544373,
                 0.8794712, 1.0045163, 1.0045163,
                 1.0045163, 1.0045163, 0.7538928,
                 0.7538928, 0.6872518, 0.7538928,
                 0.8788477,
                 1.0037913, 1.0037913, 1.0037913)
  expect_equivalent(we.obs, we.expect)

  # metric PD
  pd.obs <- terra::values(t$PD)
  pd.expect <- matrix(data = c(6.583520, 6.583520, 6.583520,
                               6.705980, 7.704928,
                               7.704928, 7.704928, 7.704928,
                               6.583520, 6.583520,
                               5.994505, 6.583520, 7.582468,
                               7.704928, 7.704928,
                               7.704928),
                      ncol = 1, byrow = F)
  expect_equivalent(pd.obs, pd.expect)

  # metric ED
  ed.obs <- round(terra::values(t$ED), 7)
  ed.expect <- c(0.1530742, 0.1530742, 0.1530742,
                 0.1563603, 0.1735835,
                 0.1735835, 0.1735835, 0.1735835,
                 0.1530742, 0.1530742,
                 0.1371549, 0.1530742, 0.1702974,
                 0.1735835, 0.1735835, 0.1735835)
  expect_equivalent(ed.obs, ed.expect)
})

test_that("error is returned when an argument is missing", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))

  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157,
                                   -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))
  data <- phylo.pres(x, tree)
  inv.R <- phyloraster::inv.range(data$x)

  # tests
  expect_error(geo.phylo(data$x,
                         inv.R = inv.R,
                         branch.length = data$branch.length))

  expect_error(geo.phylo(data$x,
                         inv.R = inv.R,
                         n.descen = data$n.descendants))
})

test_that("arguments are calculated when is missing and the
          tree is provided", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))

  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157,
                                   -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))
  data <- phylo.pres(x, tree)
  # area.branch <- phyloraster::inv.range(data$x)

  # tests
  expect(geo.phylo(data$x, tree), ok = T)

})

test_that("names are reordened in the function geo.phylo", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157,
                                   -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))
  data <- phylo.pres(x, tree)
  inv.R <- phyloraster::inv.range(data$x)

  # tests
  expect(geo.phylo(data$x, tree, #range.BL = area.branch$range.BL,
                   inv.R = inv.R,
                   edge.path = data$edge.path[sample(1:nrow(data$edge.path)),],
                   branch.length = data$branch.length,
                   n.descen = data$n.descendants), ok = T)

})

test_that("error is returned when the raster does not have a longitude/latitude
          coordinate reference system (CRS)", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))

  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157,
                                   -23.044, -22.8563)))

  w <- terra::project(x, "EPSG:2169")

  data <- phyloraster::phylo.pres(w, tree)
  # branch.length <- data$branch.length
  # n.descen <- data$n.descendants
  # inv.R <- phyloraster::inv.range(data$x)

  # tests
  expect_error(geo.phylo(data$x,
                         data$n.descendants))
})

