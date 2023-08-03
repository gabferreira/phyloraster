test_that("returned object classes are correct", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
  data <- phylo.pres(x, tree)
  branch.length <- data$branch.length
  n.descen <- data$n.descendants
  inv.R <- phyloraster::inv.range(data$x)

  t <- geo.phylo(data$x, inv.R = inv.R,
                 branch.length = data$branch.length, n.descen = data$n.descendants)

  # tests
  expect_s4_class(t, "SpatRaster")
})

test_that("Are the returned values correct?", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))

  # getting fewer cells to test all values
  xcrop <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044, -22.8563)))

  # matching phylogenetic tree and the raster
  data <- phylo.pres(xcrop, tree)
  branch.length <- data$branch.length
  n.descen <- data$n.descendants
  inv.R <- phyloraster::inv.range(data$x)

  # metric SE
  t <- geo.phylo(data$x, inv.R = inv.R,
                 branch.length = data$branch.length, n.descen = data$n.descendants)

  se.observed <- terra::values(t$SR)
  se.expect <- matrix(data = c(12, 12, 12, 13, 14, 14, 14, 14, 12, 12,
                               11, 12, 13, 14, 14, 14),
                      ncol = 1, byrow = F)

  expect_equivalent(terra::values(t$SR), se.expect)

  # metric PE
  pe.obs <- round(terra::values(t$PE),7)

  pe.expect <- c(0.2277756,0.2277756,0.2277756,0.2345368,0.3594504,
              0.3594504,0.3594504,0.3594504,0.2276112,0.2276112,
    0.1883587,0.2276112,0.3524346,0.3591910,0.3591910,0.3591910)
  expect_equivalent(pe.obs, pe.expect)

  # metric WE

  # metric PE
  we.obs <- round(terra::values(t$WE),7)
  we.expect <- c(0.7544373, 0.7544373, 0.7544373, 0.8794712, 1.0045163, 1.0045163,
                 1.0045163, 1.0045163, 0.7538928,
                 0.7538928, 0.6872518, 0.7538928, 0.8788477,
                 1.0037913, 1.0037913, 1.0037913)
  expect_equivalent(we.obs, we.expect)

  # metric PD
  pd.obs <- terra::values(t$PD)
  pd.expect <- matrix(data = c(3.603842, 3.603842, 3.603842, 3.657917, 4.656865,
                               4.656865, 4.656865, 4.656865, 3.603842,
                               3.603842, 3.014827, 3.603842, 4.602790,
                               4.656865, 4.656865, 4.656865),
                      ncol = 1, byrow = F)
  expect_equivalent(pd.obs, pd.expect)

  # metric ED
  ed.obs <- terra::values(t$ED)
  ed.expect <- c(0.08343477, 0.08343477, 0.08343477, 0.08445505,
                               0.10715842, 0.10715842, 0.10715842,
                               0.10715842, 0.08343477, 0.08343477, 0.06751545,
                               0.08343477, 0.10613813, 0.10715842,
                               0.10715842, 0.10715842)
  expect_equivalent(ed.obs, ed.expect)
})

test_that("error is returned when an argument is missing", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))

  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
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

test_that("arguments are calculated when is missing and the tree is provided", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))

  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
  data <- phylo.pres(x, tree)
  # area.branch <- phyloraster::inv.range(data$x)

  # tests
  expect(geo.phylo(data$x, tree), ok = T)

})

test_that("names are reordened in the function geo.phylo", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
  data <- phylo.pres(x, tree)
  inv.R <- phyloraster::inv.range(data$x)

  # tests
  expect(geo.phylo(data$x, tree, #range.BL = area.branch$range.BL,
                   inv.R = inv.R,
                   branch.length = sort(data$branch.length),
                   n.descen = data$n.descendants), ok = T)

})

test_that("error is returned when the raster does not have a longitude/latitude coordinate reference system (CRS)", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))

  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157, -23.044, -22.8563)))

  w <- terra::project(x, "EPSG:2169")

  data <- phyloraster::phylo.pres(w, tree)
  branch.length <- data$branch.length
  n.descen <- data$n.descendants
  # inv.R <- phyloraster::inv.range(data$x)

  # tests
  expect_error(geo.phylo(data$x,
                         data$n.descendants))
})

