test_that("returned object classes are correct", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157,
                                   -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))

  # tests
  expect_s4_class(rast.ed(x, tree), "SpatRaster")

})

test_that("error is returned when the argument n.descendents
          is missing", {

            x <- terra::rast(system.file("extdata", "rast.presab.tif",
                                         package="phyloraster"))
            # getting fewer cells to test all values
            x <- terra::crop(x, terra::ext(c(150.0157, 150.8157,
                                             -23.044, -22.8563)))

            tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                               package="phyloraster"))
            data <- phylo.pres(x, tree)

            # metric PE
            expect_error(rast.ed(data$x, branch.length = data$branch.length))

          })

test_that("error is returned when the raster does not have a
          longitude/latitude coordinate reference system (CRS)", {

            x <- terra::rast(system.file("extdata", "rast.presab.tif",
                                         package="phyloraster"))

            # getting fewer cells to test all values
            x <- terra::crop(x, terra::ext(c(150.0157, 150.8157,
                                             -23.044, -22.8563)))

            tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                               package="phyloraster"))

            w <- terra::project(x, "EPSG:2169")

            data <- phylo.pres(w, tree)

            # metric PE
            expect_error(rast.ed(data$x, data$branch.length,
                                 data$n.decendents))

          })


