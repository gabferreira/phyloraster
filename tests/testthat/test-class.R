test_that("returned object classes are correct", {

  # load data and functions
  phylo.pres <- phylogrid::phylo.pres
  inv.range <- phylogrid::inv.range
  geo.phylo <- phylogrid::geo.phylo

  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
  datapres <- phylo.pres(ras, tree)
  range <- inv.range(datapres$x, datapres$branch.length)
  gp <- geo.phylo(x = datapres$x, area.inv = range$inv.R,
                  area.tips = range$LR, branch.length = datapres$branch.length, filename = NULL)

  # tests
  expect_s4_class(gp$PD, "SpatRaster")
  expect_s4_class(gp$PDR, "SpatRaster")
  expect_s4_class(gp$WE, "SpatRaster")
  expect_s4_class(gp$PE, "SpatRaster")
  expect_s4_class(range$inv.R, "SpatRaster")
  expect_s4_class(range$LR, "SpatRaster")
})

