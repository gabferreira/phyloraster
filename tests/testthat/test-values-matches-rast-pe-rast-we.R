test_that("PE function returns the same result as WE when the tree is a rake", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))

  # setting all branches to length equal to one
  tree2 <- glottoTrees::rescale_branches(tree, length = 1)
  data <- phyloraster::phylo.pres(x, tree2)

  # calculating WE and PE with all branch lengths = 1
  we <- phyloraster::rast.we(x)
  pe <- phyloraster::rast.pe(data$x, data$branch.length)

  # testing if the values are the same for WE and PE
  expect_equivalent(terra::values(we), terra::values(pe))

})
