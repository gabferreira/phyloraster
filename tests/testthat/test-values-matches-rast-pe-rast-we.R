test_that("PE function returns the same result as WE when the tree is a rake", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))

  # setting all branches to length equal to one
  tree2 <- glottoTrees::rescale_branches(tree, length = 1)
  data <- phylogrid::phylo.pres(x, tree2)

  # calculating WE and PE with all branch lengths = 1
  we <- phylogrid::rast.we(x)
  pe <- phylogrid::rast.pe(data$x, data$branch.length)

  # testing if the values are the same for WE and PE
  expect_equivalent(terra::values(we), terra::values(pe))

})
