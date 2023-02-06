test_that("Are the returned values correct?", {

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))

  # metric SE
  se <- terra::values(phylogrid::geo.phylo(x, metric = "richness"))[2995:2999,]
  expect_equivalent(se, c(4, 3, 2, 6, 7))

  # metric PE
  pe <- terra::values(phylogrid::geo.phylo(x, tree, metric = "phylo.endemism"))[2995:2999,]
  expect_equivalent(pe, c(7.026690e-12, 2.603572e-12, 1.183266e-12, 1.273825e-11,
                          1.142087e-11))

  # metric PD
  pd <- terra::values(phylogrid::geo.phylo(x, tree, metric = "phylo.diversity"))[2995:2999,]
  expect_equivalent(pd, c(0.985009, 0.774032, 0.383494, 2.163039, 1.892623))

  # metric ED
  ed <- terra::values(phylogrid::geo.phylo(x, tree, metric = "evol.distinct"))[2995:2999,]
  expect_equivalent(ed, c(0.01928055, 0.01544461, 0.00778700, 0.05111920, 0.04233876))

})

