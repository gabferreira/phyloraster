test_that("check if the returned object class is correct", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157,
                                   -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))
  pp <- phylo.pres(x, tree, full_tree_metr = F)
  ppF <- phylo.pres(x, tree, full_tree_metr = T)

  # tests
  expect_s4_class(pp$x, "SpatRaster")
  expect_type(pp$branch.length, "double")
  expect_s4_class(ppF$x, "SpatRaster")
  expect_type(ppF$branch.length, "double")
})

test_that("Test that error is returned with wrong input class", {

  # load data and functions
  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))

  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157,
                                   -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))
  pp <- phylo.pres(x, tree)

  # tests
  expect_error(phylo.pres(x, pp$branch.length))
})

test_that("Test that error when species names do not match between
          the raster and the tree", {

  # load data and functions
  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))

  names(x) <- c( "1Litoria_revelata", "1Litoria_rothii",
                 "1Litoria_longirostris",
                 "1Litoria_dorsalis", "1Litoria_rubella",
                 "1Litoria_nigrofrenata",
                 "1Litoria_nasuta",  "1Litoria_tornieri",
                 "1Litoria_inermis",
                 "1Litoria_pallida", "1Litoria_latopalmata",
                 "1Litoria_bicolor", "1Litoria_fallax",
                 "1Litoria_genimaculata",
                 "1Litoria_andiirrmalin", "1Litoria_wilcoxii",
                 "1Litoria_jungguy",
                 "1Litoria_caerulea", "1Litoria_gracilenta",
                 "1Litoria_chloris",
                 "1Litoria_xanthomera", "1Cyclorana_brevipes",
                 "1Cyclorana_novaehollandiae", "1Cyclorana_manya",
                 "1Cyclorana_cultripes", "1Litoria_alboguttata",
                 "1Cyclorana_longipes", "1Nyctimystes_dayi",
                 "1Litoria_nannotis", "1Litoria_lorica",
                 "1Litoria_rheocola", "1Litoria_nyakalensis",
                 "1Litoria_infrafrenata")

  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157,
                                   -23.044, -22.8563)))


  # tests
  expect_error(phylo.pres(x, tree))
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

  pp.obs <- suppressWarnings(phylo.pres(xcrop[[1:15]], tree,
                                        full_tree_metr = T))
  descen.expect <- c(0,  1,  1,  1,  1,  2,  1,  5,  1,  1,  1,  1,  1,  2,  1,
                     4,  5, 6, 11,  1,  1,  2, 13,  1,  1,  2)

  branch.expect <- c(0.173157, 0.072386, 0.175442, 0.589016, 0.589015, 0.193871,
                     0.395144, 0.395144, 0.589015, 0.273621, 0.490836, 0.201929,
                     0.288907, 0.092767, 0.196139, 0.066625, 0.129514, 0.129514,
                     0.196139, 0.325841, 0.511002, 0.511002, 0.011052, 0.095197,
                     0.398434, 0.020969, 0.017805, 0.466543, 0.043192, 0.423350,
                     0.322782, 0.100568, 0.100568, 0.093810, 0.390538, 0.268077,
                     0.122460, 0.068385, 0.054075, 0.054075, 0.251337, 0.253980,
                     0.253980, 0.043002, 0.210977, 0.210978, 0.032253, 0.178724,
                     0.178724, 0.154130, 0.749621, 0.401051, 0.348570, 0.348570,
                     0.081570, 0.267000, 0.267000, 0.998948)

  branch.alt.expect <- c(rep(0.01724138, 58))

  # test
  expect_equivalent(pp.obs$n.descendants, descen.expect)
  expect_equivalent(pp.obs$branch.length, branch.expect)
  expect_equivalent(pp.obs$branch.length.alt, branch.alt.expect)

})

test_that("names in the raster and in the node path matrix
          are in the same order", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif",
                               package="phyloraster"))
  # getting fewer cells to test all values
  x <- terra::crop(x, terra::ext(c(150.0157, 150.8157,
                                   -23.044, -22.8563)))

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))
  datapres <- phylo.pres(x, tree, full_tree_metr = FALSE)

  # test
  expect_equal(all.equal(rownames(datapres$edge.path),
                         tree$tip.label), TRUE)
})

test_that("error is returned when x class is
          different of SpatRaster", {

  # load data
  x <- matrix(nrow = 5, ncol = 5)

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))

  # tests
  expect_error(phylo.pres(x, tree))

})

