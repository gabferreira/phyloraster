test_that("check if the returned object class is correct", {

  # load data
  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
  pp <- phylo.pres(ras, tree)

  # tests
  expect_s4_class(pp$x, "SpatRaster")
  expect_type(pp$branch.length, "double")
})

test_that("Test that error is returned with wrong input class", {

  # load data and functions
  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
  pp <- phylo.pres(ras, tree)

  # tests
  expect_error(phylo.pres(ras, pp$branch.length))
})

test_that("Test that error when species names do not match between the raster and the tree", {

  # load data and functions
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))

  names(x) <- c( "1Litoria_revelata" ,  "1Litoria_rothii" ,  "1Litoria_longirostris" ,
                 "1Litoria_dorsalis" ,  "1Litoria_rubella" , "1Litoria_nigrofrenata" ,
                 "1Litoria_nasuta",  "1Litoria_tornieri",  "1Litoria_inermis" ,
                 "1Litoria_pallida",   "1Litoria_latopalmata" ,
                 "1Litoria_bicolor","1Litoria_fallax" ,"1Litoria_genimaculata",
                 "1Litoria_andiirrmalin" , "1Litoria_wilcoxii" , "1Litoria_jungguy"  ,
                 "1Litoria_caerulea"  , "1Litoria_gracilenta" ,       "1Litoria_chloris" ,
                 "1Litoria_xanthomera"  ,      "1Cyclorana_brevipes"  ,
                 "1Cyclorana_novaehollandiae", "1Cyclorana_manya"   ,
                 "1Cyclorana_cultripes" ,      "1Litoria_alboguttata" ,
                 "1Cyclorana_longipes"    ,    "1Nyctimystes_dayi"     ,
                 "1Litoria_nannotis"    ,      "1Litoria_lorica"    ,
                 "1Litoria_rheocola"  ,        "1Litoria_nyakalensis" ,
                 "1Litoria_infrafrenata")

  # tests
  expect_error(phylo.pres(x, tree))
})

test_that("Are the returned values correct?", {

  # load data
  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))

  # getting fewer cells to test all values
  r <- terra::rast()
  terra::ext(r) <- c(150.0157, 150.8157, -23.044, -22.8563)
  xcrop <- terra::crop(x, terra::ext(r))

  pp.obs <- suppressWarnings(phylo.pres(xcrop[[1:10]], tree, pruning = "subtree"))
  descen.expect <- c(12, 12, 13, 13, 12, 14, 15, 16, 16, 15)
  # (terra::values(pp.obs[[1]]))[,10]
  rast.expect <- matrix(data = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))

  branch.expect <- c(0.589016, 0.589015, 0.395144, 0.395144, 0.589015,
                     0.288907, 0.196139, 0.129514, 0.129514, 0.196139)
  # test
  expect_equivalent(pp.obs$n.descendants, descen.expect)
  expect_equivalent(terra::values(pp.obs$x), rast.expect)
  expect_equivalent(pp.obs$branch.length, branch.expect)
})

test_that("names in the raster and in the branch lenght are in the same order", {

  # load data
  ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
  datapres <- phylo.pres(ras, tree, pruning = "subtree")

  # test
  expect_equal(all.equal(names(datapres$x), tree$tip.label), TRUE)
})

test_that("error is returned when x class is different of SpatRaster", {

  # load data
  require(raster)
  ras <- raster::brick(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))

  # tests
  expect_error(phylo.pres(ras, tree))

})
