test_that("check if the returned object class is correct", {

  # load data
  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))

  ed <- species.ed(tree)

  # tests
  expect_s3_class(ed, "data.frame")
})

test_that("Are the returned values of species.ed correct?", {

  tree <- ape::read.tree(system.file("extdata", "tree.nex",
                                     package="phyloraster"))

  # values of ED from the package picante
  # ed.picante <- picante::evol.distinct(tree, type = "fair.proportion")
  ed.picante <- c(0.6440047, 0.6440037, 0.5470682, 0.5470682,
                  0.6440037, 0.5563398, 0.3947966, 0.3252204,
                  0.2919079, 0.2919079, 0.3252204, 0.6872423, 0.6872423,
                  0.5076379, 0.4788422, 0.3174512,
                  0.3174512, 0.4506342, 0.2719152, 0.2377227, 0.2377227,
                  0.3298920, 0.3298920, 0.2976395,
                  0.2976405, 0.2815130, 0.2815130, 0.7860100, 0.4852217,
                  0.4852217, 0.4444367, 0.4444367, 0.9995006)
  ed.phyloraster <- phyloraster::species.ed(tree)

  # metric ED
  expect_equivalent(round(ed.picante, 7), round(ed.phyloraster$ED, 7))

})
