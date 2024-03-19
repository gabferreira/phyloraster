test_that("tip lengths are correctly calculated", {
  library(ape)
  set.seed(1)
  tree <- rtree(n=5)

  bls <- setNames(c(0.0617863, 0.6870228, 0.7698414, 0.4976992, 0.7176185 ),
                  paste0("t", c(2,1,3,4,5)))

  edge.info <- tip.root.path(tree)

  expect_equal(round(species.tip.length(tree), 7),
               bls)
  expect_equal(round(species.tip.length(edge.info = edge.info), 7),
               bls)

})
