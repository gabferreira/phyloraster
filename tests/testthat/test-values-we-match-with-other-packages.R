test_that("results of the analyses replicate those of other packages", {

  library(epm)
  library(phyloraster)
  # library(phyloregion)

  x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))

  # phyloraster
  pg <- phyloraster::rast.we(x)

  # # phyloregion
  # fdir <- system.file("extdata/rasters", package="phyloraster")
  # files <- file.path(fdir, dir(fdir))
  # com <- phyloregion::raster2comm(files)
  # pr <- phyloregion::weighted_endemism(com$comm_dat)
  # m <- merge(com$poly_shp, data.frame(grids=names(pr), WE=pr), by="grids")
  # m <- m[!is.na(m$WE),]
  # x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
  # r <- terra::rasterize(terra::vect(m), x, field = "WE")

  # epm
  datepm <- epm::createEPMgrid(x, resolution = 0.01)
  ep <- epm::gridMetrics(datepm, metric = "weightedEndemism")

  testthat::expect_equivalent(round(terra::values(pg), 10), round(terra::values(ep$grid$weightedEndemism)), 10)
  # testthat::expect_equivalent(round(values(pg), 10), round(values(r$WE)), 10)
})
