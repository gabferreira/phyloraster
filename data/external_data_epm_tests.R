# # generating data for testthat
# x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#
# # WE
# library(epm)
# datepm <- epm::createEPMgrid(x, resolution = 0.01)
# we <- epm::gridMetrics(datepm, metric = "weightedEndemism")
#
# # PE and PD
# tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
# data <- phylo.pres(x, tree)
#
# datepm <- epm::createEPMgrid(x, resolution = 0.01)
# data$tree <- as(data$tree, "phylo")
# datepm <- epm::addPhylo(datepm, data$tree)
# pe <- epm::gridMetrics(datepm, metric = "phyloWeightedEndemism")
# pd <- epm::gridMetrics(datepm, metric = "pd")
#
#
# ## saving in inst/extdata
# if(!dir.exists("inst/extdata")) dir.create("inst/extdata", recursive = T)
# terra::writeRaster(we$grid$weightedEndemism, "inst/extdata/epm_WE.tif")
# terra::writeRaster(pe$grid$phyloWeightedEndemism, "inst/extdata/epm_PE.tif")
# terra::writeRaster(pd$grid$pd, "inst/extdata/epm_PD.tif")
