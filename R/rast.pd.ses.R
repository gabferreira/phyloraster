#' Phylogenetic diversity (PD. Faith, 1992) standardized for specie richness
#'
#' The function calculates the PD corrected for species richness. In each null
#' model run, richness is kept constant and branch lengths of each species are
#' randomized. The function provides the mean, standard deviation of all null
#' models and also calculates the standardized effect size (SES).
#'
#' @param x SpatRaster. A presence-absence SpatRaster.
#' @param branch.length numeric. A Named numeric vector containing the branch length of each specie. See phylo.pres function.
#' @param aleats numeric. A number indicating how many times the calculation should be repeated.
#' @param filename character. Output filename.
#' @param random character. A character indicating what type of randomization. Could be by tip, site, specie or full spat(site and specie).
#' @return SpatRaster
#' @export
#' @author Gabriela Alves-Ferreira and Neander Heming
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological conservation, 61(1), 1-10.
#' @example
#' \dontrun{
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
#' data <- phylogrid::phylo.pres(ras, tree)
#' t <- phylogrid::rast.pe.ses(ras, data$branch.length, aleats = 10, random = "fullspat")
#' plot(t)
#' }
#'
rast.pd.ses <- function(x, branch.length, aleats,
                        random = c("tip", "site", "specie", "fullspat"), filename = NULL){

  aleats <- aleats # number of null models
  temp <- vector("list", length = aleats) # to create a temporary vector with the raster number

  ## Null model (bootstrap structure)
  if(random == "tip"){

    pd.rand <- list() # to store the rasters in the loop
    pd.rand2 <- list() # list to store the rasters
    bl.random <- branch.length # to store the branch length in the loop

    for(i in 1:aleats){

      bl.random[] <- sample(branch.length) # randomize branch lengths

      temp[[i]] <- paste0(tempfile(), i, ".tif") # temporary names to rasters

      pd.rand[[i]] <- terra::app(x, fun = .vec.pd,
                                 branch.length = bl.random, filename = temp[[i]])
      ### selecionando apenas o primeiro raster (PD) e excluindo o segundo (PDR)
      pd.rand2[[i]] <- pd.rand[[i]][[1]] # only the first layer for each specie
    }

    pd.rand2 <- terra::rast(pd.rand2) # to transform a list in raster

  } else if (random == "site"){

    pd.rand <- list() # to store the rasters in the loop
    pd.rand2 <- list() # list to store the rasters

    for(i in 1:aleats){

      # temporary names to rasters
      temp[[i]] <- paste0(tempfile(), i, ".tif")

      ### embaralha por lyr - ordem dos sítios de cada espécie separada
      pres.site.null <- spat.rand(x, aleats = 1, random = "site")

      # calculate pd
      pd.rand[[i]] <- terra::app(pres.site.null, fun = .vec.pd,
                                 branch.length = branch.length, filename = temp[[i]])

      ### selecionando apenas o primeiro raster (PD) e excluindo o segundo (PDR)
      pd.rand2[[i]] <- pd.rand[[i]][[1]] # only the first layer for each specie
    }

    pd.rand2 <- terra::rast(pd.rand2) # to transform a list in raster

  } else if (random == "specie") {

    ### randomize by cells - species in each site
    pd.rand <- list() # to store the rasters in the loop
    pd.rand2 <- list() # list to store the rasters

    for(i in 1:aleats){

      temp[[i]] <- paste0(tempfile(), i, ".tif") # temporary names to rasters
      sp.rand <- spat.rand(x, aleats = 1, random = "specie")
      pd.rand[[i]] <- terra::app(sp.rand,
                                 fun = .vec.pd,
                                 branch.length = branch.length, filename = temp[[i]])

      ### selecionando apenas o primeiro raster (PD) e excluindo o segundo (PDR)
      pd.rand2[[i]] <- pd.rand[[i]][[1]] # only the first layer for each specie
    }

    pd.rand2 <- terra::rast(pd.rand2) # to transform a list in raster

  } else if (random == "fullspat") {

    pd.rand <- list() # to store the rasters in the loop
    pd.rand2 <- list() # list to store the rasters
    fr <- terra::freq(x)

    for(i in 1:aleats){

      temp[[i]] <- paste0(tempfile(), i, ".tif") # temporary names to rasters

      ### randomize sites and species
      pres.null <- terra::app(x, fun=.lyr.sample, fr=fr)

      pd.rand[[i]] <- terra::app(pres.null, fun = .vec.pd,
                                 branch.length = branch.length, filename = temp[[i]])

      ### selecionando apenas o primeiro raster (PD) e excluindo o segundo (PDR)
      pd.rand2[[i]] <- pd.rand[[i]][[1]] # only the first layer for each specie
    }

    pd.rand2 <- terra::rast(pd.rand2) # to transform a list in raster

  } else {
    stop("Choose a valid randomization method! The methods currently available are: 'site', 'specie', 'fullspat'.")
  }

  ### PD rand mean
  pd.rand.mean <- terra::mean(pd.rand2, na.rm = TRUE) # mean pd
  ### PD rand SD
  pd.rand.sd <- terra::stdev(pd.rand2, na.rm = TRUE) # sd pd

  ### PD observed
  {
    pd.obs <- phylogrid::rast.pd(x, branch.length)
    pd.obs <- pd.obs$PD # selecting only PD
  }

  ### Concatenate rasters
  pd <- c(pd.obs, pd.rand.mean, pd.rand.sd)

  ## Calculating the standard effect size (SES)
  {
    ses <- function(x){
      (x[1] - x[2])/sqrt(x[3])
    }
    pd.ses <- terra::app(c(pd.obs, pd.rand.mean, pd.rand.sd),
                         ses)
    names(pd.ses) <- "SES"
  }

  out <- c(pd, pd.ses)
  names(out) <- c("PD Observed", "Mean", "SD", "SES")

  if (!is.null(filename)){ # to save the rasters when the path is provide
    out <- terra::writeRaster(out, filename, ...)
  }

  unlink(temp) # delete the archive that will not be used

  return(out)
}

#' Phylogenetic diversity (PD. Faith, 1992) standardized for specie richness- version 2
#'
#' The function calculates the PD corrected for species richness. In each null
#' model run, richness is kept constant and branch lengths of each species are
#' randomized. The function provides the mean, standard deviation of all null
#' models and also calculates the standardized effect size (SES).
#'
#' @param x SpatRaster. A presence-absence SpatRaster with the layers ordered according to the branch.length order. See phylo.pres function.
#' @param branch.length numeric. A Named numeric vector containing the branch length of each specie. See phylo.pres function.
#' @param aleats numeric. A number indicating how many times the calculation should be repeated.
#' @param filename character. Output filename.
#' @param random character. A character indicating what type of randomization. Could be by tip, site, specie or full spat(site and specie).
#' @return SpatRaster
# #' @export
#' @author Neander Heming and Gabriela Alves-Ferreira
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological conservation, 61(1), 1-10.
# #' @examples
rast.pd.ses.v2 <- function(x, branch.length, aleats,
                           random = c("tip", "site", "specie", "fullspat"), filename = NULL){
  {
    aleats <- aleats # number of null models
    temp <- vector("list", length = aleats) # to create a temporary vector with the raster number
  }

  ### Função para reamostrar de acordo com a frequência observada
  fr <- terra::freq(x)

  lyr.sample <- function(x, fr){
    sapply(x, function(x, fr){
      if(is.na(x)){
        return(NA)
      } else {
        return(sample(fr$value, 1, prob = fr$count))
      }
    }, fr = fr)
  }

  ## Null model (bootstrap structure)
  if(random == "tip"){

    pd.rand <- list() # to store the rasters in the loop
    bl.random <- branch.length # to store the branch length in the loop
    for(i in 1:aleats){
      bl.random[] <- sample(branch.length) # randomize branch lengths
      ## check if the values are different
      # branch.length == bl.random

      temp[[i]] <- paste0(tempfile(), i, ".tif") # temporary names to rasters

      pd.rand[[i]] <- terra::app(x, fun = .vec.pd,
                                 branch.length = bl.random, filename = temp[[i]])
    }

    ### selecionando apenas o primeiro raster (PD) e excluindo o segundo (PDR)
    pd.rand2 <- list() # list to store the rasters
    for(j in 1:length(pd.rand)){
      pd.rand2[[j]] <- pd.rand[[j]][[1]] # only the first layer for each specie
    }
    pd.rand2 <- terra::rast(pd.rand2) # to transform a list in raster
  } else if (random == "site"){

    pd.rand <- list() # to store the rasters in the loop
    pd.rand2 <- list() # list to store the rasters
    for(i in 1:aleats){
      ### embaralha por lyr - ordem dos sítios de cada espécie separada

      temp[[i]] <- paste0(tempfile(), i, ".tif") # temporary names to rasters
      pres.site.null <- terra::rast(lapply(1:terra::nlyr(x),
                                           function(i, r, fr){
                                             terra::app(r[[i]], fun=lyr.sample, fr=fr[fr$layer==i,])
                                           }, r = x, fr = fr))

      pd.rand[[i]] <- terra::app(pres.site.null, fun = .vec.pd,
                                 branch.length = branch.length, filename = temp[[i]])
      ### selecionando apenas o primeiro raster (PD) e excluindo o segundo (PDR)
      pd.rand2[[i]] <- pd.rand[[i]][[1]] # only the first layer for each specie
    }
    pd.rand2 <- terra::rast(pd.rand2) # to transform a list in raster

  } else if (random == "specie") {

    ### embaralha por célula - espécies em cada sítio
    pd.rand <- list() # to store the rasters in the loop
    pd.rand2 <- list() # list to store the rasters
    for(i in 1:aleats){
      temp[[i]] <- paste0(tempfile(), i, ".tif") # temporary names to rasters

      pd.rand[[i]] <- terra::app(terra::app(x, sample),
                                 fun = .vec.pd,
                                 branch.length = branch.length, filename = temp[[i]])
      ### selecionando apenas o primeiro raster (PD) e excluindo o segundo (PDR)
      pd.rand2[[i]] <- pd.rand[[i]][[1]] # only the first layer for each specie
    }
    pd.rand2 <- terra::rast(pd.rand2) # to transform a list in raster

  } else {
    pd.rand <- list() # to store the rasters in the loop
    pd.rand2 <- list() # list to store the rasters
    for(i in 1:aleats){
      ### embaralha tudo!
      pres.all.null <- terra::app(x, fun=lyr.sample, fr=fr)
      temp[[i]] <- paste0(tempfile(), i, ".tif") # temporary names to rasters

      pd.rand[[i]] <- terra::app(pres.all.null, fun = .vec.pd,
                                 branch.length = branch.length, filename = temp[[i]])
      ### selecionando apenas o primeiro raster (PD) e excluindo o segundo (PDR)
      pd.rand2[[i]] <- pd.rand[[i]][[1]] # only the first layer for each specie
    }
    pd.rand2 <- terra::rast(pd.rand2) # to transform a list in raster
  }

  ### PD rand mean
  pd.rand.mean <- terra::mean(pd.rand2, na.rm = TRUE, filename = filename) # mean pd

  ### PD rand SD
  pd.rand.sd <- terra::stdev(pd.rand2, na.rm = TRUE, filename = filename) # sd pd

  ### PD observed
  {
    pd.obs <- phylogrid::rast.pd(x, branch.length, filename = filename)
    pd.obs <- pd.obs$PD # selecting only PD
  }

  ### Concatenate rasters
  pd <- c(pd.obs, pd.rand.mean, pd.rand.sd)

  ## Calculating the standard effect size (SES)
  {
    ses <- function(x){
      (x[1] - x[2])/sqrt(x[3])
    }
    pd.ses <- terra::app(c(pd.obs, pd.rand.mean, pd.rand.sd),
                         ses)
    names(pd.ses) <- "SES"
  }

  out <- c(pd, pd.ses)
  names(out) <- c("PD Observed", "Mean", "SD", "SES")

  if (!is.null(filename)){ # to save the rasters when the path is provide
    out <- terra::writeRaster(out, filename, ...)
  }

  unlink(temp) # delete the archive that will not be used

  return(out)
}
