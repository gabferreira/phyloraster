#' Phylogenetic diversity (PD. Faith, 1992) standardized for specie richness
#'
#' The function calculates the PD corrected for species richness. In each null
#' model run, richness is kept constant and branch lengths of each species are
#' randomized. The function provides the mean, standard deviation of all null
#' models and also calculates the standardized effect size (SES).
#'
#' @param pres.reord SpatRaster. A presence-absence SpatRaster with the layers ordered according to the branch.length order. See phylo.pres function.
#' @param branch.length numeric. A Named numeric vector containing the branch length of each specie. See phylo.pres function.
#' @param aleats numeric. A number indicating how many times the calculation should be repeated.
#' @return SpatRaster
#' @export
#' @author Neander Heming and Gabriela Alves-Ferreira
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological conservation, 61(1), 1-10.
# #' @examples
rast.pd.ses <- function(pres.reord, branch.length, aleats, filename = NULL){
  {
    aleats <- aleats # number of null models
    temp <- vector("list", length = aleats) # to create a temporary vector with the raster number
  }

  ## Null model (bootstrap structure)
  {
    pd.rand <- list() # to store the rasters in the loop
    bl.random <- branch.length # to store the branch length in the loop
    for(i in 1:aleats){
      bl.random[] <- sample(branch.length) # randomize branch lengths
      ## check if the values are different
      # branch.length == bl.random

      temp[[i]] <- paste0(tempfile(), i, ".tif") # temporary names to rasters

      pd.rand[[i]] <- terra::app(pres.reord, fun = .vec.pd,
                                 branch.length = bl.random, filename = temp[[i]])
    }

    ### selecionando apenas o primeiro raster (PD) e excluindo o segundo (PDR)
    pd.rand2 <- list() # list to store the rasters
    for(j in 1:length(pd.rand)){
      pd.rand2[[j]] <- pd.rand[[j]][[1]] # only the first layer for each specie
    }
    pd.rand2 <- terra::rast(pd.rand2) # to transform a list in raster
  }

  ### PD rand mean
  pd.rand.mean <- terra::mean(pd.rand2, na.rm = TRUE) # mean pd

  ### PD rand SD
  pd.rand.sd <- terra::stdev(pd.rand2, na.rm = TRUE) # sd pd

  ### PD observed
  {
    pd.obs <- phylogrid::rast.pd(pres.reord = pres.reord,
                                 branch.length = branch.length)
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

pd.ses <- rast.pd.ses(teste$pres.reord, teste$branch.length, aleats = 10)
plot(pd.ses)




rast.pd.ses2 <- function(pres.reord, branch.length, aleats,
                         random = c("tip", "site", "specie", "all"), filename = NULL){
  {
    aleats <- aleats # number of null models
    temp <- vector("list", length = aleats) # to create a temporary vector with the raster number
  }

  ### Função para reamostrar de acordo com a frequência observada
  fr <- freq(pres.reord)
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

        pd.rand[[i]] <- terra::app(pres.reord, fun = .vec.pd,
                                   branch.length = bl.random, filename = temp[[i]])
      }

      ### selecionando apenas o primeiro raster (PD) e excluindo o segundo (PDR)
      pd.rand2 <- list() # list to store the rasters
      for(j in 1:length(pd.rand)){
        pd.rand2[[j]] <- pd.rand[[j]][[1]] # only the first layer for each specie
      }
      pd.rand2 <- terra::rast(pd.rand2) # to transform a list in raster
  } else if (random == "site"){

    ### embaralha por lyr - ordem dos sítios de cada espécie separada
    pres.site.null <- rast(lapply(1:nlyr(pres.reord),
                           function(i, r, fr){
                             terra::app(r[[i]], fun=lyr.sample, fr=fr[fr$layer==i,])
                           }, r = pres.reord, fr = fr))

    pd.rand <- list() # to store the rasters in the loop
    for(l in 1:aleats){
      temp[[l]] <- paste0(tempfile(), l, ".tif") # temporary names to rasters

      pd.rand[[l]] <- terra::app(pres.site.null, fun = .vec.pd,
                                 branch.length = branch.length, filename = temp[[l]])
    }

    ### selecionando apenas o primeiro raster (PD) e excluindo o segundo (PDR)
    pd.rand2 <- list() # list to store the rasters
    for(m in 1:length(pd.rand)){
      pd.rand2[[m]] <- pd.rand[[m]][[1]] # only the first layer for each specie
    }
    pd.rand2 <- terra::rast(pd.rand2) # to transform a list in raster
  } else if (random == "specie") {

    ### embaralha por célula - espécies em cada sítio
    pres.sp.null <- terra::app(pres.reord, sample)

    pd.rand <- list() # to store the rasters in the loop
    for(n in 1:aleats){
      temp[[n]] <- paste0(tempfile(), n, ".tif") # temporary names to rasters

      pd.rand[[n]] <- terra::app(pres.site.null, fun = .vec.pd,
                                 branch.length = branch.length, filename = temp[[n]])
    }

    ### selecionando apenas o primeiro raster (PD) e excluindo o segundo (PDR)
    pd.rand2 <- list() # list to store the rasters
    for(o in 1:length(pd.rand)){
      pd.rand2[[o]] <- pd.rand[[o]][[1]] # only the first layer for each specie
    }
    pd.rand2 <- terra::rast(pd.rand2) # to transform a list in raster
  } else (random == "all"){
    ### embaralha tudo!
    r4.null <- terra::app(r2, fun=lyr.sample, fr=fr2)

    pd.rand <- list() # to store the rasters in the loop
    for(p in 1:aleats){
      temp[[p]] <- paste0(tempfile(), p, ".tif") # temporary names to rasters

      pd.rand[[p]] <- terra::app(pres.site.null, fun = .vec.pd,
                                 branch.length = branch.length, filename = temp[[p]])
    }

    ### selecionando apenas o primeiro raster (PD) e excluindo o segundo (PDR)
    pd.rand2 <- list() # list to store the rasters
    for(q in 1:length(pd.rand)){
      pd.rand2[[q]] <- pd.rand[[q]][[1]] # only the first layer for each specie
    }
    pd.rand2 <- terra::rast(pd.rand2) # to transform a list in raster
  }


    ### PD rand mean
    pd.rand.mean <- terra::mean(pd.rand2, na.rm = TRUE) # mean pd

    ### PD rand SD
    pd.rand.sd <- terra::stdev(pd.rand2, na.rm = TRUE) # sd pd

    ### PD observed
    {
      pd.obs <- phylogrid::rast.pd(pres.reord = pres.reord,
                                   branch.length = branch.length)
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

