#' Phylogenetic diversity (PD. Faith, 1992) standardized for specie richness
#'
#' The function calculates the PD corrected for species richness. In each null model run, richness is kept constant and branch lengths of each species are randomized. The function provides the mean, standard deviation of all null models and also calculates the standardized effect size (SES).
#'
#' @param pres.reord SpatRaster. A presence-absence SpatRaster with the layers ordered according to the branch.length order. See phylo.pres function.
#' @param branch.length numeric. A Named numeric vector containing the branch length of each specie. See phylo.pres function.
#' @param aleats numeric. A number indicating how many times the calculation should be repeated.
#' @return SpatRaster
#' @export
#' @references Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological conservation, 61(1), 1-10.
# #' @examples
rast.pd.ses <- function(pres.reord, branch.length, aleats){
  {
    aleats <- aleats # number of null models
    temp <- vector("list", length = aleats) # to create a temporary vector with the raster number
    temp2 <- vector("list", length = 4) # to create a temporary vector
    temp2[[1]] <- paste0(tempfile(), "pd.obs.tif") # to save the pd obs
    temp2[[2]] <- paste0(tempfile(), "pd.rand.mean.tif") # to save the mean pd
    temp2[[3]] <- paste0(tempfile(), "pd.rand.sd.tif") # to save the sd pd
    temp2[[4]] <- paste0(tempfile(), "pd.ses.tif") # to save the pd ses
  }

  ## Null model (bootstrap structure)
  {
    pd.rand <- list() # to store the rasters in the loop

    for(i in 1:aleats){
      bl.random <- sample(branch.length, replace = T) # aleatorize the branch lenght
      ## check if the values are differents
      # branch.length == bl.random

      for(j in 1:aleats){
        temp[[j]] <- paste0(tempfile(), j, ".tif") # directory to store the rasters
      }
      pd.rand[[i]] <- terra::app(pres.reord, fun = .vec.pd,
                                 branch.length = bl.random, filename = temp[[i]])
    }
    pd.rand <- terra::rast(pd.rand) # to transform a list in raster
  }

  ## PD observed
  pd.obs <- phylogrid::rast.pd(pres.reord = pres.reord,
                               branch.length = branch.length, filename = temp2[[1]])
  pd.obs <- pd.obs$PD # selecting only PD

  ## PD rand mean and PD rand SD
  pd.rand.mean <- terra::mean(pd.rand, na.rm = TRUE, overwrite = TRUE, filename = temp2[[2]]) # mean pd
  pd.rand.sd <- terra::stdev(pd.rand, na.rm = TRUE, overwrite = TRUE, filename = temp2[[3]]) # sd pd

  unlink(temp) # delete the archive that will not be used
  # unlink(temp2) # delete the archive that will not be used

  ## Calculating the standard effect size (SES)
  ses <- function(x, y, z){
    (x - y)/z
  }
  pd.ses <- ses(x = pd.obs, y = pd.rand.mean, z = pd.rand.sd)
  names(pd.ses) <- "SES"

  out <- c(pd.obs, pd.rand.mean, pd.rand.sd, pd.ses)
  names(out) <- c("PD Observed", "Mean", "SD", "SES" )
  return(out)
}

# adicionar o PVALUE
### pd_obs_p: P-value (quantile) of observed PD vs. null communities \eqn{= mpd_obs_rank / iter + 1}
