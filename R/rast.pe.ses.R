#' Calculate phylogenetic endemism (PE. Rosauer et al. 2009) for a raster
#'
#' Calculate phylogenetic endemism following Rosauer et al. (2009) using rasters as input and output.
#'
#' @param pres.reord SpatRaster. A presence-absence SpatRaster with the layers ordered according to the tree order
#' @param branch.length numeric. A Named numeric vector containing the branch length of each specie
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @return SpatRaster
# #' @export
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C., & Cook, L. G. (2009). Phylogenetic endemism: a new approach for identifying geographical concentrations of evolutionary history. Molecular ecology, 18(19), 4061-4072.
#' @examples
.rast.pe.B <- function(pres.reord, branch.length, filename = NULL, ...){
  {
    # criei a função .rast.pe.B que não retorna mensagem de erro quando
    # nomes sao diferentes, ao contrario da rast.pe. Assim é possivel rodar o
    # modelo nulo sem gerar erro pq os nomes sao diferentes no bl e pres.reord.

    area.branch <- phylogrid::inv.range(pres.reord = pres.reord,
                                        branch.length = branch.length)

    rpe <- terra::app(area.branch$LR,
                      function(x){
                        if(all(is.na(x))){
                          return(NA)
                        }
                        sum(x, na.rm = T)
                      })
    rpe <- terra::app(rpe, function(x, m){ # to reescale values from 0 to 1
      (x/m)
    }, m = terra::minmax(rpe)[2,])
  }

  names(rpe) <- c("PE")

if (!is.null(filename)){ # to save the rasters when the path is provide
  rpe <- terra::writeRaster(rpe, filename, ...)
}
return(rpe)
}

#' Phylogenetic endemism (PE. Rosauer et al. 2009) standardized for specie richness
#'
#' The function calculates the PE corrected for species richness. In each null model run, richness is kept constant and branch lengths of each species are randomized. The function provides the mean and standard deviation of all null models and also calculates the standardized effect size (SES).
#'
#' @param pres.reord SpatRaster. A presence-absence SpatRaster with the layers ordered according to the branch.length order. See phylo.pres function.
#' @param branch.length numeric. A Named numeric vector containing the branch length of each specie. See phylo.pres function.
#' @param aleats numeric. A number indicating how many times the calculation should be repeated.
#' @return SpatRaster
#' @export
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C., & Cook, L. G. (2009). Phylogenetic endemism: a new approach for identifying geographical concentrations of evolutionary history. Molecular ecology, 18(19), 4061-4072.
# #' @examples
rast.pe.ses <- function(pres.reord, branch.length, aleats){
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
    pe.rand <- list() # to store the rasters in the loop

    for(i in 1:aleats){
      bl.random <- sample(branch.length, replace = T) # aleatorize the branch lenght
      ## check if the values are differents
      # branch.length == bl.random

      for(j in 1:aleats){
        temp[[j]] <- paste0(tempfile(), j, ".tif") # directory to store the rasters
      }
      pe.rand[[i]] <- .rast.pe.B(pres.reord, branch.length = bl.random,
                                         filename = temp[[i]])
    }
    pe.rand <- terra::rast(pe.rand) # to transform a list in raster
  }

  ## PD observed
  pe.obs <- phylogrid::rast.pe(pres.reord = pres.reord,
                               branch.length = branch.length, filename = temp2[[1]])

  ## PD rand mean and PD rand SD
  pe.rand.mean <- terra::mean(pe.rand, na.rm = TRUE, overwrite = TRUE, filename = temp2[[2]]) # mean pd
  pe.rand.sd <- terra::stdev(pe.rand, na.rm = TRUE, overwrite = TRUE, filename = temp2[[3]]) # sd pd

  unlink(temp) # delete the archive that will not be used
  # unlink(temp2) # delete the archive that will not be used

  ## Calculating the standard effect size (SES)
  ses <- function(x, y, z){
    (x - y)/z
  }
  pe.ses <- ses(x = pe.obs, y = pe.rand.mean, z = pe.rand.sd)
  names(pe.ses) <- "SES"

  out <- c(pe.obs, pe.rand.mean, pe.rand.sd, pe.ses)
  names(out) <- c("PE Observed", "Mean", "SD", "SES" )
  return(out)
}

tes <- rast.pe.ses(pres.reord, branch.length, aleats = 10)
plot(tes)
