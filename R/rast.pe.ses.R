#' Calculate phylogenetic endemism (PE. Rosauer et al. 2009) for a raster
#'
#' Calculate phylogenetic endemism following Rosauer et al. (2009) using rasters as input and output.
#'
#' @param x SpatRaster. A presence-absence SpatRaster with the layers ordered according to the tree order
#' @param branch.length numeric. A Named numeric vector containing the branch length of each specie
#' @param filename character. Output filename.
#' @param ... additional arguments to be passed passed down from a calling function.
#' @return SpatRaster
# #' @export
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C., & Cook, L. G. (2009). Phylogenetic endemism: a new approach for identifying geographical concentrations of evolutionary history. Molecular ecology, 18(19), 4061-4072.
#' @examples
.rast.pe.B <- function(x, branch.length, filename = NULL, ...){

  {
    # criei a função .rast.pe.B que não retorna mensagem de erro quando
    # nomes sao diferentes, ao contrario da rast.pe. Assim é possivel rodar o
    # modelo nulo sem gerar erro pq os nomes sao diferentes no bl e x

    area.branch <- phylogrid::inv.range(x, branch.length)

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

#' Phylogenetic endemism standardized for specie richness
#'
#' The function calculates the PE corrected for species richness. In each null model run, richness is kept constant and branch lengths of each species are randomized. The function provides the mean and standard deviation of all null models and also calculates the standardized effect size (SES).
#'
#' @param x SpatRaster. A presence-absence SpatRaster with the layers ordered according to the branch.length order. See phylo.pres function.
#' @param branch.length numeric. A Named numeric vector containing the branch length of each specie. See phylo.pres function.
#' @param aleats numeric. A number indicating how many times the calculation should be repeated.
#' @param filename character. Output filename.
#' @param random character. A character indicating what type of randomization. Could be by tip, site, specie or full spat(site and specie).
#' @return SpatRaster
#' @author Gabriela Alves-Ferreira and Neander Heming
#' @export
#' @references Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C., & Cook, L. G. (2009). Phylogenetic endemism: a new approach for identifying geographical concentrations of evolutionary history. Molecular ecology, 18(19), 4061-4072.
#' @examples
#' \dontrun{
#' ras <- terra::rast(system.file("extdata", "rast.presab.tif", package="phylogrid"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phylogrid"))
#' data <- phylo.pres(ras, tree)
#' t <- rast.pe.ses(data$x, data$branch.length, aleats = 10, random = "fullspat")
#' plot(t)
#' }
#'
rast.pe.ses <- function(x, branch.length, aleats,
                        random = c("area.size", "site", "specie", "fullspat"),
                        filename = NULL){

  aleats <- aleats # number of null models
  temp <- vector("list", length = aleats) # to create a temporary vector with the raster number

  ## Null model (bootstrap structure)
  if(random == "tip"){
    pe.rand <- list() # to store the rasters in the loop
    bl.random <- branch.length
    for(i in 1:aleats){
      bl.random[] <- sample(branch.length, replace = T) # aleatorize the branch lenght
      ## check if the values are differents
      # branch.length == bl.random
      temp[[i]] <- paste0(tempfile(), i, ".tif") # directory to store the rasters

      pe.rand[[i]] <- .rast.pe.B(x, branch.length = bl.random,
                                 filename = temp[[i]])
    }


    pe.rand <- terra::rast(pe.rand) # to transform a list in raster
  } else if (random == "site"){

    pe.rand <- list() # to store the rasters in the loop

    for(i in 1:aleats){

      # temporary names to rasters
      temp[[i]] <- paste0(tempfile(), i, ".tif")

      ### embaralha por lyr - ordem dos sítios de cada espécie separada
      pres.site.null <- spat.rand(x, aleats = 1, random = "site")

      # calculate pe
      pe.rand[[i]] <- .rast.pe.B(pres.site.null, branch.length = branch.length,
                                 filename = temp[[i]])
    }

    pe.rand <- terra::rast(pe.rand) # to transform a list in raster

  } else if (random == "specie") {

    ### randomize by cells - species in each site
    pe.rand <- list() # to store the rasters in the loop

    for(i in 1:aleats){

      temp[[i]] <- paste0(tempfile(), i, ".tif") # temporary names to rasters
      sp.rand <- spat.rand(x, aleats = 1, random = "specie")

      pe.rand[[i]] <- .rast.pe.B(sp.rand, branch.length = branch.length,
                                 filename = temp[[i]])
    }

    pe.rand <- terra::rast(pe.rand) # to transform a list in raster

  } else if (random == "fullspat") {

    pe.rand <- list() # to store the rasters in the loop
    fr <- terra::freq(x)

    for(i in 1:aleats){

      temp[[i]] <- paste0(tempfile(), i, ".tif") # temporary names to rasters

      ### randomize sites and species
      pres.null <- terra::app(x, fun=.lyr.sample, fr=fr)

      pe.rand[[i]] <- .rast.pe.B(pres.null, branch.length = branch.length,
                                 filename = temp[[i]])

    }

    pe.rand <- terra::rast(pe.rand) # to transform a list in raster

  } else {
    stop("Choose a valid randomization method! The methods currently available are: 'tip','site', 'specie', 'fullspat'.")
  }

  ## PE observed
  pe.obs <- phylogrid::rast.pe(x, branch.length = branch.length, filename = filename)

  ## PD rand mean and PD rand SD
  pe.rand.mean <- terra::mean(pe.rand, na.rm = TRUE, filename = filename) # mean pd
  pe.rand.sd <- terra::stdev(pe.rand, na.rm = TRUE, filename = filename) # sd pd


  unlink(temp) # delete the archive that will not be used

  ## Calculating the standard effect size (SES)
  {
    ses <- function(x){
      (x[1] - x[2])/x[3]
    }
    pe.ses <- terra::app(c(pe.obs, pe.rand.mean, pe.rand.sd), fun = ses)
  }
  names(pe.ses) <- "SES"
  out <- c(pe.obs, pe.rand.mean, pe.rand.sd, pe.ses)
  names(out) <- c("PE Observed", "Mean", "SD", "SES" )

  if (!is.null(filename)){ # to save the rasters when the path is provide
    out <- terra::writeRaster(out, filename, ...)
  }

  return(out)
}
