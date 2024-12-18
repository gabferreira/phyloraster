#' Classify Phylogenetic Endemism using rasters
#'
#' @description Use the results of rast.pe.ses() to identify centers of paleo-,
#' neo-, super-, and mixed- endemism following the CANAPE scheme of
#' Mishler et al., 2014.
#'
#' @param x SpatRaster. A SpatRaster object with the following layers in this
#' specific order:
#'  - pe.obs.p.upper : Upper p-value comparing the observed phylogenetic
#'  endemism and the randomized phylogenetic endemism values
#'  - pe.alt.obs.p.upper : Upper p-value comparing the alternate phylogenetic
#'  endemism and the randomized alternate phylogenetic endemism
#'  - rpe.obs.p.upper : Upper p-value comparing the relative phylogenetic
#'  endemism and the randomized relative phylogenetic endemism
#'  - rpe.obs.p.lower : Lower p-value comparing the relative phylogenetic
#'  endemism and the randomized relative phylogenetic endemism
#'
#' @return SpatRaster
#' @export
#' @references Mishler, B., Knerr, N., González-Orozco, C. et al. (2014)
#' Phylogenetic measures of biodiversity and neo- and paleo-endemism in
#' Australian Acacia. Nat Commun, 5: 4473. doi:10.1038/ncomms5473
#'
#' @keywords internal

.end.type <- function(x) {

  pe.obs.p.upper <- x[1]
  pe.alt.obs.p.upper <- x[2]
  rpe.obs.p.upper <- x[3]
  rpe.obs.p.lower <- x[4]

  # Classification of paleo, neo, super and mixed endemism
  if (is.na(pe.obs.p.upper) | is.na(pe.alt.obs.p.upper)
      | is.na(rpe.obs.p.upper) | is.na(rpe.obs.p.lower)) {
    return(NA)
  } else if ((pe.obs.p.upper > 0.95 | pe.alt.obs.p.upper > 0.95) &
             rpe.obs.p.upper > 0.975) {
    return(1)  # paleo
  } else if ((pe.obs.p.upper > 0.95 | pe.alt.obs.p.upper > 0.95) &
             rpe.obs.p.lower > 0.975) {
    return(2)  # neo
  } else if (pe.obs.p.upper > 0.99 & pe.alt.obs.p.upper > 0.99) {
    return(3)  # super
  } else if (pe.obs.p.upper > 0.95 | pe.alt.obs.p.upper > 0.95) {
    return(4)  # mixed
  } else {
    return(5)  # not significant
  }
}

#' Classify Phylogenetic Endemism using rasters
#'
#' @description Use the results of rast.pe.ses() to identify centers of paleo-,
#' neo-, super-, and mixed- endemism following the CANAPE scheme of
#' Mishler et al., 2014.
#'
#' @param pe.obs.p.upper SpatRaster. Upper p-value comparing the observed
#'  phylogenetic endemism and the randomized phylogenetic endemism values
#' @param pe.alt.obs.p.upper SpatRaster. Upper p-value comparing the alternate
#'  phylogenetic endemism and the randomized alternate phylogenetic endemism
#' @param rpe.obs.p.upper SpatRaster. Upper p-value comparing the relative
#'  phylogenetic endemism and the randomized relative phylogenetic endemism
#' @param rpe.obs.p.lower SpatRaster. Lower p-value comparing the relative
#'  phylogenetic endemism and the randomized relative phylogenetic endemism
#' @param filename character. Output filename
#' @param overwrite	 logical. If TRUE, filename is overwritten
#'
#' @return SpatRaster
#'
#' @seealso \code{\link{rast.pe.ses}},
#' \code{\link[SESraster]{SESraster}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(SESraster)
#' library(terra)
#' library(phyloraster)
#' x <- rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
#' tree <- ape::read.tree(system.file("extdata", "tree.nex",
#' package="phyloraster"))
#' data <- phylo.pres(x, tree)
#' ses <- rast.pe.ses(x = data$x, data$tree,
#'                 aleats = 5, metric = "all")
#' # CANAPE
#' canape <- canape.rast(ses$p.upper.PE, ses$p.upper.PE.alt,
#'                    ses$p.upper.RPE, ses$p.lower.RPE)
#' unique(canape)
#' plot(canape)
#'}
#'
canape.rast <- function(pe.obs.p.upper, pe.alt.obs.p.upper,
                        rpe.obs.p.upper, rpe.obs.p.lower, filename = NULL,
                        overwrite = FALSE){

  # if (!terra::is.lonlat(x)) {
  #   stop("Geographic coordinates are needed for the calculations.")
  # }

  # concatenate
  r <- c(pe.obs.p.upper, pe.alt.obs.p.upper,
         rpe.obs.p.upper, rpe.obs.p.lower)

  # Classify endemism
  endem.type.raster <- terra::app(r, .end.type)

  # SpatRaster Factor
  # 1:5 "paleo", "neo", "super", "mixed", "not significant"
  cls <- data.frame(id = c(1:5), cover = c("Paleo", "Neo",
                                           "Super", "Mixed",
                                           "Not significant"))
  levels(endem.type.raster) <- cls # as factor

  if(!is.null(filename)){ # to save the rasters when the output filename is
    # provided
    endem.type.raster <- terra::writeRaster(endem.type.raster, filename)
  }

  return(endem.type.raster)
}

