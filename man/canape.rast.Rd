% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/canape.rast.R
\name{canape.rast}
\alias{canape.rast}
\title{Classify Phylogenetic Endemism using rasters}
\usage{
canape.rast(
  pe.obs.p.upper,
  pe.alt.obs.p.upper,
  rpe.obs.p.upper,
  rpe.obs.p.lower,
  filename = NULL,
  overwrite = FALSE
)
}
\arguments{
\item{pe.obs.p.upper}{SpatRaster. Upper p-value comparing the observed
phylogenetic endemism and the randomized phylogenetic endemism values}

\item{pe.alt.obs.p.upper}{SpatRaster. Upper p-value comparing the alternate
phylogenetic endemism and the randomized alternate phylogenetic endemism}

\item{rpe.obs.p.upper}{SpatRaster. Upper p-value comparing the relative
phylogenetic endemism and the randomized relative phylogenetic endemism}

\item{rpe.obs.p.lower}{SpatRaster. Lower p-value comparing the relative
phylogenetic endemism and the randomized relative phylogenetic endemism}

\item{filename}{character. Output filename}

\item{overwrite}{logical. If TRUE, filename is overwritten}
}
\value{
SpatRaster
}
\description{
Use the results of rast.pe.ses() to identify centers of paleo-,
neo-, super-, and mixed- endemism following the CANAPE scheme of
Mishler et al., 2014.
}
\examples{
\donttest{
library(SESraster)
library(terra)
library(phyloraster)
x <- rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
tree <- ape::read.tree(system.file("extdata", "tree.nex",
package="phyloraster"))
data <- phylo.pres(x, tree)
ses <- rast.pe.ses(x = data$x, data$tree,
                aleats = 5, metric = "all")
# CANAPE
canape <- canape.rast(ses$p.upper.PE, ses$p.upper.PE.alt,
                   ses$p.upper.RPE, ses$p.lower.RPE)
unique(canape)
plot(canape)
}

}
\seealso{
\code{\link{rast.pe.ses}},
\code{\link[SESraster]{SESraster}}
}
