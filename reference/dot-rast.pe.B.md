# Calculate phylogenetic endemism for a raster

Calculate phylogenetic endemism using rasters as input and output.

## Usage

``` r
.rast.pe.B(
  x,
  inv.R,
  branch.length,
  branch.length.alt,
  metric = c("pe", "pe.alt", "rpe", "all")[1],
  filename = "",
  overwrite = TRUE,
  ...
)
```

## Arguments

- x:

  SpatRaster. A SpatRaster containing presence-absence data (0 or 1) for
  a set of species. The layers (species) will be sorted according to the
  tree order. See the phylo.pres function.

- inv.R:

  SpatRaster. Inverse of range size. See
  [`inv.range`](https://gabferreira.github.io/phyloraster/reference/inv.range.md)

- branch.length:

  numeric. A Named numeric vector of branch length for each species. See
  [`phylo.pres`](https://gabferreira.github.io/phyloraster/reference/phylo.pres.md)

- branch.length.alt:

  numeric. Branch length calculated by using an alternative phylogeny
  with non-zero branch lengths converted to a constant value (1) and
  rescaled so the sum of all branch lengths is 1.

- metric:

  character. Names of the biodiversity metrics to calculate. Available
  options are: "pe", "pe.alt", "rpe", or "all". See details.

- filename:

  character. Output filename

- overwrite:

  logical. If TRUE, filename is overwritten

- ...:

  additional arguments passed for terra::app

## Value

SpatRaster

## Details

Metrics available are:

- pe: Phylogenetic endemism (Rosauer et al., 2009)

- pe.alt: Alternate Phylogenetic endemism (Mishler et al., 2014)

- rpe: Relative Phylogenetic endemism (Mishler et al., 2014)

- all: Calculate all available metrics Alternate phylogenetic endemism
  (PE.alt, Mishler et al., 2014) is calculated using an alternate
  phylogeny with non-zero branch lengths converted to a constant value
  (here we use 1) and rescaled so the sum of all branch lengths is 1.
  Relative phylogenetic endemism (RPE, Mishler et al., 2014) is the
  ratio of phylogenetic endemism (PE, Rosauer et al., 2009) measured on
  the original tree versus PE measured on a alternate tree (PE.alt).

## References

Mishler, B. D., Knerr, N., González-Orozco, C. E., Thornhill, A. H.,
Laffan, S. W. and Miller, J. T. 2014. Phylogenetic measures of
biodiversity and neo- and paleo-endemism in Australian Acacia. – Nat.
Commun. 5: 4473.

Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C., &
Cook, L. G. (2009). Phylogenetic endemism: a new approach for
identifying geographical concentrations of evolutionary history.
Molecular ecology, 18(19), 4061-4072.

## Author

Gabriela Alves-Ferreira and Neander Heming
