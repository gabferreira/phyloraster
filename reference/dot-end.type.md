# Classify Phylogenetic Endemism using rasters

Use the results of rast.pe.ses() to identify centers of paleo-, neo-,
super-, and mixed- endemism following the CANAPE scheme of Mishler et
al., 2014.

## Usage

``` r
.end.type(x)
```

## Arguments

- x:

  SpatRaster. A SpatRaster object with the following layers in this
  specific order:

  - pe.obs.p.upper : Upper p-value comparing the observed phylogenetic
    endemism and the randomized phylogenetic endemism values

  - pe.alt.obs.p.upper : Upper p-value comparing the alternate
    phylogenetic endemism and the randomized alternate phylogenetic
    endemism

  - rpe.obs.p.upper : Upper p-value comparing the relative phylogenetic
    endemism and the randomized relative phylogenetic endemism

  - rpe.obs.p.lower : Lower p-value comparing the relative phylogenetic
    endemism and the randomized relative phylogenetic endemism

## Value

SpatRaster

## References

Mishler, B., Knerr, N., González-Orozco, C. et al. (2014) Phylogenetic
measures of biodiversity and neo- and paleo-endemism in Australian
Acacia. Nat Commun, 5: 4473. doi:10.1038/ncomms5473
