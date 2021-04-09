
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CROCS (Changepoints for a Range of ComplexitieS)

CROCS is an extension of the CROPS [(Haynes et
al. 2017)](https://www.tandfonline.com/doi/abs/10.1080/10618600.2015.1116445)
and sequential search [(Hocking et
al. 2018)](https://arxiv.org/abs/1810.00117) algorithm which allows to
compute all optimal changepoint segmentations of data sequences for all
penalty values accross a peak range. This package implements the CROCS
algorithm as well as segmentation models for peak calling. They have
been described in our study [(Liehrmann et
al. 2020)](https://arxiv.org/abs/2012.06848) in a genomic context. On
top of the gfpop method [(Runge et
al. 2020)](https://arxiv.org/abs/2002.03646), it proposes a framework
to design your peak caller. The user can choose among several noise
assumptions `CROCS::lossFactory`, transformations
`CROCS::transformationFactory`, peak shape assumptions
`CROCS::graphFactory` and peak start/end post-processing rules
`CROCS::postProcessingRuleFactory` to build a peak caller
`CROCS::peakCallerFactory`. The `CROCS::CROCS` output is the first step
to compute an objective function that can be optimized in the supervised
learning procedure. We give two examples of its use in the package’s
[Vignette](vignettes/CROCS-vignette.html).

<!-- badges: start -->

<!-- badges: end -->

## Quick Start

We install the package from Github:

``` r
#devtools::install_github("aLiehrmann/CROCS")
library(CROCS)
```
