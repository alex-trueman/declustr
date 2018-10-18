# declustr

`declustr` is a package for analysing the spatial patterns of point data in 
2D and 3D. For spatially clustered data there are various approaches
implimented for calculating declustering weights.

The current implimentation has functions for analysing and declustering 2D point
patterns using Voronoi tesselations (polygonal declustering) and functions from 
the `spatstat` point pattern analysis package (Baddeley, et al. 2015).

## Installation

`declustr` is not available on CRAN but can be installed from Github using the
`devtools` package.

``` r
# install.packages("devtools")
devtools::install_github("truemoid/declustr")
```

## Functions

`polydeclust2d` calcuates decluster weights for 2D point sample data using
Voronoi tesselations of the point data. Weights are calculated as the area of
influence of each sample divided by the total area of the domain and then
normalized so that the sum of all weights is 1. This normalization is the
default but can be turned off with the `normalize` argument.

The sample domain is defined by an input pixel mask, which is a regular grid of 
x and y coordinates that define the surface area of the domain. This pixel mask
can optionally be expanded by a set number of grid units. By default the 
expansion is 0; however, occassionally samples on the edge of the domain cause
errors in the tessellation. In these cases, setting the expansion to a small
number (e.g., 1 m) can resolve errors.

If the domain mask is extrapolated too too far beyond the limits of the sample
data the weights assigned to edge samples can be too high. To counter this issue
you could adjust the mask. Alternatively, the function argument `estdomain` when
set to TRUE will estimate a more compact domain using the Ripley-Rasson (1977)
method.

`plotable_tess` returns a dataframe of points defining the edges of the
tessellation constructed by `polydeclust2d`. These points can be used to
visualize the tessellations.

`point_spacing_2d` computes the avererage distance to the four nearest samples
in the mask domain. The number of nearest samples can be adjusted. The returned
object is a dataframe of regular grid points (as defined by the input point
pattern object) with column of distances.

## Referneces

Ripley, B.D. and Rasson, J.P. (1977) Finding the edge of a Poisson forest.
Journal of Applied Probability, 14, 483 – 491.

A. Baddeley, E. Rubak and R.Turner. Spatial Point Patterns: Methodology and
Applications with R. Chapman and Hall/CRC Press, 2015.
