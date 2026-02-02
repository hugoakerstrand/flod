# flod

``` r
library(flod)
#> NB! `flod` is under heavy development and subject to sweeping changes
library(flowCore)
```

## Getting data into your working environment

Read in the .fcs files using
[`flowCore::read.flowSet()`](https://rdrr.io/pkg/flowCore/man/read.flowSet.html).
The flowSet is an S4 class and a convenient way of storing related .fcs
files in a single object. However, being an S4 object, it is naturally
restricted in its scope to the registered methods and, thus, foregoes
potentially relevant functions for EDA, transformation, and modelling
from other packages.

``` r
data(GvHD, package = "flowCore")
GvHD
#> A flowSet with 35 experiments.
#> 
#> An object of class 'AnnotatedDataFrame'
#>   rowNames: s5a01 s5a02 ... s10a07 (35 total)
#>   varLabels: Patient Visit ... name (5 total)
#>   varMetadata: labelDescription
#> 
#> column names(8): FSC-H SSC-H ... FL4-H Time
```

Turn it into a `tibble` using
[`flod::to_tibble()`](https://hugoakerstrand.github.io/flod/reference/to_tibble.md)

``` r
flod::to_tibble(GvHD)
#> # A tibble: 35 × 4
#>    sample_names exprs              keywords           meta_data$Patient $Visit
#>    <chr>        <list>             <list>             <fct>             <fct> 
#>  1 s5a01        <dbl [3,420 × 8]>  <named list [170]> 5                 1     
#>  2 s5a02        <dbl [3,405 × 8]>  <named list [170]> 5                 2     
#>  3 s5a03        <dbl [3,435 × 8]>  <named list [170]> 5                 3     
#>  4 s5a04        <dbl [8,550 × 8]>  <named list [170]> 5                 4     
#>  5 s5a05        <dbl [10,410 × 8]> <named list [170]> 5                 5     
#>  6 s5a06        <dbl [3,750 × 8]>  <named list [170]> 5                 6     
#>  7 s5a07        <dbl [13,979 × 8]> <named list [170]> 5                 7     
#>  8 s6a01        <dbl [2,205 × 8]>  <named list [167]> 6                 1     
#>  9 s6a02        <dbl [13,350 × 8]> <named list [167]> 6                 2     
#> 10 s6a03        <dbl [11,610 × 8]> <named list [167]> 6                 3     
#> # ℹ 25 more rows
#> # ℹ 3 more variables: meta_data$Days <int>, $Grade <int>, $name <chr>
```

The function defaults to all slots of the flowSet, or accepts a user
specified selection:

``` r
to_tibble(GvHD, c("exprs", "keywords"))
#> # A tibble: 35 × 2
#>    exprs              keywords          
#>    <list>             <list>            
#>  1 <dbl [3,420 × 8]>  <named list [170]>
#>  2 <dbl [3,405 × 8]>  <named list [170]>
#>  3 <dbl [3,435 × 8]>  <named list [170]>
#>  4 <dbl [8,550 × 8]>  <named list [170]>
#>  5 <dbl [10,410 × 8]> <named list [170]>
#>  6 <dbl [3,750 × 8]>  <named list [170]>
#>  7 <dbl [13,979 × 8]> <named list [170]>
#>  8 <dbl [2,205 × 8]>  <named list [167]>
#>  9 <dbl [13,350 × 8]> <named list [167]>
#> 10 <dbl [11,610 × 8]> <named list [167]>
#> # ℹ 25 more rows
```

See
[`?to_tibble`](https://hugoakerstrand.github.io/flod/reference/to_tibble.md)
for a list of accepted keywords.
