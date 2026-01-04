# flod

``` r
library(flod)
#> NB! `flod` is under heavy development and subject to sweeping changes
```

The flowSet is an S4 class and a convenient way of storing related .fcs
files in a single object. However, being an S4 object, it is naturally
restricted in its scope to the registered methods and, thus, foregoes
potentially relevant functions for EDA, transformation, and modelling
from other packages.
