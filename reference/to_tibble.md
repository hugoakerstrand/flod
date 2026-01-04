# Convert a flowSet into a tibble

`to_tibble()` creates a
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
from slots of a
[`flowCore::flowSet()`](https://rdrr.io/pkg/flowCore/man/flowSet-class.html);
an S4 class that stores .fcs file data.

## Usage

``` r
to_tibble(data, ...)
```

## Arguments

- data:

  A flowSet object.

- ...:

  What slots to extract? \<kw: default = "all"\> // One or more of
  `"all"`, `"sample_names"`, `"exprs"`, `"keywords"`, or `"meta_data"`.

## Value

A tibble containing the requested data components. The function will
error if `data` is not a flowSet or if invalid keywords are provided.

Additional keywords are ignored if combined with `all`.
