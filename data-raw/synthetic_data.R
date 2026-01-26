# Generate Synthetic Flow Cytometry Data
library(tidyverse)
library(duckdb)
library(connections)

# Make experiments with synthetic data
set.seed(20260109)
experiment1 <- map(1:8, ~ synthetic_sample(n = 20000)) |>
  {
    \(x) tibble(exprs = x)
  }()

set.seed(20260109)
experiment2 <- map(1:5, ~ synthetic_sample(n = 20000)) |>
  {
    \(x) {
      x[[6]] <- synthetic_sample(live_pct = 0.1)
      x
    }
  }() |>
  {
    \(x) tibble(exprs = x)
  }()

# Common list for writing to output
common_list <- list(experiment1 = experiment1, experiment2 = experiment2)

# Create DuckDB database
dbdir <- "inst/extdata/synthetic_data.duckdb"
# con <- connection_open(duckdb(), dbdir = dbdir)
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = dbdir)

# Write data to DuckDB
iwalk(common_list, function(experiment, name) {
  duckdb::dbWriteTable(con, name, experiment, overwrite = TRUE)
})

# Close connection
# connections::connection_close(con)
DBI::dbDisconnect(con, shutdown = TRUE)
