# Generate Synthetic Flow Cytometry Data
library(tidyverse)
library(duckdb)
library(connections)

# Baseline experiments - i.e. not a lot of debris, dead cells, or duplicates
set.seed(20260109)
normal_experiments <- map(1:7, function(name, id) {
  map(1:8, ~ synthetic_sample(n = 15000)) |>
    {
      \(x) tibble(exprs = x)
    }()
}) |>
  setNames(c(paste0("experiment", 1:7)))

# Experiment 8 has sample 8 with low levels of viable cells (i.e. many dead cells)
set.seed(20260109)
experiment8 <- map(1:7, ~ synthetic_sample(n = 15000)) |>
  {
    \(x) {
      x[[8]] <- synthetic_sample(live_pct = 0.1)
      x
    }
  }() |>
  {
    \(x) tibble(exprs = x)
  }() |>
  list() |>
  setNames("experiment8")

# Common list for writing to output
common_list <- c(normal_experiments, experiment8)

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
