# Generate Synthetic Flow Cytometry Data
# Creates realistic 4-color flow cytometry data for 3 samples

library(tidyverse)
library(duckdb)
library(duckplyr)
library(connections)

set.seed(20260109)  # For reproducibility

generate_sample <- function(
  n_events,
  debris_pct,
  dead_pct,
  singlet_pct,
  fsc_mean,
  fsc_sd,
  ssc_mean,
  ssc_sd,
  fl2_positive_pct,
  fl3_positive_pct,
  outlier_pct
) {
# Calculate population sizes - ensure they sum to n_events
  n_debris <- round(n_events * debris_pct)
  n_dead <- round(n_events * dead_pct)
  # Calculate n_live as remainder to guarantee exact sum
  n_live <- n_events - n_debris - n_dead
  
  # Handle edge case where percentages sum to > 1
  if (n_live < 0) {
    n_live <- runif(1, 100, 600)
    n_dead <- n_events - n_debris - n_live
  }
  
  n_singlets <- round(n_live * singlet_pct)
  n_doublets <- n_live - n_singlets
  
  # Initialize vectors
  fsc_a <- numeric(n_events)
  ssc_a <- numeric(n_events)
  fsc_h <- numeric(n_events)
  fl1_a <- numeric(n_events)
  fl2_a <- numeric(n_events)
  fl3_a <- numeric(n_events)
  population <- character(n_events)
  
  current_idx <- 1
  
  # Generate debris
  if (n_debris > 0) {
    end_idx <- current_idx + n_debris - 1
    debris <- generate_debris(n_debris)
    fsc_a[current_idx:end_idx] <- debris$fsc_a
    ssc_a[current_idx:end_idx] <- debris$ssc_a
    fsc_h[current_idx:end_idx] <- debris$fsc_h
    fl1_a[current_idx:end_idx] <- generate_live_cells(n_debris)
    fl2_a[current_idx:end_idx] <- rlnorm(n_debris, log(200), 0.4)
    fl3_a[current_idx:end_idx] <- rlnorm(n_debris, log(200), 0.4)
    population[current_idx:end_idx] <- "debris"
    current_idx <- end_idx + 1
  }
  
  # Generate dead cells
  if (n_dead > 0) {
    end_idx <- current_idx + n_dead - 1
    dead_singlets <- generate_singlets(n_dead, fsc_mean * 0.85, fsc_sd)
    fsc_a[current_idx:end_idx] <- dead_singlets$fsc_a
    fsc_h[current_idx:end_idx] <- dead_singlets$fsc_h
    ssc_a[current_idx:end_idx] <- rlnorm(n_dead, log(ssc_mean * 0.9), ssc_sd)
    fl1_a[current_idx:end_idx] <- generate_dead_cells(n_dead)
    fl2_a[current_idx:end_idx] <- rlnorm(n_dead, log(400), 0.3)
    fl3_a[current_idx:end_idx] <- rlnorm(n_dead, log(400), 0.3)
    population[current_idx:end_idx] <- "dead"
    current_idx <- end_idx + 1
  }
  
  # Generate live singlets
  if (n_singlets > 0) {
    end_idx <- current_idx + n_singlets - 1
    singlets <- generate_singlets(n_singlets, fsc_mean, fsc_sd)
    fsc_a[current_idx:end_idx] <- singlets$fsc_a
    fsc_h[current_idx:end_idx] <- singlets$fsc_h
    ssc_a[current_idx:end_idx] <- rlnorm(n_singlets, log(ssc_mean), ssc_sd)
    fl1_a[current_idx:end_idx] <- generate_live_cells(n_singlets)
    fl2_a[current_idx:end_idx] <- generate_marker(n_singlets, fl2_positive_pct)
    fl3_a[current_idx:end_idx] <- generate_marker(n_singlets, fl3_positive_pct)
    population[current_idx:end_idx] <- "live_singlet"
    current_idx <- end_idx + 1
  }
  
  # Generate doublets
  if (n_doublets > 0) {
    end_idx <- current_idx + n_doublets - 1
    doublets <- generate_doublets(n_doublets, fsc_mean, fsc_sd)
    fsc_a[current_idx:end_idx] <- doublets$fsc_a
    fsc_h[current_idx:end_idx] <- doublets$fsc_h
    ssc_a[current_idx:end_idx] <- rlnorm(n_doublets, log(ssc_mean * 1.1), ssc_sd)
    fl1_a[current_idx:end_idx] <- generate_live_cells(n_doublets)
    fl2_a[current_idx:end_idx] <- generate_marker(n_doublets, fl2_positive_pct)
    fl3_a[current_idx:end_idx] <- generate_marker(n_doublets, fl3_positive_pct)
    population[current_idx:end_idx] <- "doublet"
  }
  
  # Mark outliers
  n_outliers <- round(n_events * outlier_pct)
  outlier_idx <- sample(1:n_events, n_outliers)
  fsc_a[outlier_idx] <- fsc_a[outlier_idx] * runif(n_outliers, 0.1, 3)
  ssc_a[outlier_idx] <- ssc_a[outlier_idx] * runif(n_outliers, 0.1, 3)
  population[outlier_idx] <- population[outlier_idx]
  
  # Create data frame
  tibble(
    event_id = 1:n_events,
    `FSC-A` = pmax(0, fsc_a),
    `SSC-A` = pmax(0, ssc_a),
    `FSC-H` = pmax(0, fsc_h),
    `FL1-A` = pmax(0, fl1_a),
    `FL2-A` = pmax(0, fl2_a),
    `FL3-A` = pmax(0, fl3_a),
    population = population
  )
}

# Define sample configurations as a tibble
sample_configs <- tibble(
  n_events = rep(20000, 8),
  debris_pct = runif(8, 0.001, 0.5),
  dead_pct = c(rexp(7, 25), 0.9),
  singlet_pct = runif(8, 0.85, 0.99),
  fsc_mean = runif(8, 85000, 100000),
  fsc_sd = runif(8, 0.25, 0.4),
  ssc_mean = runif(8, 55000, 60000),
  ssc_sd = runif(8, 0.35, 0.4),
  fl2_positive_pct = c(runif(6, 0.25, 0.75), 0.05, 0.3),
  fl3_positive_pct = c(runif(6, 0.25, 0.75), 0.05, 0.3),
  outlier_pct = runif(8, 0, 0.005)
)

# Generate all samples using pmap
all_samples <- sample_configs |> 
  pmap(generate_sample) |> 
  {\(x) tibble(exprs = x)}()

# Create DuckDB database
dbdir <- "inst/extdata/flow_cytometry_data.duckdb"
con <- connection_open(duckdb(), dbdir = dbdir)

# Write data to DuckDB
duckdb::dbWriteTable(con, "experiment1", all_samples, overwrite = TRUE)

# Show sample of data
tbl(con, "experiment1") |> 
  as_duckdb_tibble() |> 
  {\(x) map(x[["exprs"]], head)}()

# Close connection
connections::connection_close(con)
