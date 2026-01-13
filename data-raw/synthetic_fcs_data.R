# Generate Synthetic Flow Cytometry Data
# Creates realistic 4-color flow cytometry data for 3 samples

library(tidyverse)
library(duckdb)
library(connections)

set.seed(20260109)  # For reproducibility

# Helper function to generate correlated FSC-A and FSC-H for singlets
generate_singlets <- function(n, fsc_mean, fsc_sd) {
  fsc_a <- rlnorm(n, log(fsc_mean), fsc_sd)
  # FSC-H should be proportional to FSC-A for singlets with some noise
  fsc_h <- fsc_a * rnorm(n, mean = 0.85, sd = 0.05)
  list(fsc_a = fsc_a, fsc_h = fsc_h)
}

# Helper function to generate doublets (higher FSC-A for same FSC-H)
generate_doublets <- function(n, fsc_mean, fsc_sd) {
  fsc_h <- rlnorm(n, log(fsc_mean * 0.7), fsc_sd)
  fsc_a <- fsc_h * rnorm(n, mean = 1.5, sd = 0.1)  # Higher FSC-A
  list(fsc_a = fsc_a, fsc_h = fsc_h)
}

# Helper function to generate debris (low FSC/SSC)
generate_debris <- function(n) {
  fsc_a <- rlnorm(n, log(5000), 0.5)
  ssc_a <- rlnorm(n, log(3000), 0.5)
  fsc_h <- fsc_a * rnorm(n, mean = 0.7, sd = 0.2)
  list(fsc_a = fsc_a, ssc_a = ssc_a, fsc_h = fsc_h)
}

# Helper function to generate live cells (low FL1-A)
generate_live_cells <- function(n) {
  rlnorm(n, log(500), 0.4)
}

# Helper function to generate dead cells (high FL1-A)
generate_dead_cells <- function(n) {
  rlnorm(n, log(50000), 0.3)
}

# Helper function to generate marker populations
generate_marker <- function(n, pct_positive) {
  n_pos <- round(n * pct_positive)
  n_neg <- n - n_pos

  # Negative population (low fluorescence)
  neg <- rlnorm(n_neg, log(300), 0.3)
  # Positive population (high fluorescence)
  pos <- rlnorm(n_pos, log(15000), 0.4)

  sample(c(neg, pos))  # Shuffle
}

generate_sample <- function(
  sample_id,
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
  # Calculate population sizes
  n_debris <- round(n_events * debris_pct)
  n_dead <- round(n_events * dead_pct)
  n_live <- n_events - n_debris - n_dead
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
  population[outlier_idx] <- paste0(population[outlier_idx], "_outlier")
  
  # Create data frame
  df <- tibble(
    sample_id = sample_id,
    event_id = 1:n_events,
    `FSC-A` = pmax(0, fsc_a),
    `SSC-A` = pmax(0, ssc_a),
    `FSC-H` = pmax(0, fsc_h),
    `FL1-A` = pmax(0, fl1_a),
    `FL2-A` = pmax(0, fl2_a),
    `FL3-A` = pmax(0, fl3_a),
    population = population
  )
  
  # Gate 1: id_live (FSC-A vs FL1-A)
  # Calculate thresholds from live populations only
  live_populations <- c("live_singlet", "doublet")
  live_events <- df$population %in% live_populations
  
  fsc_threshold_live <- quantile(df$`FSC-A`[live_events], 0.99)
  fl1_threshold <- 5000  # Live/dead discriminator
  
  # Create id_live gate - excludes debris, dead cells, and outliers
  df$id_live <- !(df$population %in% c("debris", "debris_outlier", "dead", "dead_outlier")) &
                 df$`FL1-A` < fl1_threshold &
                 df$`FSC-A` <= fsc_threshold_live
  
# Gate 2: id_size (FSC-A vs SSC-A, applied to id_live population)
# Use Mahalanobis distance for elliptical gating (standard in flow cytometry)
if (sum(df$id_live) > 2) {  # Need at least 3 points for covariance
  # Extract FSC-A and SSC-A from id_live population
  live_data <- df[df$id_live, c("FSC-A", "SSC-A")]
  
  # Calculate center and covariance matrix
  center <- colMeans(live_data)
  cov_mat <- cov(live_data)
  
  # Calculate Mahalanobis distance for all events
  mahal_dist <- mahalanobis(
    df[, c("FSC-A", "SSC-A")],
    center = center,
    cov = cov_mat
  )
  
  # Use chi-square quantile for threshold (2 df for 2D data)
  # 0.99 quantile creates a gate that includes 99% of the live population
  threshold <- qchisq(0.99, df = 2)
  
  # Create id_size gate - elliptical boundary
  df$id_size <- df$id_live & (mahal_dist <= threshold)
} else {
  # Fallback if not enough live events
  df$id_size <- df$id_live
}
  
  df
}

# Define sample configurations as a tibble
sample_configs <- tibble(
  sample_id = c("Sample1", "Sample2", "Sample3"),
  n_events = c(20000, 20000, 20000),
  debris_pct = c(0.02, 0.025, 0.05),
  dead_pct = c(0.05, 0.06, 0.90),
  singlet_pct = c(0.95, 0.96, 0.85),
  fsc_mean = c(100000, 95000, 90000),
  fsc_sd = c(0.25, 0.25, 0.3),
  ssc_mean = c(60000, 58000, 55000),
  ssc_sd = c(0.35, 0.35, 0.4),
  fl2_positive_pct = c(runif(1, 0.25, 0.75), 0.05, 0.3),
  fl3_positive_pct = c(runif(1, 0.25, 0.75), 0.05, 0.3),
  outlier_pct = c(0.005, 0.005, 0.005)
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
