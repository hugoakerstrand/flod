# Generate Synthetic Flow Cytometry Data
# Creates realistic 4-color flow cytometry data for 3 samples

required_packages <- c("tidyverse", "duckdb")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cran.rstudio.com/", quiet = FALSE)
  } else {
    cat(sprintf("%s is already installed.\n", pkg))
  }
}

cat("\nAll dependencies installed successfully!\n")


library(tidyverse)
library(duckdb)

set.seed(42)  # For reproducibility

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

# Generate Sample 1: Normal sample with positive markers
generate_sample1 <- function(n_events = 20000) {
  # Population breakdown - mutually exclusive populations
  n_debris <- round(n_events * 0.02)  # 2% debris
  n_dead <- round(n_events * 0.05)  # 5% dead
  n_live <- n_events - n_debris - n_dead  # Remaining are live
  n_singlets <- round(n_live * 0.95)  # 95% of live cells are singlets
  n_doublets <- n_live - n_singlets

  # Initialize vectors
  fsc_a <- numeric(n_events)
  ssc_a <- numeric(n_events)
  fsc_h <- numeric(n_events)
  fl1_a <- numeric(n_events)
  fl2_a <- numeric(n_events)
  fl3_a <- numeric(n_events)

  current_idx <- 1

  # Generate debris (low FSC/SSC)
  if (n_debris > 0) {
    end_idx <- current_idx + n_debris - 1
    debris <- generate_debris(n_debris)
    fsc_a[current_idx:end_idx] <- debris$fsc_a
    ssc_a[current_idx:end_idx] <- debris$ssc_a
    fsc_h[current_idx:end_idx] <- debris$fsc_h
    fl1_a[current_idx:end_idx] <- generate_live_cells(n_debris)
    fl2_a[current_idx:end_idx] <- rlnorm(n_debris, log(200), 0.4)
    fl3_a[current_idx:end_idx] <- rlnorm(n_debris, log(200), 0.4)
    current_idx <- end_idx + 1
  }

  # Generate dead cells (high FL1-A, but still in normal FSC/SSC range)
  if (n_dead > 0) {
    end_idx <- current_idx + n_dead - 1
    dead_singlets <- generate_singlets(n_dead, 80000, 0.3)
    fsc_a[current_idx:end_idx] <- dead_singlets$fsc_a
    fsc_h[current_idx:end_idx] <- dead_singlets$fsc_h
    ssc_a[current_idx:end_idx] <- rlnorm(n_dead, log(50000), 0.4)
    fl1_a[current_idx:end_idx] <- generate_dead_cells(n_dead)
    fl2_a[current_idx:end_idx] <- rlnorm(n_dead, log(400), 0.3)
    fl3_a[current_idx:end_idx] <- rlnorm(n_dead, log(400), 0.3)
    current_idx <- end_idx + 1
  }

  # Generate live singlets (majority population)
  if (n_singlets > 0) {
    end_idx <- current_idx + n_singlets - 1
    singlets <- generate_singlets(n_singlets, 100000, 0.25)
    fsc_a[current_idx:end_idx] <- singlets$fsc_a
    fsc_h[current_idx:end_idx] <- singlets$fsc_h
    ssc_a[current_idx:end_idx] <- rlnorm(n_singlets, log(60000), 0.35)
    fl1_a[current_idx:end_idx] <- generate_live_cells(n_singlets)
    fl2_a[current_idx:end_idx] <- generate_marker(n_singlets, pct_positive = runif(1, 0.25, 0.75))
    fl3_a[current_idx:end_idx] <- generate_marker(n_singlets, pct_positive = runif(1, 0.25, 0.75))
    current_idx <- end_idx + 1
  }

  # Generate doublets
  if (n_doublets > 0) {
    end_idx <- current_idx + n_doublets - 1
    doublets <- generate_doublets(n_doublets, 100000, 0.25)
    fsc_a[current_idx:end_idx] <- doublets$fsc_a
    fsc_h[current_idx:end_idx] <- doublets$fsc_h
    ssc_a[current_idx:end_idx] <- rlnorm(n_doublets, log(70000), 0.35)
    fl1_a[current_idx:end_idx] <- generate_live_cells(n_doublets)
    fl2_a[current_idx:end_idx] <- generate_marker(n_doublets, pct_positive = 0.5)
    fl3_a[current_idx:end_idx] <- generate_marker(n_doublets, pct_positive = 0.5)
  }

  # Add some outliers (very high/low events)
  n_outliers <- round(n_events * 0.005)
  outlier_idx <- sample(1:n_events, n_outliers)
  fsc_a[outlier_idx] <- fsc_a[outlier_idx] * runif(n_outliers, 0.1, 3)
  ssc_a[outlier_idx] <- ssc_a[outlier_idx] * runif(n_outliers, 0.1, 3)

  tibble(
    sample_id = "Sample1",
    event_id = 1:n_events,
    `FSC-A` = pmax(0, fsc_a),
    `SSC-A` = pmax(0, ssc_a),
    `FSC-H` = pmax(0, fsc_h),
    `FL1-A` = pmax(0, fl1_a),
    `FL2-A` = pmax(0, fl2_a),
    `FL3-A` = pmax(0, fl3_a)
  )
}

# Generate Sample 2: Normal sample with background markers
generate_sample2 <- function(n_events = 20000) {
  # Similar structure to Sample1 but with slight variation - mutually exclusive populations
  n_debris <- round(n_events * 0.025)  # 2.5% debris
  n_dead <- round(n_events * 0.06)  # 6% dead
  n_live <- n_events - n_debris - n_dead  # Remaining are live
  n_singlets <- round(n_live * 0.96)  # 96% of live cells are singlets
  n_doublets <- n_live - n_singlets

  fsc_a <- numeric(n_events)
  ssc_a <- numeric(n_events)
  fsc_h <- numeric(n_events)
  fl1_a <- numeric(n_events)
  fl2_a <- numeric(n_events)
  fl3_a <- numeric(n_events)

  current_idx <- 1

  # Debris
  if (n_debris > 0) {
    end_idx <- current_idx + n_debris - 1
    debris <- generate_debris(n_debris)
    fsc_a[current_idx:end_idx] <- debris$fsc_a
    ssc_a[current_idx:end_idx] <- debris$ssc_a
    fsc_h[current_idx:end_idx] <- debris$fsc_h
    fl1_a[current_idx:end_idx] <- generate_live_cells(n_debris)
    fl2_a[current_idx:end_idx] <- rlnorm(n_debris, log(200), 0.4)
    fl3_a[current_idx:end_idx] <- rlnorm(n_debris, log(200), 0.4)
    current_idx <- end_idx + 1
  }

  # Dead cells
  if (n_dead > 0) {
    end_idx <- current_idx + n_dead - 1
    dead_singlets <- generate_singlets(n_dead, 85000, 0.3)
    fsc_a[current_idx:end_idx] <- dead_singlets$fsc_a
    fsc_h[current_idx:end_idx] <- dead_singlets$fsc_h
    ssc_a[current_idx:end_idx] <- rlnorm(n_dead, log(48000), 0.4)
    fl1_a[current_idx:end_idx] <- generate_dead_cells(n_dead)
    fl2_a[current_idx:end_idx] <- rlnorm(n_dead, log(400), 0.3)
    fl3_a[current_idx:end_idx] <- rlnorm(n_dead, log(400), 0.3)
    current_idx <- end_idx + 1
  }

  # Live singlets
  if (n_singlets > 0) {
    end_idx <- current_idx + n_singlets - 1
    singlets <- generate_singlets(n_singlets, 95000, 0.25)
    fsc_a[current_idx:end_idx] <- singlets$fsc_a
    fsc_h[current_idx:end_idx] <- singlets$fsc_h
    ssc_a[current_idx:end_idx] <- rlnorm(n_singlets, log(58000), 0.35)
    fl1_a[current_idx:end_idx] <- generate_live_cells(n_singlets)
    # Background/null signal for FL2-A and FL3-A (all low)
    fl2_a[current_idx:end_idx] <- rlnorm(n_singlets, log(350), 0.35)
    fl3_a[current_idx:end_idx] <- rlnorm(n_singlets, log(380), 0.35)
    current_idx <- end_idx + 1
  }

  # Doublets
  if (n_doublets > 0) {
    end_idx <- current_idx + n_doublets - 1
    doublets <- generate_doublets(n_doublets, 95000, 0.25)
    fsc_a[current_idx:end_idx] <- doublets$fsc_a
    fsc_h[current_idx:end_idx] <- doublets$fsc_h
    ssc_a[current_idx:end_idx] <- rlnorm(n_doublets, log(65000), 0.35)
    fl1_a[current_idx:end_idx] <- generate_live_cells(n_doublets)
    fl2_a[current_idx:end_idx] <- rlnorm(n_doublets, log(350), 0.35)
    fl3_a[current_idx:end_idx] <- rlnorm(n_doublets, log(380), 0.35)
  }

  # Outliers
  n_outliers <- round(n_events * 0.005)
  outlier_idx <- sample(1:n_events, n_outliers)
  fsc_a[outlier_idx] <- fsc_a[outlier_idx] * runif(n_outliers, 0.1, 3)
  ssc_a[outlier_idx] <- ssc_a[outlier_idx] * runif(n_outliers, 0.1, 3)

  tibble(
    sample_id = "Sample2",
    event_id = 1:n_events,
    `FSC-A` = pmax(0, fsc_a),
    `SSC-A` = pmax(0, ssc_a),
    `FSC-H` = pmax(0, fsc_h),
    `FL1-A` = pmax(0, fl1_a),
    `FL2-A` = pmax(0, fl2_a),
    `FL3-A` = pmax(0, fl3_a)
  )
}

# Generate Sample 3: Failed sample (95% dead)
generate_sample3 <- function(n_events = 20000) {
  # Mutually exclusive populations - failed sample has more debris and dead cells
  n_debris <- round(n_events * 0.05)  # 5% debris (more than normal)
  n_dead <- round(n_events * 0.90)  # 90% dead (most of the sample)
  n_live <- n_events - n_debris - n_dead  # Very few live cells remaining
  n_singlets <- round(n_live * 0.85)  # Fewer singlets in failed sample
  n_doublets <- n_live - n_singlets

  fsc_a <- numeric(n_events)
  ssc_a <- numeric(n_events)
  fsc_h <- numeric(n_events)
  fl1_a <- numeric(n_events)
  fl2_a <- numeric(n_events)
  fl3_a <- numeric(n_events)

  current_idx <- 1

  # Debris (more in failed sample)
  if (n_debris > 0) {
    end_idx <- current_idx + n_debris - 1
    debris <- generate_debris(n_debris)
    fsc_a[current_idx:end_idx] <- debris$fsc_a
    ssc_a[current_idx:end_idx] <- debris$ssc_a
    fsc_h[current_idx:end_idx] <- debris$fsc_h
    fl1_a[current_idx:end_idx] <- rlnorm(n_debris, log(20000), 0.6)  # Variable FL1
    fl2_a[current_idx:end_idx] <- rlnorm(n_debris, log(200), 0.4)
    fl3_a[current_idx:end_idx] <- rlnorm(n_debris, log(200), 0.4)
    current_idx <- end_idx + 1
  }

  # Dead cells (majority)
  if (n_dead > 0) {
    end_idx <- current_idx + n_dead - 1
    dead_singlets <- generate_singlets(n_dead, 70000, 0.35)
    fsc_a[current_idx:end_idx] <- dead_singlets$fsc_a
    fsc_h[current_idx:end_idx] <- dead_singlets$fsc_h
    ssc_a[current_idx:end_idx] <- rlnorm(n_dead, log(45000), 0.45)
    fl1_a[current_idx:end_idx] <- generate_dead_cells(n_dead)
    fl2_a[current_idx:end_idx] <- rlnorm(n_dead, log(400), 0.4)
    fl3_a[current_idx:end_idx] <- rlnorm(n_dead, log(400), 0.4)
    current_idx <- end_idx + 1
  }

  # Live singlets (very few)
  if (n_singlets > 0) {
    end_idx <- current_idx + n_singlets - 1
    singlets <- generate_singlets(n_singlets, 90000, 0.3)
    fsc_a[current_idx:end_idx] <- singlets$fsc_a
    fsc_h[current_idx:end_idx] <- singlets$fsc_h
    ssc_a[current_idx:end_idx] <- rlnorm(n_singlets, log(55000), 0.4)
    fl1_a[current_idx:end_idx] <- generate_live_cells(n_singlets)
    # Some marker signal even in failed sample
    fl2_a[current_idx:end_idx] <- generate_marker(n_singlets, pct_positive = 0.3)
    fl3_a[current_idx:end_idx] <- generate_marker(n_singlets, pct_positive = 0.3)
    current_idx <- end_idx + 1
  }

  # Doublets
  if (n_doublets > 0) {
    end_idx <- current_idx + n_doublets - 1
    doublets <- generate_doublets(n_doublets, 90000, 0.3)
    fsc_a[current_idx:end_idx] <- doublets$fsc_a
    fsc_h[current_idx:end_idx] <- doublets$fsc_h
    ssc_a[current_idx:end_idx] <- rlnorm(n_doublets, log(65000), 0.4)
    fl1_a[current_idx:end_idx] <- generate_live_cells(n_doublets)
    fl2_a[current_idx:end_idx] <- generate_marker(n_doublets, pct_positive = 0.3)
    fl3_a[current_idx:end_idx] <- generate_marker(n_doublets, pct_positive = 0.3)
  }

  tibble(
    sample_id = "Sample3",
    event_id = 1:n_events,
    `FSC-A` = pmax(0, fsc_a),
    `SSC-A` = pmax(0, ssc_a),
    `FSC-H` = pmax(0, fsc_h),
    `FL1-A` = pmax(0, fl1_a),
    `FL2-A` = pmax(0, fl2_a),
    `FL3-A` = pmax(0, fl3_a)
  )
}

# Generate all samples
cat("Generating Sample 1 (normal with positive markers)...\n")
sample1 <- generate_sample1()

cat("Generating Sample 2 (normal with background markers)...\n")
sample2 <- generate_sample2()

cat("Generating Sample 3 (failed sample, 95% dead)...\n")
sample3 <- generate_sample3()

# Combine all samples
all_samples <- bind_rows(sample1, sample2, sample3)

cat("\nData summary:\n")
print(all_samples %>%
  group_by(sample_id) %>%
  summarise(
    n_events = n(),
    mean_FSC_A = mean(`FSC-A`),
    mean_FL1_A = mean(`FL1-A`),
    mean_FL2_A = mean(`FL2-A`),
    mean_FL3_A = mean(`FL3-A`)
  ))

# Create DuckDB database
cat("\nCreating DuckDB database...\n")
con <- dbConnect(duckdb(), dbdir = "flow_cytometry_data.duckdb")

# Write data to DuckDB
dbWriteTable(con, "flow_data", all_samples, overwrite = TRUE)

# Verify the data
cat("\nVerifying data in DuckDB:\n")
result <- dbGetQuery(con, "SELECT sample_id, COUNT(*) as n_events FROM flow_data GROUP BY sample_id")
print(result)

# Show sample of data
cat("\nSample of data from each sample:\n")
sample_preview <- dbGetQuery(con, "
  SELECT * FROM flow_data
  WHERE event_id <= 5
  ORDER BY sample_id, event_id
")
print(sample_preview)

# Create summary statistics table
cat("\nCreating summary statistics table...\n")
dbExecute(con, "
  CREATE OR REPLACE TABLE flow_summary AS
  SELECT
    sample_id,
    COUNT(*) as total_events,
    AVG(\"FSC-A\") as mean_fsc_a,
    STDDEV(\"FSC-A\") as sd_fsc_a,
    AVG(\"SSC-A\") as mean_ssc_a,
    AVG(\"FL1-A\") as mean_fl1_a,
    AVG(\"FL2-A\") as mean_fl2_a,
    AVG(\"FL3-A\") as mean_fl3_a,
    -- Estimate live/dead based on FL1-A threshold
    SUM(CASE WHEN \"FL1-A\" < 5000 THEN 1 ELSE 0 END) as estimated_live_cells,
    SUM(CASE WHEN \"FL1-A\" >= 5000 THEN 1 ELSE 0 END) as estimated_dead_cells
  FROM flow_data
  GROUP BY sample_id
")

summary_stats <- dbGetQuery(con, "SELECT * FROM flow_summary")
cat("\nSummary statistics:\n")
print(summary_stats)

# List all tables in the database
cat("\nTables in database:\n")
print(dbListTables(con))

# Close connection
dbDisconnect(con, shutdown = TRUE)